"""Classical waiting-time model for a deterministic repeater chain.

The definitions in this module follow Eqs. (C.ZN) and (C.WN) of the
manuscript.  Keeping this classical stochastic layer separate from the quantum
state calculation makes both parts independently testable.
"""

from __future__ import annotations

import math
from itertools import combinations, permutations
from typing import Iterable

import numpy as np
from numpy.typing import ArrayLike, NDArray


def completion_rounds(success_rounds: ArrayLike) -> NDArray[np.int64] | np.int64:
    """Return ``Z_N = max_i X_i`` along the last array axis."""

    rounds = np.asarray(success_rounds, dtype=np.int64)
    if rounds.ndim == 0 or rounds.shape[-1] < 1:
        raise ValueError("success_rounds must contain at least one link")
    if np.any(rounds < 1):
        raise ValueError("geometric success rounds must be positive")
    return np.max(rounds, axis=-1)


def accumulated_storage_rounds(
    success_rounds: ArrayLike,
) -> NDArray[np.int64] | np.int64:
    """Return the accumulated storage rounds from manuscript Eq. (C.WN)."""

    rounds = np.asarray(success_rounds, dtype=np.int64)
    if rounds.ndim == 0 or rounds.shape[-1] < 1:
        raise ValueError("success_rounds must contain at least one link")
    if np.any(rounds < 1):
        raise ValueError("geometric success rounds must be positive")

    final_round = np.max(rounds, axis=-1)
    end_memories = (final_round - rounds[..., 0]) + (final_round - rounds[..., -1])
    repeater_memories = np.abs(np.diff(rounds, axis=-1)).sum(axis=-1)
    return end_memories + repeater_memories


def expected_completion_round(success_probabilities: Iterable[float]) -> float:
    """Exact ``E[max_i X_i]`` for independent geometric link times.

    Each geometric variable has support ``{1, 2, ...}``.  Inclusion--exclusion
    gives a short exact expression and also covers heterogeneous links.
    """

    probabilities = np.asarray(tuple(success_probabilities), dtype=float)
    if probabilities.ndim != 1 or probabilities.size < 1:
        raise ValueError("at least one success probability is required")
    if np.any((probabilities <= 0.0) | (probabilities > 1.0)):
        raise ValueError("success probabilities must lie in (0, 1]")

    failure_probabilities = 1.0 - probabilities
    expectation = 0.0
    indices = range(probabilities.size)
    for subset_size in range(1, probabilities.size + 1):
        sign = 1.0 if subset_size % 2 else -1.0
        for subset in combinations(indices, subset_size):
            joint_failure = float(np.prod(failure_probabilities[list(subset)]))
            expectation += sign / (1.0 - joint_failure)
    return expectation


def expected_continuum_storage(repeater_nodes: int) -> float:
    """Return the exact mean of the continuum-scaled storage variable.

    For ``m = repeater_nodes + 1`` independent unit-rate exponential
    link times, ``E[max_i Y_i] = H_m`` and
    ``E[|Y_i-Y_{i-1}|] = 1``.  Manuscript Eq. (C.WN) therefore gives
    ``E[Wbar_N] = 2 H_(N+1) + N - 2``.
    """

    if repeater_nodes < 1:
        raise ValueError("repeater_nodes must be positive")
    harmonic = math.fsum(1.0 / k for k in range(1, repeater_nodes + 2))
    return 2.0 * harmonic + repeater_nodes - 2.0


def continuum_storage_laplace(repeater_nodes: int, s: float) -> float:
    """Evaluate the arbitrary-N Laplace formula of Appendix D.

    The exact finite permutation sum scales factorially and is intended for
    analytical checks at modest node counts. Monte Carlo sampling is more
    efficient for long chains.
    """

    if repeater_nodes < 1:
        raise ValueError("repeater_nodes must be positive")
    if s < 0:
        raise ValueError("s must be nonnegative")

    node_count = repeater_nodes
    link_count = node_count + 1
    total = 0.0
    path_edges = [(index - 1, index) for index in range(1, node_count + 1)]
    for ordering in permutations(range(link_count)):
        maximum_index = ordering[-1]
        edges = [(maximum_index, 0), (maximum_index, node_count), *path_edges]
        lower_set: set[int] = set()
        conditional_transform = 1.0
        for spacing_index in range(1, node_count + 1):
            lower_set.add(ordering[spacing_index - 1])
            cut_edges = sum(
                (left in lower_set) != (right in lower_set)
                for left, right in edges
            )
            spacing_rate = link_count - spacing_index
            conditional_transform *= spacing_rate / (
                spacing_rate + s * cut_edges
            )
        total += conditional_transform
    return total / math.factorial(link_count)


def sample_waiting_times(
    success_probabilities: float | Iterable[float],
    *,
    link_count: int | None = None,
    trials: int = 100_000,
    seed: int | None = None,
) -> tuple[NDArray[np.int64], NDArray[np.int64]]:
    """Sample exact no-cutoff ``Z_N`` and ``W_N`` values.

    A scalar probability requires ``link_count``.  An iterable specifies one
    probability per link and must not be combined with ``link_count``.
    """

    if trials < 1:
        raise ValueError("trials must be positive")

    if np.isscalar(success_probabilities):
        if link_count is None or link_count < 1:
            raise ValueError("link_count must be positive for a scalar probability")
        probabilities = np.full(link_count, float(success_probabilities))
    else:
        if link_count is not None:
            raise ValueError("link_count is inferred from a probability sequence")
        probabilities = np.asarray(tuple(success_probabilities), dtype=float)

    if probabilities.ndim != 1 or probabilities.size < 1:
        raise ValueError("at least one success probability is required")
    if np.any((probabilities <= 0.0) | (probabilities > 1.0)):
        raise ValueError("success probabilities must lie in (0, 1]")

    rng = np.random.default_rng(seed)
    rounds = rng.geometric(probabilities, size=(trials, probabilities.size))
    return completion_rounds(rounds), accumulated_storage_rounds(rounds)


def sample_continuum_waiting_times(
    repeater_nodes: int,
    *,
    trials: int = 100_000,
    seed: int | None = None,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Sample the strict small-success-probability continuum variables.

    The scaled link-generation times are independent unit-rate exponential
    variables. The returned arrays are ``Zbar_N=max_i(Y_i)`` and the
    continuum storage functional in manuscript Eq. (D.4).
    """

    if repeater_nodes < 0:
        raise ValueError("repeater_nodes cannot be negative")
    if trials < 1:
        raise ValueError("trials must be positive")

    rng = np.random.default_rng(seed)
    times = rng.exponential(
        scale=1.0,
        size=(trials, repeater_nodes + 1),
    )
    completion = np.max(times, axis=-1)
    end_memories = (
        completion - times[..., 0]
        + completion - times[..., -1]
    )
    repeater_memories = np.abs(np.diff(times, axis=-1)).sum(axis=-1)
    return completion, end_memories + repeater_memories


def coherence_factor(storage_rounds: ArrayLike, repetition_rate: float, coherence_time: float) -> float:
    """Estimate ``E[exp(-W_N/(f T_coh))]`` from storage-round samples."""

    if repetition_rate <= 0.0 or coherence_time <= 0.0:
        raise ValueError("repetition_rate and coherence_time must be positive")
    storage = np.asarray(storage_rounds, dtype=float)
    return float(np.exp(-storage / (repetition_rate * coherence_time)).mean())
