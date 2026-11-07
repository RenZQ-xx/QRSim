"""Validate the continuum storage-time approximation at manuscript points.

The exact finite-probability geometric and limiting exponential waiting times
are coupled through the same exponential variates. This common-random-number
construction makes the small systematic difference much easier to resolve
than two independent Monte Carlo calculations.
"""

from __future__ import annotations

import argparse
import csv
import math
from pathlib import Path

import numpy as np

from QRSim.analytical import homogeneous_chain
from QRSim.waiting_time import accumulated_storage_rounds


from reproduce_figures import MANUSCRIPT_PARAMETERS

PARAMETERS = dict(MANUSCRIPT_PARAMETERS)

# Maximum-fidelity points for N=4,5,7,8 and the N=6 point quoted in the
# abstract. All points belong to the 1000-km scan in manuscript Fig. 7.
OPERATING_POINTS = (
    (4, 1000.0, 0.050),
    (5, 1000.0, 0.034),
    (6, 1000.0, 0.036),
    (7, 1000.0, 0.021),
    (8, 1000.0, 0.018),
)


def storage_function(times: np.ndarray) -> np.ndarray:
    """Evaluate the storage functional for continuous link times."""

    completion = np.max(times, axis=-1)
    return (
        completion - times[..., 0]
        + completion - times[..., -1]
        + np.abs(np.diff(times, axis=-1)).sum(axis=-1)
    )


def validate_point(
    repeater_nodes: int,
    total_distance_km: float,
    excitation_probability: float,
    *,
    trials: int,
    seed: int,
) -> dict[str, float | int]:
    """Return a paired finite-p/continuum comparison for one operating point."""

    result = homogeneous_chain(
        repeater_nodes=repeater_nodes,
        total_distance_km=total_distance_km,
        excitation_probability=excitation_probability,
        **PARAMETERS,
    )
    probability = result.elementary_link.heralding_probability
    rng = np.random.default_rng(seed)
    exponential_times = rng.exponential(
        size=(trials, repeater_nodes + 1),
    )

    # If Y is unit exponential and lambda=-log(1-p), then
    # floor(Y/lambda)+1 is exactly geometric on {1,2,...}.
    geometric_rounds = (
        np.floor(exponential_times / -math.log1p(-probability)).astype(np.int64)
        + 1
    )
    finite_storage = accumulated_storage_rounds(geometric_rounds)
    continuum_storage = storage_function(exponential_times)

    finite_samples = np.exp(
        -finite_storage / (result.repetition_rate * PARAMETERS["coherence_time"])
    )
    continuum_samples = np.exp(
        -continuum_storage
        / (
            probability
            * result.repetition_rate
            * PARAMETERS["coherence_time"]
        )
    )
    paired_difference = continuum_samples - finite_samples
    kappa_finite = float(finite_samples.mean())
    kappa_continuum = float(continuum_samples.mean())
    difference = float(paired_difference.mean())
    difference_se = float(paired_difference.std(ddof=1) / math.sqrt(trials))

    coherence_prefactor = (
        PARAMETERS["swapping_quality"] ** repeater_nodes
        * (
            result.elementary_link.signal
            - result.elementary_link.error
        )
        ** (repeater_nodes + 1)
        / 2
    )
    return {
        "repeater_nodes": repeater_nodes,
        "total_distance_km": total_distance_km,
        "excitation_probability": excitation_probability,
        "heralding_probability": probability,
        "trials": trials,
        "seed": seed,
        "kappa_finite_p": kappa_finite,
        "kappa_continuum": kappa_continuum,
        "continuum_minus_finite": difference,
        "paired_standard_error": difference_se,
        "fidelity_shift": coherence_prefactor * difference,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--trials", type=int, default=1_000_000)
    parser.add_argument("--seed", type=int, default=20_260_728)
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("results/validation/continuum_validation.csv"),
    )
    args = parser.parse_args()
    if args.trials < 1:
        parser.error("--trials must be positive")

    rows = [
        validate_point(
            nodes,
            distance,
            excitation_probability,
            trials=args.trials,
            seed=args.seed + nodes,
        )
        for nodes, distance, excitation_probability in OPERATING_POINTS
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

    for row in rows:
        print(
            "N={repeater_nodes}: p_hr={heralding_probability:.6g}, "
            "Delta kappa={continuum_minus_finite:+.3e} +/- "
            "{paired_standard_error:.1e}, Delta F={fidelity_shift:+.3e}".format(
                **row
            )
        )
    print(args.output.resolve())


if __name__ == "__main__":
    main()
