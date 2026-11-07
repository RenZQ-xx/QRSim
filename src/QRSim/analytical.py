"""Closed-form model used by the current QRSim manuscript.

This module implements the symmetric elementary-link coefficients in Appendix
A and the homogeneous-chain fidelity in Eq. (FidelityN).  It deliberately does
not contain plotting code.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np

from QRSim.MonteCarloSim import MonteCarloSim
from QRSim.waiting_time import expected_completion_round


@dataclass(frozen=True)
class BellDiagonalWeights:
    """Bell-basis probabilities ordered around the target ``|psi+>``.

    The four entries are the Pauli-error probabilities relative to the target:
    identity (``psi+``), phase flip (``psi-``), bit flip (``phi+``), and a
    combined bit/phase flip (``phi-``).
    """

    psi_plus: float
    psi_minus: float
    phi_plus: float
    phi_minus: float

    @property
    def probabilities(self) -> np.ndarray:
        return np.array(
            [self.psi_plus, self.psi_minus, self.phi_plus, self.phi_minus],
            dtype=float,
        )

    @property
    def parity_contrast(self) -> float:
        """Odd- versus even-excitation Bell-sector contrast."""

        return self.psi_plus + self.psi_minus - self.phi_plus - self.phi_minus

    @property
    def phase_contrast_plus(self) -> float:
        return self.psi_plus - self.psi_minus + self.phi_plus - self.phi_minus

    @property
    def phase_contrast_minus(self) -> float:
        return self.psi_plus - self.psi_minus - self.phi_plus + self.phi_minus

    def density_matrix(self) -> np.ndarray:
        """Return the state in the ``ee, eg, ge, gg`` computational basis."""

        return np.array(
            [
                [
                    (self.phi_plus + self.phi_minus) / 2,
                    0,
                    0,
                    (self.phi_plus - self.phi_minus) / 2,
                ],
                [
                    0,
                    (self.psi_plus + self.psi_minus) / 2,
                    (self.psi_plus - self.psi_minus) / 2,
                    0,
                ],
                [
                    0,
                    (self.psi_plus - self.psi_minus) / 2,
                    (self.psi_plus + self.psi_minus) / 2,
                    0,
                ],
                [
                    (self.phi_plus - self.phi_minus) / 2,
                    0,
                    0,
                    (self.phi_plus + self.phi_minus) / 2,
                ],
            ],
            dtype=float,
        )


@dataclass(frozen=True)
class ElementaryLinkResult:
    signal_weight: float
    error_weight: float
    noise_ee: float
    noise_eg: float
    noise_ge: float
    noise_gg: float
    heralding_probability: float
    fidelity: float

    @property
    def signal(self) -> float:
        return 2 * self.signal_weight / self.heralding_probability

    @property
    def error(self) -> float:
        return 2 * self.error_weight / self.heralding_probability

    @property
    def population_ee(self) -> float:
        return 2 * self.noise_ee / self.heralding_probability

    @property
    def population_eg_noise(self) -> float:
        return 2 * self.noise_eg / self.heralding_probability

    @property
    def population_gg(self) -> float:
        return 2 * self.noise_gg / self.heralding_probability

    @property
    def delta(self) -> float:
        return 1 - 2 * (self.population_ee + self.population_gg)

    @property
    def population_single_excitation(self) -> float:
        return self.population_eg_noise + (self.signal + self.error) / 2

    @property
    def coherence(self) -> float:
        return (self.signal - self.error) / 2

    def density_matrix(self) -> np.ndarray:
        """Return the normalized component-derived elementary X state."""

        b = self.population_single_excitation
        h = self.coherence
        return np.array(
            [
                [self.population_ee, 0, 0, 0],
                [0, b, h, 0],
                [0, h, b, 0],
                [0, 0, 0, self.population_gg],
            ],
            dtype=float,
        )

    def strict_pauli_twirl(self) -> BellDiagonalWeights:
        """Return the exact Bell-diagonal Pauli twirl of the X state."""

        even_weight = (self.population_ee + self.population_gg) / 2
        return BellDiagonalWeights(
            psi_plus=self.fidelity,
            psi_minus=self.population_single_excitation - self.coherence,
            phi_plus=even_weight,
            phi_minus=even_weight,
        )


@dataclass(frozen=True)
class ChainResult:
    fidelity: float
    rate: float
    coherence_factor: float
    elementary_link: ElementaryLinkResult
    repetition_rate: float


def fidelity_matched_bell_state(
    elementary_fidelity: float,
    model: str,
) -> BellDiagonalWeights:
    """Construct a controlled Bell-diagonal reduction at fixed fidelity.

    These states change only the assumed distribution of the infidelity.  The
    caller retains the component model's heralding probability, waiting-time
    statistics, memory factor, and local swapping quality.
    """

    if not 0 <= elementary_fidelity <= 1:
        raise ValueError("elementary_fidelity must lie in [0, 1]")
    infidelity = 1 - elementary_fidelity
    if model == "werner":
        other = infidelity / 3
        return BellDiagonalWeights(elementary_fidelity, other, other, other)
    if model == "phase_only":
        return BellDiagonalWeights(elementary_fidelity, infidelity, 0.0, 0.0)
    if model == "bit_only":
        return BellDiagonalWeights(elementary_fidelity, 0.0, infidelity, 0.0)
    raise ValueError("model must be 'werner', 'phase_only', or 'bit_only'")


def bell_diagonal_chain_fidelity(
    state: BellDiagonalWeights,
    *,
    repeater_nodes: int,
    coherence_factor: float,
    swapping_quality: float,
) -> float:
    """Propagate a Bell-diagonal link through the homogeneous chain.

    Bell-state errors compose by convolution on the two-bit Pauli group.  The
    three nontrivial Fourier components are the contrasts below.  Pure memory
    dephasing leaves the parity contrast unchanged and multiplies both
    phase-sensitive contrasts by the same accumulated-storage factor.
    """

    if repeater_nodes < 0:
        raise ValueError("repeater_nodes cannot be negative")
    if not 0 <= coherence_factor <= 1:
        raise ValueError("coherence_factor must lie in [0, 1]")
    if not 0 <= swapping_quality <= 1:
        raise ValueError("swapping_quality must lie in [0, 1]")
    probabilities = state.probabilities
    if not np.isclose(probabilities.sum(), 1.0, atol=1e-12):
        raise ValueError("Bell-state probabilities must sum to one")
    if np.any(probabilities < -1e-12):
        raise ValueError("Bell-state probabilities cannot be negative")

    link_count = repeater_nodes + 1
    quality = swapping_quality**repeater_nodes
    return float(
        (
            1
            + quality
            * (
                state.parity_contrast**link_count
                + coherence_factor
                * (
                    state.phase_contrast_plus**link_count
                    + state.phase_contrast_minus**link_count
                )
            )
        )
        / 4
    )


def reduced_model_fidelity(
    elementary_fidelity: float,
    model: str,
    *,
    repeater_nodes: int,
    coherence_factor: float,
    swapping_quality: float,
) -> float:
    """Return a fidelity-matched reduced model's end-to-end fidelity."""

    return bell_diagonal_chain_fidelity(
        fidelity_matched_bell_state(elementary_fidelity, model),
        repeater_nodes=repeater_nodes,
        coherence_factor=coherence_factor,
        swapping_quality=swapping_quality,
    )


def no_click_probability(photon_number: int, detector_efficiency: float, dark_count_probability: float) -> float:
    if photon_number < 0:
        raise ValueError("photon_number cannot be negative")
    if not 0 <= detector_efficiency <= 1 or not 0 <= dark_count_probability <= 1:
        raise ValueError("detector parameters must lie in [0, 1]")
    return (1 - dark_count_probability) * (
        1 - detector_efficiency * (1 - dark_count_probability)
    ) ** photon_number


def symmetric_elementary_link(
    *,
    excitation_probability: float,
    transmission: float,
    detector_efficiency: float,
    dark_count_probability: float,
    visibility: float,
    phase_coherence: float = 1.0,
) -> ElementaryLinkResult:
    """Evaluate the manuscript's symmetric single-click elementary link."""

    for name, value in {
        "excitation_probability": excitation_probability,
        "transmission": transmission,
        "detector_efficiency": detector_efficiency,
        "dark_count_probability": dark_count_probability,
        "visibility": visibility,
        "phase_coherence": phase_coherence,
    }.items():
        if not 0 <= value <= 1:
            raise ValueError(f"{name} must lie in [0, 1]")

    pe = excitation_probability
    eta = transmission
    c = visibility
    no = [no_click_probability(n, detector_efficiency, dark_count_probability) for n in range(3)]
    click = [1 - value for value in no]

    signal = click[1] * no[0] * (1 - pe) * pe * eta * np.sqrt(c)
    error = no[1] * click[0] * (1 - pe) * pe * eta * np.sqrt(c)
    n_ee = no[0] * click[0] * (1 - pe) ** 2
    n_eg = (
        no[0] * click[0] * (1 - eta) * (1 - pe) * pe
        + (no[0] * click[1] + no[1] * click[0])
        * (1 - np.sqrt(c)) / 2 * eta * (1 - pe) * pe
    )
    n_gg = (
        no[0] * click[0] * pe**2 * (1 - eta) ** 2
        + (no[0] * click[1] + no[1] * click[0]) * eta * (1 - eta) * pe**2
        + (no[0] * click[2] + no[2] * click[0]) * (1 + c) / 4 * eta**2 * pe**2
        + no[1] * click[1] * (1 - c) / 2 * eta**2 * pe**2
    )
    heralding = 2 * (signal + error + n_ee + 2 * n_eg + n_gg)
    # Average a residual relative optical phase, after removing its mean.
    # This redistributes the two Bell weights but leaves all populations,
    # two-photon event probabilities, and the heralding probability intact.
    original_signal, original_error = signal, error
    signal = ((original_signal + original_error)
              + phase_coherence * (original_signal - original_error)) / 2
    error = ((original_signal + original_error)
             - phase_coherence * (original_signal - original_error)) / 2
    fidelity = (signal + n_eg) / (heralding / 2) if heralding else 0.0
    return ElementaryLinkResult(
        signal_weight=float(signal),
        error_weight=float(error),
        noise_ee=float(n_ee),
        noise_eg=float(n_eg),
        noise_ge=float(n_eg),
        noise_gg=float(n_gg),
        heralding_probability=float(heralding),
        fidelity=float(fidelity),
    )


def fiber_transmission(
    arm_distance_km: float,
    *,
    transmission_without_fiber: float,
    attenuation_db_per_km: float,
) -> float:
    if arm_distance_km < 0 or attenuation_db_per_km < 0:
        raise ValueError("distance and attenuation cannot be negative")
    return transmission_without_fiber * 10 ** (-attenuation_db_per_km * arm_distance_km / 10)


def repetition_rate(
    elementary_link_distance_km: float,
    *,
    light_speed_km_per_s: float = 200_000.0,
    local_cycle_time_s: float = 50e-6,
) -> float:
    """Synchronized round rate for a midpoint-heralded elementary link."""

    if elementary_link_distance_km < 0 or light_speed_km_per_s <= 0 or local_cycle_time_s < 0:
        raise ValueError("invalid timing parameter")
    return 1 / (elementary_link_distance_km / light_speed_km_per_s + local_cycle_time_s)


def continuum_coherence_factor(
    repeater_nodes: int,
    *,
    heralding_probability: float,
    attempt_rate: float,
    coherence_time: float,
) -> float:
    """Evaluate the manuscript's small-heralding-probability ``kappa_N``.

    For ``N=1,2,3`` this uses the explicit densities in Eqs. (C.Pi1--C.Pi3).
    Larger chains use deterministically regenerated Monte Carlo samples of the
    strict continuum limiting variable.
    """

    if repeater_nodes < 0:
        raise ValueError("repeater_nodes cannot be negative")
    if coherence_time == float("inf") or repeater_nodes == 0:
        return 1.0
    if heralding_probability <= 0 or attempt_rate <= 0 or coherence_time <= 0:
        raise ValueError("probability, rate, and coherence time must be positive")

    chi = heralding_probability * attempt_rate * coherence_time
    inverse_chi = 1 / chi
    if repeater_nodes == 1:
        return chi / (chi + 2)
    if repeater_nodes == 2:
        return (
            (1 / 9) / (inverse_chi + 1 / 4)
            + (2 / 3) / (inverse_chi + 1 / 2)
            - (7 / 9) / (inverse_chi + 1)
        )
    if repeater_nodes == 3:
        return (
            (4 / 15) / (inverse_chi + 1 / 4)
            + (3 / 16) / (inverse_chi + 1 / 2)
            + (1 / 16) / (inverse_chi + 1 / 2) ** 2
            - (7 / 6) / (inverse_chi + 1)
            + (57 / 80) / (inverse_chi + 3 / 2)
        )

    scaled_storage = _scaled_storage_samples(repeater_nodes)
    return float(np.exp(-scaled_storage / chi).mean())


@lru_cache(maxsize=None)
def _scaled_storage_samples(repeater_nodes: int) -> np.ndarray:
    data = MonteCarloSim.load_freq_dist(repeater_nodes=repeater_nodes)
    samples = np.asarray(data["W"], dtype=float)
    samples.setflags(write=False)
    return samples


def homogeneous_chain(
    *,
    repeater_nodes: int,
    total_distance_km: float,
    excitation_probability: float,
    detector_efficiency: float,
    dark_count_probability: float,
    transmission_without_fiber: float,
    attenuation_db_per_km: float,
    visibility: float,
    coherence_time: float,
    swapping_quality: float,
    local_cycle_time_s: float = 50e-6,
    phase_coherence: float = 1.0,
) -> ChainResult:
    """Return fidelity and exact rate for a homogeneous symmetric chain."""

    if repeater_nodes < 0:
        raise ValueError("repeater_nodes cannot be negative")
    if not 0 <= swapping_quality <= 1:
        raise ValueError("swapping_quality must lie in [0, 1]")
    link_count = repeater_nodes + 1
    link_distance = total_distance_km / link_count
    eta = fiber_transmission(
        link_distance / 2,
        transmission_without_fiber=transmission_without_fiber,
        attenuation_db_per_km=attenuation_db_per_km,
    )
    link = symmetric_elementary_link(
        excitation_probability=excitation_probability,
        transmission=eta,
        detector_efficiency=detector_efficiency,
        dark_count_probability=dark_count_probability,
        visibility=visibility,
        phase_coherence=phase_coherence,
    )
    frequency = repetition_rate(link_distance, local_cycle_time_s=local_cycle_time_s)
    expected_rounds = expected_completion_round([link.heralding_probability] * link_count)
    rate = frequency / expected_rounds
    kappa = continuum_coherence_factor(
        repeater_nodes,
        heralding_probability=link.heralding_probability,
        attempt_rate=frequency,
        coherence_time=coherence_time,
    )
    fidelity = (
        (1 + swapping_quality**repeater_nodes * link.delta**link_count) / 4
        + swapping_quality**repeater_nodes
        * kappa
        * (link.signal - link.error) ** link_count
        / 2
    )
    return ChainResult(
        fidelity=float(fidelity),
        rate=float(rate),
        coherence_factor=float(kappa),
        elementary_link=link,
        repetition_rate=float(frequency),
    )
