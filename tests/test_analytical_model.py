import math

from QRSim import HeraldedLink
from QRSim.analytical import fiber_transmission, homogeneous_chain, symmetric_elementary_link


def test_symmetric_elementary_link_matches_qutip_implementation():
    parameters = dict(
        excitation_probability=0.037,
        detector_efficiency=0.9,
        dark_count_probability=1.5e-6,
        visibility=0.9,
    )
    eta = fiber_transmission(
        62.5,
        transmission_without_fiber=0.08,
        attenuation_db_per_km=0.2,
    )
    analytical = symmetric_elementary_link(transmission=eta, **parameters)

    qutip_link = HeraldedLink(
        detector_efficiency=parameters["detector_efficiency"],
        dark_count_prob=parameters["dark_count_probability"],
        transmission_without_fiber=0.08,
        distance=125.0,
        fiber_attenuation_coefficient=0.2,
        visibility=parameters["visibility"],
    )
    qutip_link.generate_entanglement(parameters["excitation_probability"])
    assert math.isclose(analytical.heralding_probability, qutip_link.herald_prob, rel_tol=1e-12)
    assert math.isclose(analytical.fidelity, qutip_link.initial_fidelity_set["psi+"], rel_tol=1e-12)


def test_paper_scale_chain_result_is_physical():
    result = homogeneous_chain(
        repeater_nodes=6,
        total_distance_km=1000,
        excitation_probability=0.045,
        detector_efficiency=0.9,
        dark_count_probability=1.5e-6,
        transmission_without_fiber=0.09,
        attenuation_db_per_km=0.17,
        visibility=0.9,
        coherence_time=600,
        swapping_quality=0.995,
    )
    assert 0 <= result.fidelity <= 1
    assert result.rate > 0


def test_asymmetric_arm_parameters_are_applied_to_the_correct_arm():
    link = HeraldedLink(
        detector_efficiency=0.9,
        dark_count_prob=1.5e-6,
        transmission_without_fiber=[0.08, 0.12],
        distance=[30, 70],
        fiber_attenuation_coefficient=[0.2, 0.17],
        visibility=0.9,
    )
    assert math.isclose(link.transmission[0], 0.08 * 10 ** (-0.2 * 30 / 10))
    assert math.isclose(link.transmission[1], 0.12 * 10 ** (-0.17 * 70 / 10))


def test_generation_round_dephasing_is_opt_in():
    common = dict(
        detector_efficiency=0.9,
        dark_count_prob=1.5e-6,
        transmission_without_fiber=0.08,
        distance=100,
        visibility=0.9,
        coherence_time=1e-3,
    )
    manuscript_link = HeraldedLink(**common)
    dephased_link = HeraldedLink(**common, include_generation_round_dephasing=True)
    manuscript_link.generate_entanglement(0.03)
    dephased_link.generate_entanglement(0.03)
    assert manuscript_link.initial_fidelity_set["psi+"] > dephased_link.initial_fidelity_set["psi+"]
