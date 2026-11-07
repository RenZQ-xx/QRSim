import math

import numpy as np
import pytest

from QRSim import HeraldedLink
from QRSim.analytical import (
    bell_diagonal_chain_fidelity,
    fiber_transmission,
    fidelity_matched_bell_state,
    homogeneous_chain,
    reduced_model_fidelity,
    symmetric_elementary_link,
)


@pytest.mark.parametrize("visibility,phase_coherence", [(0.9, 1.0), (0.99, math.exp(-0.30**2/2)), (0.99, 0.0)])
def test_symmetric_elementary_link_matches_qutip_implementation(visibility, phase_coherence):
    parameters = dict(
        excitation_probability=0.037,
        detector_efficiency=0.9,
        dark_count_probability=1.5e-6,
        visibility=visibility,
        phase_coherence=phase_coherence,
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
        phase_coherence=phase_coherence,
    )
    qutip_link.generate_entanglement(parameters["excitation_probability"])
    assert math.isclose(analytical.heralding_probability, qutip_link.herald_prob, rel_tol=1e-12)
    assert math.isclose(analytical.fidelity, qutip_link.initial_fidelity_set["psi+"], rel_tol=1e-12)
    from QRSim import RelayChain
    chain = RelayChain(Links=[qutip_link, qutip_link], coherence_time=math.inf)
    chain.solve(swapping_quality=0.995)
    expected = (1 + 0.995 * analytical.delta**2) / 4 + 0.995 * (2*analytical.coherence)**2 / 2
    assert chain.final_fide == pytest.approx(expected, abs=1e-12)


def test_revised_baseline_and_strict_twirl():
    from scripts.reproduce_figures import MANUSCRIPT_PARAMETERS
    assert MANUSCRIPT_PARAMETERS["visibility"] == 0.99
    assert MANUSCRIPT_PARAMETERS["phase_coherence"] == math.exp(-0.30**2/2)
    result = homogeneous_chain(repeater_nodes=6, total_distance_km=1000,
                               excitation_probability=0.036, **MANUSCRIPT_PARAMETERS)
    assert result.fidelity == pytest.approx(0.5905896738759612, abs=1e-12)
    assert result.rate == pytest.approx(0.0990267741, abs=1e-9)
    rho = result.elementary_link.density_matrix()
    assert np.trace(rho) == pytest.approx(1.0)
    assert np.linalg.eigvalsh(rho).min() >= -1e-12
    twirled = bell_diagonal_chain_fidelity(
        result.elementary_link.strict_pauli_twirl(), repeater_nodes=6,
        coherence_factor=result.coherence_factor, swapping_quality=0.995)
    assert twirled == pytest.approx(result.fidelity, abs=1e-14)


def test_accessible_region_includes_baseline_efficiency_and_expands():
    from scripts.map_requested_component_region import ETA0_MAX, boundary
    assert ETA0_MAX == 0.08
    chis = np.array([0.5, 0.8, 0.9, 0.95, 0.99])
    old = boundary(chis, 201, 31, eta_max=0.02)
    expanded = boundary(chis, 201, 31)
    assert len(old) == len(expanded) == len(chis)
    for a, b in zip(old, expanded):
        assert b["delta_max"] >= a["delta_max"] - 1e-5
        assert 0 < b["eta0"] <= ETA0_MAX
        state = symmetric_elementary_link(
            excitation_probability=b["p_e"], transmission=b["eta"],
            detector_efficiency=0.9, dark_count_probability=b["p_dc"],
            visibility=b["V"])
        assert state.delta == pytest.approx(b["delta_max"], abs=1e-10)


@pytest.mark.parametrize("dc,sigma,expected_x,expected_w", [
    (0.0, 0.10, 0.9557362377168104, 0.9556099048117469),
    (2.404044477595425e-6, 0.10, 0.5411553090967312, 0.5301572904248086),
    (2.404044477595425e-6, 0.35, 0.5042129754246051, 0.49807927606924823),
])
def test_table_three_nonzero_phase_cases(dc, sigma, expected_x, expected_w):
    result = homogeneous_chain(
        repeater_nodes=2, total_distance_km=75., excitation_probability=.001,
        detector_efficiency=.9, dark_count_probability=dc,
        transmission_without_fiber=.019, attenuation_db_per_km=.2,
        visibility=.999, phase_coherence=math.exp(-sigma**2/2),
        coherence_time=600., swapping_quality=.995)
    link=result.elementary_link
    assert result.fidelity == pytest.approx(expected_x, abs=1e-12)
    werner=reduced_model_fidelity(link.fidelity,'werner',repeater_nodes=2,
        coherence_factor=result.coherence_factor,swapping_quality=.995)
    assert werner == pytest.approx(expected_w, abs=1e-12)
    assert np.linalg.eigvalsh(link.density_matrix()).min() >= -1e-12
    twirl=bell_diagonal_chain_fidelity(link.strict_pauli_twirl(),repeater_nodes=2,
        coherence_factor=result.coherence_factor,swapping_quality=.995)
    assert twirl == pytest.approx(result.fidelity,abs=1e-14)
    if dc > 0:
        assert result.rate == pytest.approx(0.07492988450602923,abs=1e-12)
    if sigma == .35:
        assert werner < .5 < result.fidelity


@pytest.mark.parametrize("phase_coherence", [-0.01, 1.01, float('nan')])
def test_qutip_link_rejects_invalid_phase_coherence(phase_coherence):
    with pytest.raises(ValueError):
        HeraldedLink(detector_efficiency=.9, dark_count_prob=1.5e-6,
                     transmission_without_fiber=.08, distance=100,
                     phase_coherence=phase_coherence)


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


@pytest.fixture
def representative_link():
    eta = fiber_transmission(
        1000 / (6 + 1) / 2,
        transmission_without_fiber=0.08,
        attenuation_db_per_km=0.2,
    )
    return symmetric_elementary_link(
        excitation_probability=0.036,
        transmission=eta,
        detector_efficiency=0.9,
        dark_count_probability=1.5e-6,
        visibility=0.9,
    )


def test_component_and_reduced_states_are_normalized_and_positive(representative_link):
    states = [representative_link.density_matrix()]
    states.append(representative_link.strict_pauli_twirl().density_matrix())
    states.extend(
        fidelity_matched_bell_state(representative_link.fidelity, model).density_matrix()
        for model in ("werner", "phase_only", "bit_only")
    )
    for state in states:
        assert np.isclose(np.trace(state), 1.0, atol=1e-12)
        assert np.allclose(state, state.T.conj(), atol=1e-12)
        assert np.linalg.eigvalsh(state).min() >= -1e-12


def test_strict_twirl_preserves_the_chain_fidelity(representative_link):
    component = homogeneous_chain(
        repeater_nodes=6,
        total_distance_km=1000,
        excitation_probability=0.036,
        detector_efficiency=0.9,
        dark_count_probability=1.5e-6,
        transmission_without_fiber=0.08,
        attenuation_db_per_km=0.2,
        visibility=0.9,
        coherence_time=600,
        swapping_quality=0.995,
    )
    twirled = representative_link.strict_pauli_twirl()
    assert not np.allclose(
        representative_link.density_matrix(),
        twirled.density_matrix(),
        atol=1e-12,
    )
    twirled_fidelity = bell_diagonal_chain_fidelity(
        twirled,
        repeater_nodes=6,
        coherence_factor=component.coherence_factor,
        swapping_quality=0.995,
    )
    assert math.isclose(component.fidelity, twirled_fidelity, abs_tol=1e-14)
    assert math.isclose(
        representative_link.delta,
        twirled.parity_contrast,
        abs_tol=1e-14,
    )
    assert math.isclose(
        representative_link.signal - representative_link.error,
        twirled.phase_contrast_plus,
        abs_tol=1e-14,
    )
    assert math.isclose(
        twirled.phase_contrast_plus,
        twirled.phase_contrast_minus,
        abs_tol=1e-14,
    )


def test_reduced_models_keep_elementary_fidelity_and_representative_values(representative_link):
    component = homogeneous_chain(
        repeater_nodes=6,
        total_distance_km=1000,
        excitation_probability=0.036,
        detector_efficiency=0.9,
        dark_count_probability=1.5e-6,
        transmission_without_fiber=0.08,
        attenuation_db_per_km=0.2,
        visibility=0.9,
        coherence_time=600,
        swapping_quality=0.995,
    )
    expected = {
        "werner": 0.5864688499,
        "phase_only": 0.6403847374,
        "bit_only": 0.6309662601,
    }
    for model, expected_fidelity in expected.items():
        state = fidelity_matched_bell_state(representative_link.fidelity, model)
        assert math.isclose(state.psi_plus, representative_link.fidelity, abs_tol=1e-15)
        value = reduced_model_fidelity(
            representative_link.fidelity,
            model,
            repeater_nodes=6,
            coherence_factor=component.coherence_factor,
            swapping_quality=0.995,
        )
        assert math.isclose(value, expected_fidelity, rel_tol=1e-9, abs_tol=1e-11)


def test_phase_and_bit_models_coincide_for_ideal_memory():
    common = dict(
        elementary_fidelity=0.83,
        repeater_nodes=5,
        coherence_factor=1.0,
        swapping_quality=0.995,
    )
    phase = reduced_model_fidelity(model="phase_only", **common)
    bit = reduced_model_fidelity(model="bit_only", **common)
    assert math.isclose(phase, bit, abs_tol=1e-15)


@pytest.mark.parametrize("model", ["werner", "phase_only", "bit_only"])
def test_bell_diagonal_formula_matches_direct_pauli_convolution(model):
    state = fidelity_matched_bell_state(0.79, model)
    distribution = state.probabilities
    convolution = np.array([1.0, 0.0, 0.0, 0.0])
    for _ in range(5):
        updated = np.zeros(4)
        for left_error, left_probability in enumerate(convolution):
            for right_error, right_probability in enumerate(distribution):
                updated[left_error ^ right_error] += left_probability * right_probability
        convolution = updated
    analytical = bell_diagonal_chain_fidelity(
        state,
        repeater_nodes=4,
        coherence_factor=1.0,
        swapping_quality=1.0,
    )
    assert math.isclose(analytical, convolution[0], abs_tol=1e-15)


@pytest.mark.parametrize("excitation_probability", [0.003, 0.036, 0.10])
@pytest.mark.parametrize("transmission", [1e-5, 1e-3, 0.08])
def test_strict_twirl_invariance_over_component_parameter_grid(
    excitation_probability,
    transmission,
):
    link = symmetric_elementary_link(
        excitation_probability=excitation_probability,
        transmission=transmission,
        detector_efficiency=0.9,
        dark_count_probability=1.5e-6,
        visibility=0.9,
    )
    repeater_nodes = 4
    coherence_factor = 0.87
    swapping_quality = 0.995
    link_count = repeater_nodes + 1
    component_fidelity = (
        1
        + swapping_quality**repeater_nodes * link.delta**link_count
        + 2
        * swapping_quality**repeater_nodes
        * coherence_factor
        * (link.signal - link.error) ** link_count
    ) / 4
    twirled_fidelity = bell_diagonal_chain_fidelity(
        link.strict_pauli_twirl(),
        repeater_nodes=repeater_nodes,
        coherence_factor=coherence_factor,
        swapping_quality=swapping_quality,
    )
    assert math.isclose(component_fidelity, twirled_fidelity, abs_tol=2e-15)


@pytest.mark.parametrize("nu", [0., .5, math.exp(-.5**2/2), 1.])
def test_optical_phase_channel_changes_only_coherence(nu):
    kwargs=dict(excitation_probability=.02,transmission=.003,
                detector_efficiency=.9,dark_count_probability=1.5e-6,visibility=.95)
    original=symmetric_elementary_link(**kwargs)
    dephased=symmetric_elementary_link(**kwargs,phase_coherence=nu)
    rho=dephased.density_matrix()
    assert np.allclose(np.diag(rho),np.diag(original.density_matrix()),atol=1e-15)
    assert math.isclose(dephased.heralding_probability,original.heralding_probability,abs_tol=1e-16)
    assert math.isclose(dephased.coherence,nu*original.coherence,abs_tol=1e-15)
    assert math.isclose(np.trace(rho),1.,abs_tol=1e-14)
    assert np.linalg.eigvalsh(rho).min()>=-1e-14


@pytest.mark.parametrize("nu", [-.1,1.1])
def test_invalid_phase_coherence_is_rejected(nu):
    with pytest.raises(ValueError):
        symmetric_elementary_link(excitation_probability=.02,transmission=.003,
            detector_efficiency=.9,dark_count_probability=1e-6,visibility=.9,phase_coherence=nu)


def test_fixed_fidelity_envelope_and_sensitivity():
    for n in (1,3,5,7):
        for f in (.5,.7,.8,.9,.999):
            delta=np.linspace(2*f-1,1,1001)
            C=(4*f-1-delta)/2
            values=.25+.995**n*(delta**(n+1)+2*C**(n+1))/4
            w=(4*f-1)/3
            lower=.25+3*.995**n*w**(n+1)/4
            upper=.25+.995**n*(1+2*(2*f-1)**(n+1))/4
            assert np.all(values>=lower-1e-14)
            assert np.all(values<=upper+1e-14)
            d=(delta[0]+delta[-1])/2;c=(4*f-1-d)/2;k=.6
            def fun(x):
                return .25+.995**n*(x**(n+1)+2*k*((4*f-1-x)/2)**(n+1))/4
            derivative=.995**n*(n+1)/4*(d**n-k*c**n)
            assert math.isclose((fun(d+1e-6)-fun(d-1e-6))/2e-6,derivative,abs_tol=1e-8)


def test_appendix_b_recursion_against_full_bell_instrument():
    """Independent 16-dimensional contraction, including unequal populations."""
    rng=np.random.default_rng(20261105)
    eye=np.eye(2);z=np.diag([1.,-1.]);x=np.array([[0.,1.],[1.,0.]])
    bells=[np.array(v)/np.sqrt(2) for v in (
        [0.,1.,1.,0.],[0.,1.,-1.,0.],[1.,0.,0.,1.],[1.,0.,0.,-1.])]
    for _ in range(40):
        states=[];coefficients=[]
        for i in range(2):
            a,b,g,d=rng.dirichlet(np.ones(4))
            h=rng.uniform(-1,1)*np.sqrt(b*g)
            states.append(np.array([[a,0,0,0],[0,b,h,0],[0,h,g,0],[0,0,0,d]]))
            coefficients.append((a,b,g,d,h))
        tensor=np.kron(*states);actual=np.zeros((4,4))
        for bell,u in zip(bells,(eye,z,x,x@z)):
            contraction=np.kron(eye,np.kron(bell.reshape(1,4),eye))
            correction=np.kron(eye,u)
            branch=contraction@tensor@contraction.T
            actual+=correction@branch@correction.T
        a,b,g,d,h=coefficients[0];aa,bb,gg,dd,hh=coefficients[1]
        expected=np.diag([a*(bb+gg)+b*(aa+dd),
            a*(aa+dd)+b*(bb+gg),g*(bb+gg)+d*(aa+dd),
            g*(aa+dd)+d*(bb+gg)])
        expected[1,2]=expected[2,1]=2*h*hh
        assert np.allclose(actual,expected,atol=1e-14)
        assert np.isclose(np.trace(actual),1)
        assert np.linalg.eigvalsh(actual).min()>=-1e-14


def test_independent_link_phase_noise_preserves_timing():
    common=dict(repeater_nodes=2,total_distance_km=75.,
        excitation_probability=.001,detector_efficiency=.9,
        dark_count_probability=2.4e-6,transmission_without_fiber=.019,
        attenuation_db_per_km=.2,visibility=.999,coherence_time=600.,
        swapping_quality=.995)
    original=homogeneous_chain(**common)
    nu=np.exp(-.5**2/2)
    phased=homogeneous_chain(**common,phase_coherence=nu)
    assert original.rate==phased.rate
    assert original.coherence_factor==phased.coherence_factor
    diagonal=.25+.995**2*original.elementary_link.delta**3/4
    assert math.isclose(phased.fidelity-diagonal,
        nu**3*(original.fidelity-diagonal),abs_tol=1e-14)
