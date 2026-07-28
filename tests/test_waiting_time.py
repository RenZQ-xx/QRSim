import math

import numpy as np

from QRSim.waiting_time import (
    accumulated_storage_rounds,
    completion_rounds,
    continuum_storage_laplace,
    expected_continuum_storage,
    expected_completion_round,
    sample_continuum_waiting_times,
    sample_waiting_times,
)


def test_one_repeater_definitions_match_manuscript():
    rounds = np.array([[1, 4], [3, 3], [7, 2]])
    np.testing.assert_array_equal(completion_rounds(rounds), [4, 3, 7])
    np.testing.assert_array_equal(accumulated_storage_rounds(rounds), [6, 0, 10])


def test_general_storage_definition_matches_eq_c_wn():
    rounds = np.array([2, 7, 4, 9])
    expected = (9 - 2) + (9 - 9) + abs(7 - 2) + abs(4 - 7) + abs(9 - 4)
    assert accumulated_storage_rounds(rounds) == expected


def test_exact_two_link_expectation():
    p0, p1 = 0.2, 0.35
    denominator = p0 + p1 - p0 * p1
    expected = 1 / p0 + 1 / p1 - 1 / denominator
    assert math.isclose(expected_completion_round([p0, p1]), expected, rel_tol=1e-14)


def test_monte_carlo_agrees_with_exact_completion_time():
    probabilities = [0.08, 0.13, 0.21]
    z, w = sample_waiting_times(probabilities, trials=300_000, seed=20260721)
    assert abs(z.mean() - expected_completion_round(probabilities)) < 0.05
    assert np.all(w >= 0)


def test_continuum_storage_mean_matches_low_node_densities():
    assert math.isclose(expected_continuum_storage(1), 2.0)
    assert math.isclose(expected_continuum_storage(2), 11 / 3)
    assert math.isclose(expected_continuum_storage(3), 31 / 6)


def test_direct_continuum_sampling_matches_exact_means():
    for nodes in (1, 2, 3, 6):
        z_bar, w_bar = sample_continuum_waiting_times(
            nodes,
            trials=300_000,
            seed=20260728 + nodes,
        )
        assert abs(w_bar.mean() - expected_continuum_storage(nodes)) < 0.02
        assert np.all(z_bar >= 0)
        assert np.all(w_bar >= 0)


def test_general_laplace_formula_reproduces_pi_1_to_pi_3():
    s = 0.7
    expected = {
        1: 0.5 / (s + 0.5),
        2: (1 / 9) / (s + 0.25) + (2 / 3) / (s + 0.5) - (7 / 9) / (s + 1),
        3: (
            (4 / 15) / (s + 0.25)
            + (3 / 16) / (s + 0.5)
            + (1 / 16) / (s + 0.5) ** 2
            - (7 / 6) / (s + 1)
            + (57 / 80) / (s + 1.5)
        ),
    }
    for nodes, value in expected.items():
        assert math.isclose(
            continuum_storage_laplace(nodes, s),
            value,
            rel_tol=1e-14,
            abs_tol=1e-14,
        )
