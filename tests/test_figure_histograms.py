import numpy as np

from scripts.reproduce_figures import full_sample_histogram_weights


def test_truncated_histogram_preserves_full_sample_probability_mass():
    samples = np.array([0.1, 0.3, 0.7, 1.2])
    bins = np.array([0.0, 0.5, 1.0])

    density, _ = np.histogram(
        samples,
        bins=bins,
        weights=full_sample_histogram_weights(samples, bins),
    )

    # Three of four samples lie in the displayed interval.  Its plotted area
    # must therefore be 3/4, rather than being renormalized to one.
    assert np.sum(density * np.diff(bins)) == 0.75
