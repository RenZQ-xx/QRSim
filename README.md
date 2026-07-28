# QRSim

QRSim implements the component-to-chain model for single-photon-interference
quantum repeaters described in `main.tex`.

## Install

From the repository root:

```powershell
python -m pip install -r requirements.txt
python -m pip install -e . --no-deps
```

Run the active analytical and waiting-time tests with:

```powershell
python -m pytest
```

## Reproduce recovered figures

The command below writes new files under `results/figures`; it does not
overwrite historical manuscript images.

```powershell
python scripts/reproduce_figures.py fig4 fig2 fig3 Figure_xx
```

- `fig2` and `fig3`: homogeneous-chain distance and rate--fidelity results.
- `Figure_xx`: a deterministic Monte Carlo check of the continuum waiting-time
  densities.
- `fig4`: uses the active manuscript parameters and solves the central
  excitation probabilities from the heralding probabilities stated in the
  caption.

Large Monte Carlo arrays are not distributed. If absent, they are regenerated
from unit-rate exponential random variables using 500,000 trials and the fixed
seed `20260728 + N`, where `N` is the number of repeater nodes.

## Validate the continuum approximation

The following command compares the continuum coherence factor with an exact
finite-success-probability geometric waiting-time simulation at the operating
points reported in the manuscript:

```powershell
python scripts/validate_continuum_approximation.py
```

It writes `results/validation/continuum_validation.csv`. The validation uses
common random numbers and a documented fixed seed so that both the estimate
and its standard error are reproducible.

Cutoff figures are intentionally outside the restoration scope because the
cutoff subsection has been removed from the active manuscript.
