# QRSim v1.0.0 reproducibility guide

This package reproduces every computational result used by the manuscript
“Analytical modeling of single-photon-interference-based quantum repeaters.”
Historical figures, cutoff simulations, IDE files, and development-only audit
notes are intentionally excluded.

## Environment

The verified environment is Python 3.11 with the exact package versions in
`requirements.txt`.

```powershell
python -m venv venv
.\venv\Scripts\python.exe -m pip install -r requirements.txt
.\venv\Scripts\python.exe -m pip install -e . --no-deps
```

On Linux or macOS, replace the final two commands with the corresponding
`venv/bin/python` paths.

## Reproduce the manuscript figures

```powershell
.\venv\Scripts\python.exe scripts\reproduce_figures.py
```

The command creates:

| Manuscript result | Generated file | Calculation |
| --- | --- | --- |
| Fig. 5 | `results/figures/fig4_reproduced.pdf` | Asymmetric one-repeater placement |
| Fig. 6 | `results/figures/fig2_reproduced.pdf` | Distance scan for `N=1,3,5` |
| Fig. 7 | `results/figures/fig3_reproduced.pdf` | 1000-km rate--fidelity scan |
| Fig. 8 | `results/figures/Figure_xx_reproduced.pdf` | Continuum-density Monte Carlo check |

CSV data for Figs. 5--7 are written alongside the figures. The corresponding
verified outputs are supplied under `reference_outputs`.

## Monte Carlo provenance

Large derived arrays are not distributed. For each repeater-node count `N`,
the figure script directly samples 500,000 sets of `N+1` independent unit-rate
exponential variables with seed `20260728 + N`. The resulting continuum
completion and storage variables are cached locally only to avoid repeated
generation during one run. Deleting `src/QRSim/datas/distribution_*.npz` and
rerunning the command regenerates them exactly.

## Validate the continuum approximation

```powershell
.\venv\Scripts\python.exe scripts\validate_continuum_approximation.py
```

This writes `results/validation/continuum_validation.csv`. It compares the
continuum coherence factor with exact finite-success-probability geometric
waiting times at the five operating points in manuscript Table I. Each point
uses 1,000,000 paired samples and seed `20260728 + N`. The expected maximum
absolute changes are

- `|Delta kappa_N| < 5.5e-6`;
- `|Delta F| < 1.3e-6`.

## Run the checks

```powershell
.\venv\Scripts\python.exe -m pytest
```

The tests verify the analytical elementary-link model against the QuTiP
component model, the X-state swapping recursion, exact completion-time
expectations, continuum moments, and the arbitrary-length Laplace formula for
the first three explicit cases.

