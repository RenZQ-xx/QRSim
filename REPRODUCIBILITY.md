# QRSim v1.1.0: PRA revision reproducibility

This version supports AV12595 and its revised numerical comparisons.
Release: https://github.com/RenZQ-xx/QRSim/releases/tag/v1.1.0
The original v1.0.1 release is preserved; v1.1.0 adds the revised baseline,
residual phase noise, reduced-model comparisons, accessible-region mapping,
component sensitivity, and seeded validation.

## Environment and commands

Use Python 3.10 or later in a writable project directory. The checked run
used Python 3.11.7 and QuTiP 5.1.1. Exact package versions are recorded in
environment_checked.txt in the package. Run from the extracted package root:

```text
python -m pip install -e ".[test]"
python -m pytest
python scripts/reproduce_figures.py fig4 fig2 fig3 Figure_xx
python scripts/validate_continuum_approximation.py
python scripts/reproduce_pra_revision.py
```

No cached Monte Carlo arrays are required: missing arrays are regenerated
in src/QRSim/datas using 500,000 unit-rate exponential samples with seed
20260728 + N. A full run takes several minutes and has no network dependency
after Python packages are installed. The commands above, not the historical
exploratory scripts, are the authoritative pipeline.

## Figure mapping

Static scheme drawings are not simulation outputs. Computational figures:

| Content | Output |
| --- | --- |
| Asymmetric placement | results/figures/fig4_reproduced.pdf |
| Distance fidelity/rate | results/figures/fig2_reproduced.pdf |
| Excitation/node-count rate-fidelity | results/figures/fig3_reproduced.pdf |
| Waiting-time histogram checks | results/figures/Figure_xx_reproduced.pdf |
| State-space paths and envelope | results/revision/state_space_envelope.pdf |
| Orange phase-extended accessible region | results/revision/component_accessible_region.pdf |
| Auxiliary maximum ambiguity (not used in the current paper) | results/revision/maximum_gap_summary.pdf |

New CSV outputs cover distance/excitation comparisons, the false-herald sweep,
representative cases, state-envelope samples, boundary curves, maximum gaps,
component sensitivity, and validation. LaTeX tables are generated from CSV.
The former standalone hardware-comparison figure is omitted from the
manuscript and current package. Fig. 9 reads markers I--III directly from
results/revision/representative_model_cases.csv, matching Table III rows.
The hardware maps also regenerate 100,000 independent check points under
results/state_space_exploration. Their density is not a physical prior.
The interface-efficiency range is 0 < eta0 <= 0.08, including the main
baseline efficiency. Other mapping ranges remain V in [0.9,1], detector
efficiency 0.9, dark-count probability in [1e-8,1e-5], link length in
[0,250] km, 0 < p_e < 0.1, and phase standard deviation in [0,0.5] rad.

## Model scope

The revised baseline uses V=0.99 and an independent zero-mean Gaussian
relative phase with sigma_phi=0.30 rad, hence phase_coherence=exp(-0.045).
All baseline parameters are owned by MANUSCRIPT_PARAMETERS in
scripts/reproduce_figures.py; validation and revision scripts import it.
The earlier V=0.9, zero-phase-noise results are retained locally under
output/revision/before_high_visibility_baseline and are not the current
manuscript benchmark. The separate high-visibility controls use
sigma_phi=0.10 rad in cases I and II and 0.35 rad in case III.
State-space envelope figures keep
their displayed ideal-memory conditions and are not hardware benchmarks.

All reduced models preserve elementary fidelity, heralding probability,
waiting statistics, pure-memory-dephasing model, and effective swap quality
in each paired comparison. Werner imposes equal non-target Bell weights;
phase-only/bit-only put the infidelity in the corresponding single error.
Strict Pauli twirling is distinct and fidelity preserving under the stated
averaged, non-postselected swapping assumptions.

Independent Gaussian relative phases attenuate only coherence by
exp(-sigma_phi**2/2); they do not replace visibility in two-photon terms.
Correlated drift is not modeled. The controlled examples use N=2, 25-km
links, p_e=0.001, eta0=0.019, V=0.999, detector efficiency 0.9,
coherence time 600 s, and swap quality 0.995. I and II differ only in
dark-count probability. II and III share the same 20% dark-only false-herald
fraction and differ only in phase noise. The dark-free reference lies outside the positive
dark-count stress-test range.

## Accuracy and seeds

Completion-time expectations and rates are finite-probability exact.
N<=3 uses explicit continuum transforms for plotted memory factors;
N>3 uses seeded continuum averages. The original finite-p validation uses
1,000,000 paired trials and seed 20260728+N. New controls use 500,000
paired trials and seed 20261105. Hardware samples use seed 20261105;
phase samples use 20261106.

Sampling uncertainty is distinct from continuum bias (Appendix I).
Sharing the memory factor makes the sampling error of the paired model
difference much smaller than that of either absolute fidelity.
Floating-point/plot-font differences can affect last digits or PDF bytes.

## Optional manuscript build

```text
python scripts/build_pra_submission.py --pdflatex /path/to/pdflatex
```

This needs REVTeX 4.2, the LaTeX packages declared in the source, and the
static figure assets. It builds clean and marked manuscripts, a
response letter, and a change list. Unchanged relocated appendix material
remains black; deletions are documented separately. No script pushes to GitHub,
publishes a release, or submits to APS.

In the split Sections II and III, unchanged original prose, equations, and
Protocol 1 remain black even when relocated. Only actual additions/rewrites,
new headings, and cross-reference repairs are blue there; other regions use
the stated revised-block convention.

Section IV compares prose sentences with the submitted analytical text
stored in docs/revision/submitted_analytical_baseline.txt. The 15 approved
retained display expressions are fingerprinted in
docs/revision/analytical_retained_displays.json: they remain black when
reformatted or regrouped, but a later formula change is marked blue rather
than silently inheriting that exception. The storage-parameter rename is
marked separately. The full submitted baseline is also used to distinguish
unchanged relocated appendix material from new additions. These reference
files are included in the local source archive.

The public numerical archive does not include manuscript text, static
drawings, referee reports, or response letters. The optional manuscript
build therefore requires the authors' separately maintained source package;
the public archive alone regenerates all computational figures and tables.
