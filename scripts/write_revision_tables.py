"""Generate manuscript tables from reproducible CSV, not hand-entered numbers."""
from pathlib import Path
import csv

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/revision'


def read(name):
    with (OUT/name).open(encoding='utf-8') as s:
        return list(csv.DictReader(s))


def main():
    cases=read('representative_model_cases.csv')
    labels=dict(dark_free='No dark counts',false_20percent='$20\\%$ false heralds',near_threshold='Threshold example')
    lines=[r'\begin{table}[tbp]',
        r'\caption{\label{tab:werner-cases}Component-model comparisons for $N=2$, $L_{\rm el}=25$ km, $p_e=0.001$, $\eta_0=0.019$, $V=0.999$, $T_{\rm coh}=600$ s, $\eta_d=0.9$, and $\lambda=0.995$. Residual relative phases are independent zero-mean Gaussian variables. Both models have the same elementary fidelity and timing in each row.}',
        r'\setlength{\tabcolsep}{2.3pt}',
        r'\begin{ruledtabular}',r'\begin{tabular}{ccccccc}',
        r'Case & $\sigma_\phi$ (rad) & $10^6p_{\rm dc}$ & $f_{\rm false}^{(0)}$ & $F_X$ & $F_W$ & $|F_W-F_X|$\\',r'\hline']
    markers=dict(dark_free='I',false_20percent='II',near_threshold='III')
    for r in cases:
        if r['case']=='paper_baseline':continue
        lines.append(f"{markers[r['case']]} & {float(r['sigma_phi_rad']):.2f} & {float(r['dark_count_probability'])*1e6:.4f} & {100*float(r['false_fraction']):.2f}\\% & {float(r['F_X']):.5f} & {float(r['F_W']):.5f} & {abs(float(r['delta_F'])):.5f}"+r'\\')
    lines += [r'\end{tabular}',r'\end{ruledtabular}',r'\end{table}']
    (OUT/'model_cases_table.tex').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    sensitivity=read('component_sensitivity.csv')
    labels=dict(baseline='Baseline',no_dark_counts='$p_{\\rm dc}=0$',perfect_visibility='$V=1$',perfect_phase='$\\sigma_\\phi=0$',ideal_memory='$T_{\\rm coh}\\to\\infty$',perfect_swaps='$\\lambda=1$',ideal_detector='$\\eta_d=1$',improved_optical_interface='$\\eta_0=0.16$')
    lines=[r'\begin{table}[tbp]',
        r'\caption{\label{tab:component-sensitivity}One-at-a-time component diagnostics at the representative working point of Fig.~\ref{fig:fig3}: 1000 km, $N=6$, $p_e=0.036$, and the parameters in Table~\ref{tab:baseline-parameters}, including $V=0.99$, $\sigma_\phi=0.30$ rad, and $\eta_0=0.08$. No reoptimization is performed. Differences are not additive.}',
        r'\begin{ruledtabular}',r'\begin{tabular}{lrrrr}',
        r'Modification & $F_X$ & $\Delta F_X$ & $R$ (Hz) & $R/R_0$\\',r'\hline']
    for r in sensitivity:
        lines.append(f"{labels[r['case']]} & {float(r['fidelity']):.5f} & {float(r['delta_fidelity']):+.5f} & {float(r['rate_hz']):.5f} & {float(r['rate_ratio']):.4f}"+r'\\')
    lines += [r'\end{tabular}',r'\end{ruledtabular}',r'\end{table}']
    (OUT/'component_sensitivity_table.tex').write_text('\n'.join(lines)+'\n',encoding='utf-8')
    validation=read('new_case_continuum_validation.csv')
    lines=[r'\begin{table}[tbp]',
        r'\caption{\label{tab:new-continuum-validation}Paired continuum/geometric checks of the $N=2$ cases in Table~\ref{tab:werner-cases}. Each uses $500\,000$ trials and seed $20261105$. Bias is continuum minus finite-$p_{\rm hr}$; the gap is $F_W-F_X$.}',
        r'\begin{ruledtabular}',r'\begin{tabular}{lrrr}',
        r'Case & $10^6\Delta F_X$ & $10^6\Delta(F_W-F_X)$ & $10^6{\rm SE}_{\Delta\kappa}$\\',r'\hline']
    labels=dict(dark_free='I',false_20percent='II',near_threshold='III')
    for r in validation:
        if r['case']=='paper_baseline':continue
        lines.append(f"{labels[r['case']]} & {float(r['continuum_bias_F_X'])*1e6:+.3f} & {float(r['continuum_bias_model_gap'])*1e6:+.3f} & {float(r['paired_se'])*1e6:.4f}"+r'\\')
    lines += [r'\end{tabular}',r'\end{ruledtabular}',r'\end{table}']
    (OUT/'new_continuum_table.tex').write_text('\n'.join(lines)+'\n',encoding='utf-8')


if __name__=='__main__':main()
