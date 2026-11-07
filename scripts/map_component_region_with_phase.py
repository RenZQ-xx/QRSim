"""Extend the hardware-accessible region by Gaussian relative phase noise.

Run map_requested_component_region.py first to regenerate the zero-phase
numerical boundary. Phase noise changes only h, not diagonal probabilities.
"""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from QRSim.analytical import symmetric_elementary_link
from map_requested_component_region import save


def main():
    out=Path(__file__).resolve().parents[1]/'results/state_space_exploration/requested_component_region'
    with (out/'reachable_boundary.csv').open(encoding='utf-8') as s:
        base=list(csv.DictReader(s))
    y=np.r_[0.,[float(r['chi']) for r in base]]
    d=np.r_[-.64,[float(r['delta_max']) for r in base]]
    nu_min=np.exp(-.5*.5**2)
    # The phase-extended right edge at chi is max D0(t),
    # t in [chi, min(chi/nu_min, chi0_max)]. Evaluate this exactly for the
    # piecewise-linear numerical input envelope, including interval endpoints.
    expanded=[];source_chi=[]
    for chi in y:
        upper=min(chi/nu_min,y[-1])
        interior=(y>=chi)&(y<=upper)
        xs=np.r_[chi,y[interior],upper]
        ds=np.interp(xs,y,d)
        best=np.argmax(ds)
        expanded.append(ds[best]);source_chi.append(xs[best])
    expanded=np.array(expanded);source_chi=np.array(source_chi)
    assert np.all(expanded>=d-1e-12)
    rows=[]
    for chi,d0,dp,src in zip(y,d,expanded,source_chi):
        sigma=0. if src==0 else np.sqrt(max(0.,-2*np.log(chi/src)))
        assert sigma<=.5+1e-10
        rows.append(dict(chi=chi,delta_max_no_phase=d0,delta_max_with_phase=dp,source_chi=src,source_sigma_rad=sigma))
    save(out/'reachable_boundary_phase_0_05rad.csv',rows)
    # Apply phase channel to the existing reproducible hardware samples.
    with (out/'hardware_samples.csv').open(encoding='utf-8') as s:
        samples=list(csv.DictReader(s))
    rng=np.random.default_rng(20261106)
    sigma=rng.uniform(0,.5,len(samples))
    phase_rows=[]
    for i,r in enumerate(samples):
        chi0=float(r['chi']);delta=float(r['delta']);nu=np.exp(-sigma[i]**2/2)
        chi=nu*chi0;f=(1+delta)*(1+chi)/4
        phase_rows.append(dict(eta0=float(r['eta0']),L_el_km=float(r['L_el_km']),p_dc=float(r['p_dc']),V=float(r['V']),p_e=float(r['p_e']),sigma_phi_rad=sigma[i],nu_phi=nu,delta=delta,chi_without_phase=chi0,chi=chi,F_init=f,delta_minus_w=delta-(4*f-1)/3))
    save(out/'hardware_samples_phase_0_05rad.csv',phase_rows)
    sampled_delta=np.array([r['delta'] for r in phase_rows])
    sampled_chi=np.array([r['chi'] for r in phase_rows])
    excess=sampled_delta-np.interp(sampled_chi,y,expanded)
    assert excess.max()<1e-4, "Phase-noisy hardware sample outside the envelope"
    print('Maximum phase-sampled excess over envelope',excess.max())
    for i in rng.choice(len(samples),40,replace=False):
        r=phase_rows[i]
        state=symmetric_elementary_link(excitation_probability=r['p_e'],transmission=r['eta0']*10**(-.2*r['L_el_km']/20),detector_efficiency=.9,dark_count_probability=r['p_dc'],visibility=r['V'])
        rho=state.density_matrix();diagonal=np.diag(rho).copy()
        rho[1,2]*=r['nu_phi'];rho[2,1]*=r['nu_phi']
        assert np.allclose(np.diag(rho),diagonal,atol=0,rtol=0)
        assert abs(np.trace(rho)-1)<1e-12
        assert np.linalg.eigvalsh(rho).min()>-1e-12
        assert abs((rho[1,1]+rho[2,2]+rho[1,2]+rho[2,1])/2-r['F_init'])<1e-12
    fig,axes=plt.subplots(1,2,figsize=(12.8,5.2),layout='constrained')
    for ax in axes:
        ax.fill_betweenx(y,-1,expanded,color='#f2c88e',alpha=.85,label=r'Accessible region: $0\leq\sigma_\phi\leq0.5$ rad')
        ax.plot(expanded,y,color='#b66b16',lw=1.5)
        x=np.linspace(0,1,1000)
        ax.plot(x,2*x/(1+x),'--',color='#d75540',lw=2,label='Werner locus')
        X,Y=np.meshgrid(np.linspace(-1,1,401),np.linspace(0,1,301))
        ctr=ax.contour(X,Y,(1+X)*(1+Y)/4,levels=[.5,.7,.8,.9],colors='0.45',linewidths=.7)
        ax.clabel(ctr,fmt='F=%.1f',fontsize=8)
        ax.set(xlabel=r'Population contrast $\Delta$',ylabel=r'Normalized coherence $\chi=h/b$')
        ax.grid(alpha=.12)
    axes[0].set(xlim=(-1,1),ylim=(0,1.005),title='(a) Full state-space range')
    axes[1].set(xlim=(.2,1),ylim=(.4,1.005),title='(b) Original plotting window')
    axes[0].legend(loc='upper left',fontsize=8,framealpha=.95)
    fig.savefig(out/'component_accessible_state_space_phase_0_05rad_orange.png',dpi=220)
    fig.savefig(out/'component_accessible_state_space_phase_0_05rad_orange.svg')
    plt.close(fig)
    for delta in (.5,.8,.9,.95):
        before=y[np.flatnonzero(d>=delta)[0]]
        after=y[np.flatnonzero(expanded>=delta)[0]]
        print('Delta',delta,'min chi without phase',before,'with phase',after)
    delta=.9;chi=4*.8/(1+delta)-1
    print('Target F=.8 Delta=.9 required chi',chi,'max accessible Delta at this chi',np.interp(chi,y,expanded))
    print('nu_min',nu_min,'40 phase-channel state checks passed')


if __name__=='__main__':main()
