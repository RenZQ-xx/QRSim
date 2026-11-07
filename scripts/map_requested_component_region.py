"""Map specified open-ended hardware ranges into the (Delta, chi) plane.

The p_e direction is eliminated analytically, not truncated in a random scan.
The remaining envelope is a convergent V/dark-count grid optimization.
"""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from QRSim.analytical import symmetric_elementary_link

ETA0_MAX = 0.08

def coefficients(eta,dc,V):
    n0=1-dc;n1=(1-dc)*(1-.9*(1-dc));n2=(1-dc)*(1-.9*(1-dc))**2
    c0=1-n0;c1=1-n1;c2=1-n2
    A=n0*c0
    B=eta*(n0*c1+n1*c0)+2*A*(1-eta)
    D=A*(1-eta)**2+(n0*c1+n1*c0)*eta*(1-eta)+(n0*c2+n2*c0)*(1+V)/4*eta**2+n1*c1*(1-V)/2*eta**2
    H=eta*np.sqrt(V)*(n0*c1-n1*c0)
    return A,B,D,H


def boundary(chis,nv,nd,eta_max=ETA0_MAX):
    vs,dcs=np.meshgrid(np.linspace(.9,1,nv),np.geomspace(1e-8,1e-5,nd))
    vs=vs.ravel();dcs=dcs.ravel()
    n0=1-dcs;n1=n0*(1-.9*(1-dcs))
    c0=1-n0
    A=n0*c0;qp=n0*(1-n1)+n1*c0;qm=n0*(1-n1)-n1*c0
    rows=[]
    for chi in chis:
        # Include the maximum-transmission intersection explicitly; a uniform V grid
        # otherwise misses the narrow high-chi edge near V=1.
        edge_v=(chi*(eta_max*qp+2*A*(1-eta_max))/(eta_max*qm))**2
        vv=np.r_[vs,np.clip(edge_v,.9,1)]
        darks=np.r_[dcs,dcs]
        AA=np.r_[A,A];QP=np.r_[qp,qp];QM=np.r_[qm,qm]
        denominator=np.sqrt(vv)*QM-chi*(QP-2*AA)
        eta=np.divide(2*chi*AA,denominator,out=np.full_like(AA,np.nan),where=denominator>0)
        valid=(eta>0)&(eta<=eta_max*(1+1e-8))
        eta=np.where(valid,np.minimum(eta,eta_max),np.nan)
        aa,bb,dd,hh=coefficients(eta,darks,vv)
        # x=p_e/(1-p_e); maximize Delta by minimizing A/x+D*x.
        xx=np.minimum(np.sqrt(aa/dd),1/9)
        ratio=(aa/xx+dd*xx)/bb
        delta=(1-ratio)/(1+ratio)
        delta[~valid]=-np.inf
        if not np.any(valid):continue
        j=np.argmax(delta)
        rows.append(dict(chi=float(chi),delta_max=float(delta[j]),eta=float(eta[j]),V=float(vv[j]),p_dc=float(darks[j]),p_e=float(xx[j]/(1+xx[j])),L_el_km=0.,eta0=float(eta[j])))
    return rows


def save(path,rows):
    with path.open('w',newline='',encoding='utf-8') as s:
        w=csv.DictWriter(s,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def main():
    out=Path(__file__).resolve().parents[1]/'results/state_space_exploration/requested_component_region'
    out.mkdir(parents=True,exist_ok=True)
    _,bmax,_,hmax=coefficients(ETA0_MAX,1e-8,1.)
    chi_max=hmax/bmax
    chi_grid=np.unique(np.r_[np.linspace(.0001,.94,401),
        np.linspace(.9401,.99999,601),np.linspace(.99999,chi_max,12)])
    fine=boundary(chi_grid,801,121)
    coarse=boundary(chi_grid,401,61)
    error=max(abs(a['delta_max']-b['delta_max']) for a,b in zip(fine,coarse))
    print('maximum envelope change on grid refinement',error)
    save(out/'reachable_boundary.csv',fine)
    for r in fine[::25]:
        state=symmetric_elementary_link(excitation_probability=r['p_e'],transmission=r['eta'],detector_efficiency=.9,dark_count_probability=r['p_dc'],visibility=r['V'])
        assert abs(state.delta-r['delta_max'])<1e-10
        assert abs(state.coherence/state.population_single_excitation-r['chi'])<1e-10
    # Deterministic hardware samples for occupied-region validation, not a prior.
    rng=np.random.default_rng(20261105);count=100000
    eta0=10**rng.uniform(-8,np.log10(ETA0_MAX),count)
    lengths=rng.uniform(0,250,count)
    eta=eta0*10**(-.2*lengths/20)
    dc=10**rng.uniform(-8,-5,count);V=rng.uniform(.9,1,count)
    pe=10**rng.uniform(-8,-1,count)
    A,B,D,H=coefficients(eta,dc,V)
    den=A*(1-pe)**2+B*pe*(1-pe)+D*pe**2
    delta=(B*pe*(1-pe)-A*(1-pe)**2-D*pe**2)/den
    chi=H/B;F=(1+delta)*(1+chi)/4
    w=(4*F-1)/3;offset=delta-w
    rows=[dict(eta0=eta0[i],L_el_km=lengths[i],p_dc=dc[i],V=V[i],p_e=pe[i],delta=delta[i],chi=chi[i],F_init=F[i],delta_minus_w=offset[i]) for i in range(count)]
    save(out/'hardware_samples.csv',rows)
    for i in rng.choice(count,40,replace=False):
        state=symmetric_elementary_link(excitation_probability=pe[i],transmission=eta[i],detector_efficiency=.9,dark_count_probability=dc[i],visibility=V[i])
        assert abs(state.delta-delta[i])<1e-12
        assert abs(state.fidelity-F[i])<1e-12
    useful=F>=.5
    for label,mask in [('all',np.ones(count,dtype=bool)),('F_init>=0.5',useful)]:
        print(label,'sample count',mask.sum(),'Delta range',delta[mask].min(),delta[mask].max(),'chi range',chi[mask].min(),chi[mask].max(),'max |Delta-w|',abs(offset[mask]).max())
    fig,axes=plt.subplots(1,2,figsize=(12.8,5.2),layout='constrained')
    ys=np.array([r['chi'] for r in fine]);ds=np.array([r['delta_max'] for r in fine])
    # At eta->0, chi->0 and p_e->0.1, Delta_max->-0.64.
    ys=np.r_[0,ys];ds=np.r_[-.64,ds]
    outside=delta-np.interp(chi,ys,ds)
    assert outside.max()<1e-4, "Hardware sample outside the numerical envelope"
    print("Maximum sampled excess over interpolated envelope",outside.max())
    for ax in axes:
        ax.fill_betweenx(ys,-1,ds,color='#bfd8e7',alpha=.75,label='Component-accessible region (closure)')
        ax.plot(ds,ys,color='#0072B2',lw=1.4,label='Numerical upper envelope')
        select=np.arange(0,count,4)
        ax.scatter(delta[select],chi[select],s=.25,c='#164a62',alpha=.08,rasterized=True)
        x=np.linspace(0,1,1000)
        ax.plot(x,2*x/(1+x),'--',color='#e76f51',lw=2,label='Werner locus')
        xx,yy=np.meshgrid(np.linspace(-1,1,401),np.linspace(0,1,301))
        ctr=ax.contour(xx,yy,(1+xx)*(1+yy)/4,levels=[.5,.7,.8,.9],colors='0.5',linewidths=.7)
        ax.clabel(ctr,fmt='F=%.1f',fontsize=8)
        ax.set(xlabel=r'Population contrast $\Delta$',ylabel=r'Normalized coherence $\chi=h/b$')
        ax.grid(alpha=.12)
    axes[0].set(xlim=(-1,1),ylim=(0,1.005),title='(a) Full state-space range')
    axes[1].set(xlim=(.2,1),ylim=(.4,1.005),title='(b) Original plotting window')
    axes[0].legend(loc='upper left',fontsize=8,framealpha=.95)
    fig.savefig(out/'component_accessible_state_space.png',dpi=220)
    fig.savefig(out/'component_accessible_state_space.svg')
    plt.close(fig)
    for target in [.5,.6842105263,.8,.9,.95,.99]:
        r=min(fine,key=lambda r:abs(r['chi']-target));print('BOUNDARY',r)
    print(out/'component_accessible_state_space.png')


if __name__=='__main__':main()
