"""Selected (a,b) figure: state-space paths and full symmetric-X envelope."""
from pathlib import Path
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from plot_fixed_delta_comparison import models


def main(output=None):
    if output is None:
        output=Path(__file__).resolve().parents[1]/'results/state_space_exploration/拟在正文中采用的图'
    output.mkdir(parents=True,exist_ok=True)
    plt.rcParams.update({'font.size':8,'axes.titlesize':9,'axes.labelsize':9,'legend.fontsize':6.8})
    fig,(left,right)=plt.subplots(1,2,figsize=(7.1,3.2),layout='constrained',gridspec_kw={'width_ratios':[1.08,1]})
    delta,chi=np.meshgrid(np.linspace(.2,1,501),np.linspace(.4,1,401))
    f_el=(1+delta)*(1+chi)/4
    heat=left.pcolormesh(delta,chi,f_el,cmap='viridis',vmin=.42,vmax=1,shading='auto',rasterized=True)
    contours=left.contour(delta,chi,f_el,levels=[.55,.65,.75,.85,.95],colors='white',linewidths=.75)
    left.clabel(contours,fmt='%.2f',fontsize=7)
    threshold=left.contour(delta,chi,f_el,levels=[.5],colors='#ffb347',linewidths=1.8)
    left.clabel(threshold,fmt={.5:r'$F_{\rm init}=0.50$'},fontsize=7,manual=[(.27,.58)])
    d=np.linspace(.25,1,500)
    left.plot(d,2*d/(1+d),'--',color='#e76f51',lw=2,label=r'Werner: $\chi=2\Delta/(1+\Delta)$')
    left.plot([.2,1],[1,1],color='black',lw=2.5,label=r'Double-excitation only: $\chi=1$')
    left.plot([1,1],[.4,1],color='#56b4e9',lw=2.5,label=r'Dephasing only: $\Delta=1$')
    left.set(xlim=(.2,1.005),ylim=(.4,1.005),xlabel=r'Population contrast $\Delta$',ylabel=r'Normalized coherence $\chi=h/b$',title=r'$N=0$: $F_{\rm end}=F_{\rm init}$')
    left.legend(loc='lower right',framealpha=.94)
    cb=fig.colorbar(heat,ax=left,pad=.02,fraction=.045,label=r'Initial fidelity $F_{\rm init}$')
    cb.set_ticks([.5,.6,.7,.8,.9,1])

    n=3;quality=.995
    f=np.linspace(.5,1,2001)
    w,t,p=models(f,n)
    right.fill_between(f,w,p,color='#c4d7e8',alpha=.65,label='Allowed symmetric X-state range',zorder=0)
    right.plot(f,w,'--',color='#e76f51',lw=2,label='Werner (lower bound)')
    right.plot(f,t,color='black',lw=2,label='Double-excitation only')
    right.plot(f,p,color='#56b4e9',lw=2,label='Dephasing only (upper bound)')
    # Exact envelope checked against a dense physical-state sweep at each F.
    fraction=np.linspace(0,1,401)
    ds=(2*f[:,None]-1)+(2-2*f[:,None])*fraction
    cs=(4*f[:,None]-1-ds)/2
    ends=.25+quality**n*(ds**(n+1)+2*cs**(n+1))/4
    assert np.all(ends>=w[:,None]-1e-12)
    assert np.all(ends<=p[:,None]+1e-12)
    f0=.8;low,_,high=models(np.array([f0]),n)
    right.annotate('',xy=(f0,high[0]),xytext=(f0,low[0]),arrowprops=dict(arrowstyle='<->',color='#38566d',lw=1.3))
    right.text(f0-.012,(low[0]+high[0])/2,r'$W(F_{\rm init})$',ha='right',va='center',fontsize=7,color='#38566d')
    right.axhline(.5,color='0.55',ls=':',lw=.9)
    right.set(xlim=(.5,1),ylim=(.25,1),xlabel=r'Initial fidelity $F_{\rm init}$',ylabel=r'End-to-end fidelity $F_{\rm end}$',title=r'$N=3$, $\kappa_N=1$, $\lambda=0.995$')
    right.grid(alpha=.15);right.legend(loc='upper left',framealpha=.95)
    for label,ax in zip(('(a)','(b)'),(left,right)):
        ax.text(-.13,1.06,label,transform=ax.transAxes,fontsize=10,fontweight='bold',va='bottom')
    stem='状态空间与预测范围_ab'
    fig.savefig(output/f'{stem}.png',dpi=220)
    fig.savefig(output/f'{stem}.svg')
    fig.savefig(output/'state_space_envelope.pdf')
    plt.close(fig)
    with (output/f'{stem}_envelope.csv').open('w',newline='',encoding='utf-8') as s:
        writer=csv.writer(s)
        writer.writerow(['F_init','N','kappa','lambda','F_min_Werner','F_double_only','F_max_dephasing','width'])
        writer.writerows((f[i],n,1,quality,w[i],t[i],p[i],p[i]-w[i]) for i in range(len(f)))
    print(output/f'{stem}.png')
    print('Physical-state envelope checks passed')


if __name__=='__main__':main()
