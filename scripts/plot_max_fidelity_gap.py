"""Width of the allowed symmetric-X prediction envelope at fixed initial F."""
from pathlib import Path
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from plot_fixed_delta_comparison import models


def main():
    out=Path(__file__).resolve().parents[1]/'results/state_space_exploration/prediction_envelope_width'
    out.mkdir(parents=True,exist_ok=True)
    f=np.linspace(.5,1,2001)
    quality=.995
    plt.rcParams.update({'font.size':8,'axes.labelsize':9,'legend.fontsize':8})
    fig,ax=plt.subplots(figsize=(3.5,2.8),layout='constrained')
    rows=[]
    largest=0.
    for n,color in zip((1,3,5,7),('#0072B2','#009E73','#D55E00','#8E44AD')):
        lower,_,upper=models(f,n)
        width=upper-lower
        assert np.all(width>=-1e-12) and abs(width[-1])<1e-12
        largest=max(largest,float(width.max()))
        positive=width>0
        ax.plot(f[positive],width[positive],color=color,lw=2.2,label=f'N = {n}')
        rows.extend((f[i],n,1,quality,lower[i],upper[i],width[i]) for i in range(len(f)))
    ax.set(xlim=(.7,1),ylim=(1e-3,.3),yscale='log',
        xlabel=r'Initial fidelity $F_{\rm init}$',
        ylabel=r'$\Delta F_{\max}$')
    from matplotlib.ticker import LogLocator, LogFormatterMathtext, NullFormatter
    ax.yaxis.set_major_locator(LogLocator(base=10,subs=(1.,)))
    ax.yaxis.set_major_formatter(LogFormatterMathtext(base=10))
    ax.yaxis.set_minor_locator(LogLocator(base=10,subs=(2.,3.,4.,5.,6.,7.,8.,9.)))
    ax.yaxis.set_minor_formatter(NullFormatter())
    ax.axhline(.01,color='0.45',ls='--',lw=1,zorder=0)
    ax.grid(which='major',alpha=.18)
    ax.legend(loc='upper right',framealpha=.95)
    fig.savefig(out/'maximum_gap_N1_3_5_7_log_focused_070_100.png',dpi=220)
    fig.savefig(out/'maximum_gap_N1_3_5_7_log_focused_070_100.svg')
    fig.savefig(out/'maximum_gap_summary.pdf')
    plt.close(fig)
    with (out/'maximum_gap_N1_3_5_7.csv').open('w',newline='',encoding='utf-8') as stream:
        writer=csv.writer(stream)
        writer.writerow(['F_init','N','kappa','lambda','F_min','F_max','maximum_gap'])
        writer.writerows(rows)


if __name__=='__main__':main()
