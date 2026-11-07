"""Single zoomed orange accessible-region panel with manuscript contours."""
from pathlib import Path
import csv
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def main():
    out=Path(__file__).resolve().parents[1]/'results/state_space_exploration/requested_component_region'
    with (out/'reachable_boundary_phase_0_05rad.csv').open(encoding='utf-8') as s:
        rows=list(csv.DictReader(s))
    chi=np.array([float(r['chi']) for r in rows])
    delta=np.array([float(r['delta_max_with_phase']) for r in rows])
    plt.rcParams.update({'font.size':8,'axes.labelsize':9,'legend.fontsize':7})
    fig,ax=plt.subplots(figsize=(3.5,3.0),layout='constrained')
    ax.fill_betweenx(chi,-1,delta,color='#f2c88e',alpha=.85,label='Accessible region')
    ax.plot(delta,chi,color='#b66b16',lw=1.5)
    x=np.linspace(.2,1,601);y=np.linspace(.4,1,501)
    X,Y=np.meshgrid(x,y);F=(1+X)*(1+Y)/4
    ctr=ax.contour(X,Y,F,levels=[.55,.65,.75,.85,.95],colors='0.4',linewidths=.85)
    ax.clabel(ctr,fmt='%.2f',fontsize=7)
    threshold=ax.contour(X,Y,F,levels=[.5],colors='#d98b16',linewidths=2)
    ax.clabel(threshold,fmt={.5:r'$F_{\rm init}=0.50$'},fontsize=7,manual=[(.27,.58)])
    ax.plot(x,2*x/(1+x),'--',color='#d75540',lw=2,label='Werner')
    # Use the canonical comparison CSV, not hand-entered marker coordinates.
    cases_path=Path(__file__).resolve().parents[1]/'results/revision/representative_model_cases.csv'
    with cases_path.open(encoding='utf-8') as stream:
        cases={r['case']:r for r in csv.DictReader(stream)}
    for case,label,offset in (
        ('dark_free','I',(-14,-17)),
        ('false_20percent','II',(18,-6)),
        ('near_threshold','III',(-18,-17)),
    ):
        row=cases[case];dx=float(row['delta']);cy=float(row['chi'])
        assert abs(float(row['F_init'])-(1+dx)*(1+cy)/4)<1e-12
        ax.plot(dx,cy,'*',ms=11,mfc='#ffbf00',mec='black',mew=.9,zorder=6,clip_on=False)
        ax.annotate(label,(dx,cy),xytext=offset,textcoords='offset points',
                    ha='center',va='top',fontsize=9,fontweight='bold',color='black',
                    arrowprops=dict(arrowstyle='-',lw=.8,color='black'),zorder=7)
    ax.set(xlim=(.2,1),ylim=(.4,1.005),xlabel=r'Population contrast $\Delta$',ylabel=r'Normalized coherence $\chi=h/b$')
    ax.legend(loc='lower right',fontsize=6.8,framealpha=.95)
    fig.savefig(out/'component_accessible_phase_orange_single.png',dpi=220)
    fig.savefig(out/'component_accessible_phase_orange_single.svg')
    fig.savefig(out/'component_accessible_phase_orange_single.pdf')
    plt.close(fig)


if __name__=='__main__':main()
