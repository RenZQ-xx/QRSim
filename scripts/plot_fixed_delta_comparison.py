"""State-space examples, not fitted SPI hardware scenarios; ideal memory."""
from pathlib import Path
import csv
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np


def models(f, n):
    m, q = n+1, .995**n
    return (.25+3*q*((4*f-1)/3)**m/4,
            .25+q*((2*f-1)**m+2*f**m)/4,
            .25+q*(1+2*(2*f-1)**m)/4)


def save_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main(selected_nodes=(1,2,3), selected_delta=None):
    output=Path(__file__).resolve().parents[1]/"results/state_space_exploration/ideal_memory/fixed_delta_comparison"
    output.mkdir(parents=True, exist_ok=True)
    suffix = "" if selected_delta is None else f"_delta{selected_delta:.2f}"
    summaries=[]
    for n in selected_nodes:
        fig, ax=plt.subplots(figsize=(7.4,5.5),layout="constrained")
        f=np.linspace(.5,1,2001)
        w,t,p=models(f,n)
        ax.plot(f,w,"--",color="#e76f51",lw=2,label="Werner")
        ax.plot(f,t,"-",color="black",lw=2,label="Double-excitation only")
        ax.plot(f,p,"-",color="#56b4e9",lw=2,label="Dephasing only")
        rows=[]
        for d,color in ((.85,"#009E73"),(.9,"#8e44ad"),(.95,"#b8860b")):
            if selected_delta is not None and not np.isclose(d,selected_delta):
                continue
            f=np.linspace(.5,(1+d)/2,5001)
            c=(4*f-1-d)/2
            x=.25+.995**n*(d**(n+1)+2*c**(n+1))/4
            # A representative normalized PSD X state exists at every point:
            # a=d_population=(1-Delta)/4, b=(1+Delta)/4, h=C/2.
            assert np.all(c>=0) and np.all(c<=(1+d)/2+1e-12)
            w,t,p=models(f,n)
            gap=np.minimum.reduce([abs(x-w),abs(x-t),abs(x-p)])
            useful=x>=.5
            j=np.argmax(np.where(useful,gap,-1))
            summaries.append(dict(N=n,delta=d,F_init=f[j],F_end_X=x[j],F_end_Werner=w[j],F_end_double_only=t[j],F_end_phase_only=p[j],minimum_gap_to_three=gap[j],region="maximum nearest-model gap with F_end_X>=0.5",kappa=1,swapping_quality=.995))
            for i in range(len(f)):
                rows.append(dict(N=n,delta=d,F_init=f[i],C=c[i],chi=2*c[i]/(1+d),F_end_X=x[i],F_end_Werner=w[i],F_end_double_only=t[i],F_end_phase_only=p[i],minimum_gap_to_three=gap[i],kappa=1,swapping_quality=.995))
            ax.plot(f,x,"-.",color=color,lw=2,label=rf"X state: $\Delta={d:.2f}$")
            ax.plot(f[-1],x[-1],"o",color=color,ms=4)
        ax.axhline(.5,color="0.6",ls=":",lw=1)
        ax.set(xlim=(.5,1),ylim=(.25,1),xlabel=r"Initial fidelity $F_{\rm init}$",ylabel=r"End-to-end fidelity $F_{\rm end}$",title=rf"$N={n}$, $\kappa_N=1$, $\lambda=0.995$")
        ax.grid(alpha=.18)
        ax.legend(loc="upper left",fontsize=9)
        fig.savefig(output/f"fixed_delta_N{n}{suffix}.png",dpi=220)
        fig.savefig(output/f"fixed_delta_N{n}{suffix}.pdf")
        plt.close(fig)
        save_csv(output/f"fixed_delta_N{n}{suffix}.csv",rows)
        if n == 3 and selected_delta is not None and np.isclose(selected_delta,.95):
            import shutil
            selected=output.parents[1]/"拟在正文中采用的图"
            selected.mkdir(parents=True,exist_ok=True)
            for extension in ("png","pdf","csv"):
                shutil.copy2(output/f"fixed_delta_N{n}{suffix}.{extension}",selected/f"图2_Finit_Fend_N3_delta095.{extension}")
    save_csv(output/f"representative_gap_summary{suffix}.csv",summaries)
    for row in summaries: print(row)


if __name__=="__main__":
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument("--nodes",type=int,nargs="+",default=[1,2,3])
    parser.add_argument("--delta",type=float,choices=[.85,.9,.95])
    args=parser.parse_args()
    main(args.nodes,args.delta)
