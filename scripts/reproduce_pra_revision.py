"""Canonical, seeded regeneration of the final PRA revision results.

Run from the repository root. Historical outputs and main.tex are untouched.
"""
from pathlib import Path
import csv
import math
import shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
from QRSim.analytical import homogeneous_chain, symmetric_elementary_link, reduced_model_fidelity
from QRSim.MonteCarloSim import MonteCarloSim
from QRSim.waiting_time import accumulated_storage_rounds
from reproduce_figures import MANUSCRIPT_PARAMETERS, reproduce_component_sensitivity
from plot_state_space_envelope_panels import main as state_space_figure
from map_requested_component_region import main as map_region
from map_component_region_with_phase import main as map_phase_region
from plot_component_region_single import main as plot_region
from plot_max_fidelity_gap import main as plot_maximum_gap
from write_revision_tables import main as write_tables

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/revision'


def save(name,rows):
    with (OUT/name).open('w',newline='',encoding='utf-8') as s:
        writer=csv.DictWriter(s,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)


def compare(parameters):
    result=homogeneous_chain(**parameters)
    link=result.elementary_link;n=parameters['repeater_nodes']
    C=2*link.coherence;w=(4*link.fidelity-1)/3
    models={name:reduced_model_fidelity(link.fidelity,name,repeater_nodes=n,coherence_factor=result.coherence_factor,swapping_quality=parameters['swapping_quality']) for name in ('werner','phase_only','bit_only')}
    row=dict(**parameters,F_init=link.fidelity,delta=link.delta,C=C,chi=link.coherence/link.population_single_excitation,
        kappa=result.coherence_factor,F_X=result.fidelity,F_W=models['werner'],
        F_phase=models['phase_only'],F_bit=models['bit_only'],delta_F=models['werner']-result.fidelity,
        rate_hz=result.rate,p_hr=link.heralding_probability)
    row['S_delta']=parameters['swapping_quality']**n*(n+1)/4*(link.delta**n-result.coherence_factor*C**n)
    row['delta_minus_w']=link.delta-w
    row['trace_distance_twirl_to_werner']=abs(link.delta-w)/2
    eta=parameters['transmission_without_fiber']*10**(-parameters['attenuation_db_per_km']*parameters['total_distance_km']/(n+1)/20)
    pe=parameters['excitation_probability'];dc=parameters['dark_count_probability']
    row['false_fraction']=2*dc*(1-dc)*(1-pe*eta)**2/link.heralding_probability
    nu=parameters.get('phase_coherence',1.)
    s=pe*eta*parameters['detector_efficiency']
    approx=(s*(1-pe)*(1+math.sqrt(parameters['visibility'])*nu)/2+dc*pe*(1-pe))/(s+dc)
    row['F_init_approx']=approx
    row['F_W_approx']=reduced_model_fidelity(approx,'werner',repeater_nodes=n,coherence_factor=result.coherence_factor,swapping_quality=parameters['swapping_quality'])
    row['structure_error']=row['F_W']-row['F_X']
    row['link_estimate_error']=row['F_W_approx']-row['F_W']
    return row,result


def validate_case(parameters,trials=500000,seed=20261105):
    row,result=compare(parameters)
    n=parameters['repeater_nodes'];p=result.elementary_link.heralding_probability
    times=np.random.default_rng(seed).exponential(size=(trials,n+1))
    maximum=times.max(axis=1)
    storage=2*maximum-times[:,0]-times[:,-1]+np.abs(np.diff(times,axis=1)).sum(axis=1)
    continuum=np.exp(-storage/(p*result.repetition_rate*parameters['coherence_time']))
    rounds=np.floor(times/-math.log1p(-p)).astype(np.int64)+1
    finite=np.exp(-accumulated_storage_rounds(rounds)/(result.repetition_rate*parameters['coherence_time']))
    paired=continuum-finite
    lmd=parameters['swapping_quality']**n;C=row['C'];w=(4*row['F_init']-1)/3
    return dict(N=n,trials=trials,seed=seed,p_hr=p,kappa_continuum_sample=float(continuum.mean()),
        kappa_continuum_se=float(continuum.std(ddof=1)/math.sqrt(trials)),
        continuum_minus_finite=float(paired.mean()),paired_se=float(paired.std(ddof=1)/math.sqrt(trials)),
        continuum_bias_F_X=float(paired.mean()*lmd*C**(n+1)/2),
        continuum_bias_model_gap=float(paired.mean()*lmd*(w**(n+1)-C**(n+1))/2),
        gap_mc_se=float(continuum.std(ddof=1)/math.sqrt(trials)*lmd*abs(w**(n+1)-C**(n+1))/2))


def main():
    OUT.mkdir(parents=True,exist_ok=True)
    base=dict(MANUSCRIPT_PARAMETERS)
    distance=[];excitation=[]
    for n in (1,3,5):
        for length in np.linspace(100,1000,100):
            row,_=compare(dict(**base,repeater_nodes=n,total_distance_km=float(length),excitation_probability=.03));distance.append(row)
    for n in range(3,9):
        for pe in np.linspace(.001,.12,120):
            row,_=compare(dict(**base,repeater_nodes=n,total_distance_km=1000.,excitation_probability=float(pe)));excitation.append(row)
    save('baseline_distance_comparison.csv',distance);save('baseline_excitation_comparison.csv',excitation)
    for label,data in [('distance',distance),('excitation',excitation)]:
        for threshold in (0.,.5):
            selected=[r for r in data if r['F_X']>=threshold]
            peak=max(selected,key=lambda r:abs(r['delta_F']))
            print(label,'threshold',threshold,'maximum gap',peak['delta_F'],'N',peak['repeater_nodes'],'L',peak['total_distance_km'],'pe',peak['excitation_probability'])
    for n in range(3,9):
        peak=max((r for r in excitation if r['repeater_nodes']==n),key=lambda r:r['F_X'])
        print('PEAK',n,peak['excitation_probability'],peak['F_X'],peak['rate_hz'])
    cases=[];validation=[]
    paper=dict(**base,repeater_nodes=6,total_distance_km=1000.,excitation_probability=.036)
    row,_=compare(paper);cases.append(dict(case='paper_baseline',**row,sigma_phi_rad=.30))
    # Cases I--III share hardware; I->II varies dark counts and II->III phase noise.
    hardware={**base,'repeater_nodes':2,'total_distance_km':75.,'excitation_probability':.001,
        'transmission_without_fiber':.019,'visibility':.999,'coherence_time':600.,'phase_coherence':math.exp(-.10**2/2)}
    eta=.019*10**(-.2*12.5/10)
    lo,hi=0.,1e-5
    for _ in range(70):
        dc=(lo+hi)/2
        link=symmetric_elementary_link(excitation_probability=.001,transmission=eta,detector_efficiency=.9,dark_count_probability=dc,visibility=.999)
        frac=2*dc*(1-dc)*(1-.001*eta)**2/link.heralding_probability
        if frac<.2:lo=dc
        else:hi=dc
    target_dc=(lo+hi)/2
    dark_rows=[]
    for factor in np.linspace(0,1.5,151):
        p={**hardware,'dark_count_probability':target_dc*factor}
        row,_=compare(p);dark_rows.append(dict(dark_count_factor=float(factor),**row))
    for label,factor,sigma in [('dark_free',0.,.10),('false_20percent',1.,.10),('near_threshold',1.,.35)]:
        p={**hardware,'dark_count_probability':target_dc*factor,
           'phase_coherence':math.exp(-sigma**2/2)}
        row,_=compare(p);cases.append(dict(case=label,**row,sigma_phi_rad=sigma))
        validation.append(dict(case=label,**validate_case(p)))
        print('CASE',label,row)
    validation.append(dict(case='paper_baseline',**validate_case(paper)))
    save('representative_model_cases.csv',cases);save('false_heralding_sweep.csv',dark_rows)
    save('new_case_continuum_validation.csv',validation)
    reproduce_component_sensitivity(OUT)
    state_space_figure(OUT)
    # Independent statistical errors of the existing cached continuum average.
    rr,result=compare(paper);data=MonteCarloSim.load_freq_dist(repeater_nodes=6)
    samples=np.exp(-np.asarray(data['W'])/(rr['p_hr']*result.repetition_rate*600))
    se=float(samples.std(ddof=1)/math.sqrt(samples.size))
    print('baseline kappa sampling SE',se,'F sampling SE',se*.995**6*rr['C']**7/2)
    write_tables()
    map_region()
    map_phase_region()
    plot_region()
    plot_maximum_gap()
    exploration=ROOT/'results/state_space_exploration'
    shutil.copy2(exploration/'requested_component_region/component_accessible_phase_orange_single.pdf',OUT/'component_accessible_region.pdf')
    shutil.copy2(exploration/'prediction_envelope_width/maximum_gap_summary.pdf',OUT/'maximum_gap_summary.pdf')
    shutil.copy2(exploration/'prediction_envelope_width/maximum_gap_N1_3_5_7.csv',OUT/'maximum_gap_N1_3_5_7.csv')
    for name in ('reachable_boundary.csv','reachable_boundary_phase_0_05rad.csv'):
        shutil.copy2(exploration/'requested_component_region'/name,OUT/name)
    print('Revision outputs:',OUT)


if __name__=='__main__':main()
