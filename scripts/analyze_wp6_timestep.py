"""Compare two strict M0/M1 protocols; never mix cross-protocol legs."""
import csv,hashlib,json
from pathlib import Path
import numpy as np
from wp6_timestep_bounds import pair_tail_bound,sensitivity_interval
from analyze_wp6 import onset
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'analysis/wp6_timestep_20261103'
BASE=ROOT/'analysis/wp6_execution_20261103';REVIEW=ROOT/'analysis/wp6_review_20261103'
def rows(p):
    with p.open() as f:return list(csv.DictReader(f))
def load(p):return json.loads(p.read_text())
def write(name,v):(OUT/name).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def csvout(name,rr):
    with (OUT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rr[0]));w.writeheader();w.writerows(rr)
def nums(rr,k):return np.array([float(r[k]) for r in rr])
def rel(a,b):return abs(b-a)/abs(a)
def widths(x):return np.diff(np.r_[x[0],(x[1:]+x[:-1])/2,x[-1]])

def main():
    dirs={'coarse':{'M0':BASE/'runs/WP6-M0-extend','M1':BASE/'runs/WP6-M1-extend'},'fine':{'M0':REVIEW/'runs/WP6-M0-dt_half','M1':OUT/'runs/WP6-M1-dt_half'}}
    events={};budgets={};weighted={};perf=[]
    for protocol,legs in dirs.items():
        events[protocol]={};budgets[protocol]={};weighted[protocol]={}
        for role,d in legs.items():
            events[protocol][role]=rows(d/'surface_events.csv');budgets[protocol][role]=rows(d/'transport_budget.csv');weighted[protocol][role]=rows(d/'post/weighted_diagnostics.csv')
            f=load(d/'finished.json');s=f['solver'];assert load(d/'post/independent_checks.json')['passed']
            perf.append(dict(protocol=protocol,role=role,accepted_end_s=s['stop']['last_accepted_time_s'],steps=s['stop']['accepted_steps'],wall_s=s['wall_s'],cpu_s=s['cpu_s'],peak_rss_kib=s['peak_rss_kib'],archive_bytes=f['archive']['bytes']))
        assert np.array_equal(nums(events[protocol]['M0'],'time_s'),nums(events[protocol]['M1'],'time_s')),'strict actual pair times differ'
    pair={'coarse':rows(REVIEW/'paired_transport.csv'),'fine':rows(OUT/'refined_pair/paired_transport.csv')}
    cuts=load(dirs['coarse']['M0']/'transport_metadata.json')['cuts'][:7]
    t0=nums(events['coarse']['M0'],'time_s');t1=nums(events['fine']['M0'],'time_s')
    # Pair protocols have independent accepted clocks. Verify the known nested
    # nominal step schedule; never equate the clocks or interpolate integrals.
    assert len(t1)==2*len(t0)-1
    # The coarse last step is shortened by 4 ns to hit TEnd. Record actual
    # clocks throughout; this schedule check is not the time matching rule.
    assert all(abs(float(r['dt_s'])-227.3717729539139)<1e-7 for r in events['coarse']['M0'][1:])
    assert all(abs(float(r['dt_s'])-113.68588647695695)<1e-10 for r in events['fine']['M0'][1:])
    exact=[i for i in range(len(t0)) if t0[i]==t1[2*i]];last_exact=exact[-1]
    records=[];summary=[]
    for cut in cuts:
        reg=cut['name'];pp={k:[r for r in v if r['region']==reg] for k,v in pair.items()}
        bb={k:{role:[r for r in v if r['region']==reg] for role,v in roles.items()} for k,roles in budgets.items()}
        rr=[]
        for i,tc in enumerate(t0):
            j=2*i;tf=t1[j];report=min(tc,tf);bound=0.;detail=None
            if tc!=tf:
                later='coarse' if tc>tf else 'fine';n=i if later=='coarse' else j
                detail=pair_tail_bound(bb[later]['M0'][n-1],bb[later]['M0'][n],bb[later]['M1'][n-1],bb[later]['M1'][n],report);bound=detail['bound_m']
            dc0=float(pp['coarse'][i]['delta_C_m']);dc1=float(pp['fine'][j]['delta_C_m']);iv=sensitivity_interval(dc0,dc1,bound)
            r=dict(region=reg,coarse_step=i,fine_step=j,coarse_accepted_time_s=float(tc),fine_accepted_time_s=float(tf),fine_minus_coarse_time_s=float(tf-tc),common_report_time_for_bound_s=float(report),exact_shared_time=bool(tc==tf),delta_C_coarse_m=dc0,delta_C_fine_m=dc1,**iv,
                relative_U_dt_upper_to_fine_signal=iv['U_dt_upper_m']/abs(dc1) if dc1 else None,
                coarse_residual_envelope_m=float(pp['coarse'][i]['residual_envelope_m']),fine_residual_envelope_m=float(pp['fine'][j]['residual_envelope_m']),
                coarse_corrected_inventory_difference_m=float(pp['coarse'][i]['corrected_inventory_difference_m']),fine_corrected_inventory_difference_m=float(pp['fine'][j]['corrected_inventory_difference_m']))
            rr.append(r);records.append(r)
        end=rr[-1];signal=abs(end['delta_C_fine_m']);dominant=[]
        for r in rr[1:]:
            envelope=max(r['coarse_residual_envelope_m'],r['fine_residual_envelope_m'],r['timing_bound_m'])
            if min(abs(r['delta_C_coarse_m']),abs(r['delta_C_fine_m']))>envelope:dominant.append(r)
        summary.append(dict(region=reg,depth_mm=(1-cut['actual_z_over_H'])*200,**{k:v for k,v in end.items() if k!='region'},
            latest_exact_time_s=float(t0[last_exact]),exact_time_delta_C_coarse_m=rr[last_exact]['delta_C_coarse_m'],exact_time_delta_C_fine_m=rr[last_exact]['delta_C_fine_m'],exact_time_U_dt_m=rr[last_exact]['U_dt_raw_m'],
            final_fine_to_coarse_ratio=end['delta_C_fine_m']/end['delta_C_coarse_m'] if end['delta_C_coarse_m'] else None,
            same_sign_fraction_all_postinitial=sum(r['delta_C_coarse_m']*r['delta_C_fine_m']>0 for r in rr[1:])/len(rr[1:]),
            resolved_sample_count=len(dominant),same_sign_fraction_above_both_closure=sum(r['delta_C_coarse_m']*r['delta_C_fine_m']>0 for r in dominant)/len(dominant) if dominant else None,
            fine_signal_over_Udt_upper=signal/end['U_dt_upper_m'] if end['U_dt_upper_m'] else None,
            paired_convergence_status='single refinement sensitivity only; no convergence order established'))
    csvout('paired_timestep_sensitivity.csv',records);csvout('seven_plane_timestep_summary.csv',summary);csvout('four_leg_performance.csv',perf)
    write('clock_comparison.json',dict(exact_shared_report_times_s=[float(t0[i]) for i in exact],latest_exact_report_time_s=float(t0[last_exact]),coarse_end_s=float(t0[-1]),fine_end_s=float(t1[-1]),end_gap_s=float(t1[-1]-t0[-1]),max_nominal_report_gap_s=float(max(abs(t1[::2]-t0))),
        coarse_final_dt_s=float(events['coarse']['M0'][-1]['dt_s']),fine_final_dt_s=float(events['fine']['M0'][-1]['dt_s']),interpretation='U_dt at exactly shared clocks is direct. At unequal endpoints, report raw difference and conservative accepted-step interval; no interpolated cumulative values created.',
        bound_formula='B = sum over later M0,M1 of |gap| (|q_adv|+|q_diff|) + |Delta C_step - q_total Delta t_actual| + 8 max(ulp(C_current),ulp(C_previous))',
        bound_scope='recorded backward-Euler conservative quadrature only; not a bound on continuous unresolved flux or U_z/U_x/U_domain'))
    growth={};checks={};compare=[]
    for protocol in dirs:
        rr=weighted[protocol]['M1'];sigma=nums(rr,'weighted_top_sigma');m0sigma=nums(weighted[protocol]['M0'],'weighted_top_sigma');detect=onset(dirs[protocol]['M1'])
        growth[protocol]=dict(accepted_end_s=float(events[protocol]['M1'][-1]['time_s']),weighted_sigma_end=float(sigma[-1]),weighted_sigma_min=float(min(sigma)),weighted_amplification=float(sigma[-1]/min(sigma)),fourier_end=float(rr[-1]['weighted_fourier']),surface_mean=float(rr[-1]['weighted_top_mean']),surface_max=float(rr[-1]['top_max']),domain_max=float(rr[-1]['domain_max']),detection=detect,
            max_M0_sigma_over_M1_sigma=float(max(m0sigma/sigma)),max_M1_odd_energy=float(max(nums(rr,'odd_neumann_energy_fraction'))))
        checks[protocol+'_M1_symmetry']=growth[protocol]['max_M1_odd_energy']<1e-6 and all(float(v['max_endpoint_difference'])<max(1e-12,1e-3*float(v['weighted_fourier'])) for v in rr)
        checks[protocol+'_background_signal_separation']=growth[protocol]['max_M0_sigma_over_M1_sigma']<=.001
    for key in ['weighted_sigma_end','weighted_amplification','fourier_end','surface_mean','surface_max','domain_max']:
        compare.append(dict(quantity=key,coarse=growth['coarse'][key],fine=growth['fine'][key],relative_difference=rel(growth['coarse'][key],growth['fine'][key]),coarse_time_s=float(t0[-1]),fine_time_s=float(t1[-1])))
    for k in ['weighted_sigma_end','weighted_amplification','fourier_end']:checks[k+'_5percent']=rel(growth['coarse'][k],growth['fine'][k])<=.05
    for k in ['surface_mean','surface_max','domain_max']:checks[k+'_0p1percent']=rel(growth['coarse'][k],growth['fine'][k])<=.001
    ta=growth['coarse']['detection']['onset_refined_s'];tb=growth['fine']['detection']['onset_refined_s'];checks['detection_consistency']=abs(tb-ta)<=max(2*227.3717729539139,.02*ta)
    coarse_cache=ROOT/'tmp/wp6_execution_20261103/runs/WP6-M1-extend/post/profiles.npz';fine_cache=ROOT/'tmp/wp6_timestep_20261103/runs/WP6-M1-dt_half/post/profiles.npz'
    cp=np.load(coarse_cache);fp=np.load(fine_cache);assert np.array_equal(cp['z_m'],fp['z_m']);wz=widths(cp['z_m']);profiles=[];profile_checks={}
    for label,i,j in [('latest_exact',last_exact,2*last_exact),('own_endpoints',len(t0)-1,len(t1)-1)]:
        a=cp['x_sigma'][i];b=fp['x_sigma'][j];error=float(np.sqrt(sum((b-a)**2*wz)/sum(a*a*wz)))
        profile_checks[label]=dict(relative_sigma_profile_L2_error=error,coarse_accepted_time_s=float(t0[i]),fine_accepted_time_s=float(t1[j]))
        checks[label+'_vertical_profile_5percent']=error<=.05
        for z,x,y in zip(cp['z_m'],a,b):profiles.append(dict(report=label,z_m=float(z),coarse_time_s=float(t0[i]),fine_time_s=float(t1[j]),coarse_sigma=float(x),fine_sigma=float(y)))
    csvout('vertical_profile_comparison.csv',profiles);csvout('growth_endpoint_comparison.csv',compare)
    # Distinct engineering exploration and physical-claim decisions.
    write('growth_and_detection.json',dict(protocols=growth,profile_comparison=profile_checks,exact_common_report_time_s=float(t0[last_exact]),endpoint_time_gap_s=float(t1[-1]-t0[-1])))
    pairinfo=load(OUT/'refined_pair/async_stopping_checks.json')
    write('decision_checks.json',dict(exploratory_A='GO_for_proposal_only' if all(checks.values()) else 'NO_GO',exploratory_checks={k:bool(v) for k,v in checks.items()},
        quantitative_B='NO_GO',quantitative_reason='U_z/U_x/U_domain unknown; no observed solubility events; a single timestep refinement cannot validate tiny transport or delay',
        original_M0_absolute_screen='failed_retained',alternative_criterion_registration='20261103-C01, before new M1 fine result, engineering scope only',
        refined_surface_event_difference=pairinfo['pair']['surface_time_difference'],unknown_errors=['U_z','U_x','U_domain'],t_eff=None,four_T_started=False))
    print(json.dumps(dict(decision=load(OUT/'decision_checks.json'),growth=growth,planes=[{k:r[k] for k in ['depth_mm','delta_C_coarse_m','delta_C_fine_m','U_dt_raw_m','U_dt_upper_m','timing_bound_m','final_fine_to_coarse_ratio']} for r in summary]),indent=2))
if __name__=='__main__':main()
