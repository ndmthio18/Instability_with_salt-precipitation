"""Exploratory 4T acceptance/growth analysis; immutable histories remain read-only."""
import csv,hashlib,json
from pathlib import Path
import numpy as np
from analyze_wp6 import onset
from wp6_4t_metrics import log_growth_fit
from wp6_paired_stopping import same_time
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'analysis/wp6_4T_exploratory_20261103';T=136423.0637723484
def rows(p):
    with p.open() as f:return list(csv.DictReader(f))
def load(p):return json.loads(p.read_text())
def nums(r,k):return np.array([float(x[k]) for x in r])
def write(n,v):(OUT/n).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def csvout(n,r):
    with (OUT/n).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(r[0]));w.writeheader();w.writerows(r)
def main():
    perf=[];metrics={};check={};weighted={};events={};checkpoints=[];bridges={}
    old={'M0':ROOT/'analysis/wp6_review_20261103/runs/WP6-M0-dt_half','M1':ROOT/'analysis/wp6_timestep_20261103/runs/WP6-M1-dt_half'}
    for role in ['M0','M1']:
        d=OUT/'runs'/('WP6-'+role+'-explore4T');f=load(d/'finished.json');s=f['solver'];w=rows(d/'post/weighted_diagnostics.csv');ev=rows(d/'surface_events.csv');r=rows(d/'transport_budget.csv');weighted[role]=w;events[role]=ev
        au=load(d/'post/field_audit.json');check[role+'_independent_quadrature']=load(d/'post/independent_checks.json')['passed'];check[role+'_initial_fields']=au['all_rho_pass'] and au['all_sw_pass'] and au['initial_field_error']<1e-12
        for saltwater,tol in [('salt',1e-14),('water',1e-11)]:
            for whole in ['', 'whole_']:check[role+'_'+whole+saltwater+'_closure']=max(abs(nums(r,whole+saltwater+'_residual_m')))<tol
        times=nums(ev,'time_s');sigma=nums(w,'weighted_top_sigma');amp=nums(w,'weighted_fourier')
        metrics[role]=dict(end_actual_s=float(times[-1]),accepted_steps=s['stop']['accepted_steps'],stop=s['stop'],surface_mean=float(w[-1]['weighted_top_mean']),surface_max=float(w[-1]['top_max']),domain_max=float(w[-1]['domain_max']),max_sigma=float(max(sigma)),max_endpoint=float(max(nums(w,'max_endpoint_difference'))),max_mirror=float(max(nums(w,'max_mirror_error'))),max_salt_residual_m=float(max(abs(nums(r,'salt_residual_m')))),max_water_residual_m=float(max(abs(nums(r,'water_residual_m')))))
        check[role+'_original_absolute_endpoint']=metrics[role]['max_endpoint']<=1e-12
        if role=='M1':
            metrics[role].update(weighted_amplification=float(sigma[-1]/min(sigma)),sigma_min=float(min(sigma)),sigma_end=float(sigma[-1]),fourier_end=float(amp[-1]),projection_end=float(w[-1]['weighted_projection']),relative_sigma_end=float(w[-1]['weighted_top_relative']),detection=onset(d),max_odd_energy=float(max(nums(w,'odd_neumann_energy_fraction'))))
            metrics[role]['amplitude_extrema']={}
            for key in ['weighted_top_sigma','weighted_fourier','weighted_projection']:
                values=abs(nums(w,key));i=int(np.argmax(values));j=int(np.argmin(values))
                metrics[role]['amplitude_extrema'][key]=dict(minimum=float(values[j]),minimum_actual_s=float(times[j]),maximum=float(values[i]),maximum_actual_s=float(times[i]),final_over_minimum=float(values[-1]/values[j]),final_over_initial=float(values[-1]/values[0]),final_over_peak=float(values[-1]/values[i]))
            check['M1_original_relative_symmetry']=metrics[role]['max_odd_energy']<1e-6 and all(float(v['max_endpoint_difference'])<max(1e-12,1e-3*float(v['weighted_fourier'])) for v in w)
        for factor in [0,1,2,3,4]:
            if factor*T>times[-1]+1e-7:continue
            i=int(np.argmin(abs(times-factor*T)));v=w[i]
            checkpoints.append(dict(role=role,nominal_T=factor,actual_time_s=float(times[i]),time_minus_nominal_s=float(times[i]-factor*T),sigma=float(v['weighted_top_sigma']),fourier=float(v['weighted_fourier']),projection=float(v['weighted_projection']),surface_mean=float(v['weighted_top_mean']),surface_max=float(v['top_max']),domain_max=float(v['domain_max']),relative_sigma=float(v['weighted_top_relative']),amplification_from_global_min=float(sigma[i]/min(sigma)) if min(sigma)>0 else None))
        perf.append(dict(role=role,steps=s['stop']['accepted_steps'],wall_s=s['wall_s'],cpu_s=s['cpu_s'],rss_kib=s['peak_rss_kib'],post_archive_wall_s=f['post_wall_s'],raw_bytes=s['solver_output_bytes'],archive_bytes=f['archive']['bytes']))
        # Exact accepted-time overlap and field hashes, not nominal T labels.
        before=rows(old[role]/'surface_events.csv');oldw=rows(old[role]/'post/weighted_diagnostics.csv');oldhash=load(old[role]/'post/field_hashes.json');newhash=load(d/'post/field_hashes.json')
        idx={float(x['time_s']):i for i,x in enumerate(ev)};pairs=[(i,idx[float(x['time_s'])]) for i,x in enumerate(before) if float(x['time_s']) in idx]
        errors={k:max(abs(float(oldw[i][k])-float(w[j][k])) for i,j in pairs) for k in ['weighted_top_mean','weighted_top_sigma','weighted_fourier','weighted_projection','whole_salt_m']}
        bridges[role]=dict(reference=str(old[role].relative_to(ROOT)),exact_shared_states=len(pairs),reference_states=len(before),identical_VTU_hashes=sum(oldhash[i]['sha256']==newhash[j]['sha256'] for i,j in pairs),maximum_absolute_differences=errors,latest_exact_overlap_s=float(ev[pairs[-1][1]]['time_s']))
    check={k:bool(v) for k,v in check.items()};core=all(v for k,v in check.items() if 'original_absolute_endpoint' not in k and k!='M1_original_relative_symmetry')
    check['core_numerical_acceptance']=core
    # compare noise only at genuinely shared accepted times from the validated interface
    paired=load(OUT/'paired/async_stopping_checks.json');ratios=[float(weighted['M0'][i]['weighted_top_sigma'])/float(weighted['M1'][j]['weighted_top_sigma']) for i,j in paired['pair']['indices']]
    check['engineering_noise_signal_margin']=max(ratios,default=float('inf'))<1e-3
    def alert_interval(values,times,threshold):
        hit=next((i for i,v in enumerate(values) if v>threshold),None)
        return None if hit is None else dict(low_s=float(times[max(0,hit-1)]),high_s=float(times[hit]),threshold=threshold,meaning='sampled diagnostic alert, not a physical onset label')
    shared_times=[float(events['M0'][i]['time_s']) for i,j in paired['pair']['indices']]
    alerts=dict(background_ratio=alert_interval(ratios,shared_times,1e-3),M1_odd_energy=alert_interval(nums(weighted['M1'],'odd_neumann_energy_fraction'),nums(events['M1'],'time_s'),1e-6),M0_absolute_endpoint=alert_interval(nums(weighted['M0'],'max_endpoint_difference'),nums(events['M0'],'time_s'),1e-12))
    write('numerical_acceptance.json',dict(checks=check,metrics=metrics,max_M0_sigma_over_M1_sigma=max(ratios,default=None),diagnostic_alert_intervals=alerts,historical_M0_absolute_screen='failed_retained',errors=dict(U_dt_4T='unknown; dt/2 verified only through 2T',U_z='unknown',U_x='unknown',U_domain='unknown'),t_eff=None))
    csvout('runtime_cost.csv',perf);csvout('T_2T_4T_checkpoints.csv',checkpoints);write('fine_2T_overlap.json',bridges)
    fits=[];times=nums(events['M1'],'time_s')
    for k in ['weighted_top_sigma','weighted_fourier','weighted_projection']:
        values=abs(nums(weighted['M1'],k))
        for a,b in [(1,2),(2,3),(3,4)]:fits.append(dict(quantity=k,window_start_T=a,window_end_T=b,**log_growth_fit(times,values,a*T,b*T)))
    write('growth_rate_fits.json',fits)
    # Depth-resolved transport evolution on actual common states.
    pairs=rows(OUT/'paired/paired_transport.csv');summ=[]
    for reg in [c['name'] for c in load(OUT/'runs/WP6-M0-explore4T/transport_metadata.json')['cuts'][:7]]:
        rr=[r for r in pairs if r['region']==reg];t=nums(rr,'time_s');i2=int(np.argmin(abs(t-2*T)));end=rr[-1]
        summ.append(dict(region=reg,paired_end_s=float(t[-1]),actual_near_2T_s=float(t[i2]),delta_C_near_2T_m=float(rr[i2]['delta_C_m']),delta_C_end_m=float(end['delta_C_m']),end_to_2T_ratio=float(end['delta_C_m'])/float(rr[i2]['delta_C_m']) if float(rr[i2]['delta_C_m']) else None,
            relative_to_initial_inventory=float(end['relative_delta_C_initial_inventory']),corrected_inventory_difference_m=float(end['corrected_inventory_difference_m']),residual_envelope_end_m=float(end['residual_envelope_m']),
            **{k:float(end[k]) for k in end if k.startswith('M0_signed') or k.startswith('M1_signed')},
            M0_net_downward_sample_fraction=float(np.mean(nums(rr,'M0_signed_cut_total')>0)),M1_net_downward_sample_fraction=float(np.mean(nums(rr,'M1_signed_cut_total')>0))))
    csvout('seven_plane_2T_to_4T.csv',summ)
    print(json.dumps(dict(acceptance=load(OUT/'numerical_acceptance.json'),bridge=bridges,growth_fits=fits,transport=summ),indent=2))
if __name__=='__main__':main()
