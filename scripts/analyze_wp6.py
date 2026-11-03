"""WP6 acceptance and paired transport analysis. Never launches a simulation."""
import argparse,csv,hashlib,importlib.util,json
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'analysis/wp6_execution_20261103';WP5=ROOT/'analysis/single_wavelength_validation_20261103'
def rows(p):return list(csv.DictReader(Path(p).open()))
def nums(rr,k):return np.array([float(r[k]) for r in rr])
def write(name,v): (OUT/name).write_text(json.dumps(v,indent=2)+'\n')
def csvout(name,rr):
    with (OUT/name).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=list(rr[0]));w.writeheader();w.writerows(rr)
def onset(d):
    spec=importlib.util.spec_from_file_location('frozen',ROOT/'analysis/pilot0_integrated_20261102/provenance/extract_onset_2025.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m.extract_onset(str(next(d.glob('*sigmatop.csv'))))
def fields_bridge():
    from collect_post_onset import remote
    code=r'''import json,hashlib,xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
result={}
for role in ['M0','M1']:
 old=Path('/temp/zhixin/paper2_wp5_20261103_r02/runs')/('R24-'+role+'-full')
 new=Path('/temp/zhixin/paper2_wp6_execution_20261103_r01/runs')/('WP6-'+role+'-bridge')
 expected=json.loads((old/'post/field_hashes.json').read_text());actual=json.loads((new/'post/field_hashes.json').read_text());assert len(expected)==len(actual)==601
 identical=0;errors={};records=[]
 for a,b in zip(expected,actual):
  x=Path(a['path']).read_bytes();y=Path(b['path']).read_bytes();ha=hashlib.sha256(x).hexdigest();hb=hashlib.sha256(y).hexdigest();assert ha==a['sha256'] and hb==b['sha256']
  identical+=ha==hb;records.append(dict(frame=a['frame'],reference_sha=ha,new_sha=hb))
  if ha!=hb:
   aa={v.attrib.get('Name','unnamed_'+str(i)):np.fromstring(v.text or '',sep=' ') for i,v in enumerate(ET.fromstring(x).iter('DataArray'))}
   bb={v.attrib.get('Name','unnamed_'+str(i)):np.fromstring(v.text or '',sep=' ') for i,v in enumerate(ET.fromstring(y).iter('DataArray'))}
   assert aa.keys()==bb.keys()
   for k in aa:assert aa[k].shape==bb[k].shape;errors[k]=max(errors.get(k,0),float(np.max(abs(aa[k]-bb[k]),initial=0)))
 result[role]=dict(frames=601,byte_identical_frames=identical,max_field_absolute_difference=errors,hashes=records)
print(json.dumps(result))
'''
    result=remote(code,timeout=180);write('bridge_field_comparison.json',result);print({k:{a:b for a,b in v.items() if a!='hashes'} for k,v in result.items()})
def main(stage):
    checks={};metrics={};budget={};weighted={};onsets={};perf=[]
    for role in ['M0','M1']:
        d=OUT/'runs'/f'WP6-{role}-{stage}';f=json.loads((d/'finished.json').read_text());s=f['solver'];r=rows(d/'transport_budget.csv');w=rows(d/'post/weighted_diagnostics.csv');budget[role]=r;weighted[role]=w
        checks[role+'/finished']=s['returncode']==0 and s['stop']['reason'] in ['requested_horizon','surface_solubility_threshold','interior_solubility_threshold']
        checks[role+'/independent_quad_and_budget']=json.loads((d/'post/independent_checks.json').read_text())['passed']
        audit=json.loads((d/'post/field_audit.json').read_text())
        checks[role+'/rho_sw_initial_field']=audit['all_rho_pass'] and audit['all_sw_pass'] and audit['initial_field_error']<1e-12
        checks[role+'/whole_budget']=max(abs(nums(r,'whole_salt_residual_m')))<1e-14 and max(abs(nums(r,'whole_water_residual_m')))<1e-11
        locked=json.loads((OUT/'amendment_20261103_A01.json').read_text())['cut_coordinates']
        actual=json.loads((d/'transport_metadata.json').read_text())['cuts'][:7]
        checks[role+'/fixed_cuts']=all(abs(c['actual_z_over_H']*.2-z)<1e-12 for c,z in zip(actual,locked)) and len(actual)==7
        events=rows(d/'surface_events.csv');checks[role+'/accepted_dt']=all(abs(float(v['dt_s'])-227.3717729539139)<1e-7 for v in events[1:])
        checks[role+'/salt_closure']=max(abs(nums(r,'salt_residual_m')))<1e-14;checks[role+'/water_closure']=max(abs(nums(r,'water_residual_m')))<1e-11
        metrics[role]=dict(stop=s['stop'],max_salt_residual_m=max(abs(nums(r,'salt_residual_m'))),max_water_residual_m=max(abs(nums(r,'water_residual_m'))),max_sigma=max(nums(w,'weighted_top_sigma')),final_surface_mean=float(w[-1]['weighted_top_mean']),final_surface_max=float(w[-1]['top_max']),final_domain_max=float(w[-1]['domain_max']),max_endpoint_difference=max(nums(w,'max_endpoint_difference')))
        if role=='M0':checks['M0/absolute_endpoint_screen']=metrics[role]['max_endpoint_difference']<=1e-12
        if role=='M1':
            onsets[role]=onset(d);sigma=nums(w,'weighted_top_sigma');metrics[role]['weighted_growth_since_minimum']=sigma[-1]/min(sigma);metrics[role]['relative_sigma_end']=float(w[-1]['weighted_top_relative']);metrics[role]['odd_energy_max']=max(nums(w,'odd_neumann_energy_fraction'))
            checks[role+'/symmetry']=metrics[role]['odd_energy_max']<1e-6 and all(float(v['max_endpoint_difference'])<max(1e-12,1e-3*float(v['weighted_fourier'])) for v in w)
        perf.append(dict(run=d.name,accepted_steps=s['stop']['accepted_steps'],wall_s=s['wall_s'],cpu_s=s['cpu_s'],peak_rss_kib=s['peak_rss_kib'],newton_iterations=s['newton_iterations'],linear_iterations='',linear_status='unknown/not emitted',solver_output_bytes=s['solver_output_bytes'],archive_bytes=f['archive']['bytes'],post_wall_s=f['post_wall_s']))
        if stage=='bridge':
            old=WP5/'runs'/f'R24-{role}-full';a=rows(old/'post/weighted_diagnostics.csv')
            checks[role+'/accepted_times']=len(a)==len(w) and np.max(abs(nums(a,'time_s')-nums(w,'time_s')))<1e-7
            for key in ['weighted_top_mean','weighted_top_sigma','weighted_fourier','weighted_projection','whole_salt_m','whole_water_m','top_max']:
                delta=float(max(abs(nums(a,key)-nums(w,key))));metrics[role]['bridge_'+key+'_max_abs']=delta
                tol=1e-12 if key in ['weighted_top_sigma','weighted_fourier','weighted_projection'] else max(abs(nums(a,key)))*1e-10
                checks[role+'/'+key]=delta<=tol
            b=rows(old/'transport_budget.csv')
            for region in ['upper_0.900000','upper_0.750000','upper_0.500000']:
                aa=[v for v in b if v['region']==region];bb=[v for v in r if v['region']==region]
                for key in ['cut_salt_total_down_m_s','cut_salt_adv_down_m_s','cut_salt_diff_down_m_s','inventory_salt_m','inventory_water_m','cumulative_cut_salt_m','salt_residual_m','water_residual_m']:
                    delta=float(max(abs(nums(aa,key)-nums(bb,key))));metrics[role][region+'/'+key+'_difference']=delta
                    checks[role+'/'+region+'/'+key]=delta<(1e-18 if 'm_s' in key else 1e-12 if 'water' in key else 1e-15)
            if role=='M1':checks['onset_unchanged']=abs(onsets[role]['onset_refined_s']-onset(old)['onset_refined_s'])<1e-5
    checks['same_accepted_time_series']=len(weighted['M0'])==len(weighted['M1']) and max(abs(nums(weighted['M0'],'time_s')-nums(weighted['M1'],'time_s')))<1e-7
    pairs=[];summary=[]
    for reg in [c['name'] for c in json.loads((OUT/'runs'/f'WP6-M0-{stage}/transport_metadata.json').read_text())['cuts']]:
        a=[r for r in budget['M0'] if r['region']==reg];b=[r for r in budget['M1'] if r['region']==reg];idx0={round(float(r['time_s']),6):r for r in a};idx1={round(float(r['time_s']),6):r for r in b};times=sorted(set(idx0)&set(idx1));init=float(b[0]['inventory_salt_m'])-float(a[0]['inventory_salt_m']);checks[reg+'/initial_pair_inventory']=abs(init)<1e-15;rr=[]
        for t in times:
            x,y=idx0[t],idx1[t];assert abs(float(x['time_s'])-float(y['time_s']))<1e-7
            dc=float(y['cumulative_cut_salt_m'])-float(x['cumulative_cut_salt_m']);di=float(y['inventory_salt_m'])-float(x['inventory_salt_m'])-init;dr=float(y['salt_residual_m'])-float(x['salt_residual_m']);env=abs(float(y['salt_residual_m']))+abs(float(x['salt_residual_m']))
            q=dict(region=reg,time_s=float(x['time_s']),initial_difference_m=init,corrected_inventory_difference_m=di,delta_C_m=dc,relative_delta_C_initial_inventory=dc/float(x['inventory_initial_salt_m']),delta_C_adv_m=float(y['cumulative_cut_adv_m'])-float(x['cumulative_cut_adv_m']),delta_C_diff_m=float(y['cumulative_cut_diff_m'])-float(x['cumulative_cut_diff_m']),M0_signed_cut_total=float(x['cut_salt_total_down_m_s']),M1_signed_cut_total=float(y['cut_salt_total_down_m_s']),residual_envelope_m=env,identity_error_m=di+dc-dr,candidate_threshold_m=max(1e-13,10*env),instantaneous_candidate=abs(dc)>max(1e-13,10*env))
            pairs.append(q);rr.append(q)
        checks[reg+'/pair_identity']=max(abs(v['identity_error_m']) for v in rr)<1e-18
        # Persistent sign/threshold in a physical interval >=10 base steps, not 10 samples after dt changes.
        triggered=False
        for j,q in enumerate(rr):
            start=q['time_s']-10*227.3717729539139
            if start<0:continue
            window=[v for v in rr[:j+1] if v['time_s']>=start-1e-7]
            if window[-1]['time_s']-window[0]['time_s']<10*227.3717729539139-1e-7:continue
            for sign in [-1,1]:
                if sum(v['instantaneous_candidate'] and sign*v['delta_C_m']>0 for v in window)/len(window)>=.8:triggered=True
        summary.append(dict(region=reg,initial_difference_m=init,delta_C_end_m=rr[-1]['delta_C_m'],max_abs_delta_C_m=max(abs(v['delta_C_m']) for v in rr),relative_delta_C_end=rr[-1]['relative_delta_C_initial_inventory'],max_residual_envelope_m=max(v['residual_envelope_m'] for v in rr),persistent_candidate=triggered,classification='candidate_requires_sensitivity' if triggered else 'unresolved'))
    if stage=='bridge':
        fields=json.loads((OUT/'bridge_field_comparison.json').read_text())
        for role,r in fields.items():checks[role+'/actual_2d_fields']=r['frames']==601 and all(v<=1e-12 for v in r['max_field_absolute_difference'].values())
    checks={k:bool(v) for k,v in checks.items()}
    numerical=all(v for k,v in checks.items() if k!='M0/absolute_endpoint_screen')
    result=dict(stage=stage,passed=all(checks.values()),numerical_checks_passed=numerical,domain_gate=('passed_for_this_screen' if checks['M0/absolute_endpoint_screen'] else 'review_required_background_endpoint_exceeds_frozen_absolute_floor'),checks=checks,metrics=metrics,onsets=onsets,paired_summary=summary,errors=dict(U_dt='unknown',U_z='unknown',U_x='unknown',U_domain='unknown'),t_eff=None)
    write(stage+'_checks.json',result);csvout(stage+'_paired_transport.csv',pairs);csvout(stage+'_pair_summary.csv',summary);csvout(stage+'_performance.csv',perf)
    print(json.dumps(dict(passed=result['passed'],failed=[k for k,v in checks.items() if not v],metrics=metrics,paired_summary=summary),indent=2));assert numerical
    if stage=='bridge':assert result['passed']
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('stage',choices=['fields','bridge','extend']);a=p.parse_args();fields_bridge() if a.stage=='fields' else main(a.stage)
