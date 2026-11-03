"""Compare completed R24/R48 pairs only at accepted physical times."""
import csv,json,math
from pathlib import Path
import numpy as np
from analyze_matched_pair import compare_region
from wp6_paired_stopping import same_time,validate_leg,paired_times
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'analysis/wp7_shortwave_validation_20261103';R24=ROOT/'analysis/wp6_4T_exploratory_20261103';T=136423.0637723484
RUNS={'R24':{'M0':R24/'runs/WP6-M0-explore4T','M1':R24/'runs/WP6-M1-explore4T'},'R48':{'M0':OUT/'runs/WP7-R48-M0-4T','M1':OUT/'runs/WP7-R48-M1-4T'}}
SPATIAL={'R24':{'M0':R24/'spatial_M0','M1':R24/'spatial'},'R48':{'M0':OUT/'spatial_M0','M1':OUT/'spatial_M1'}}
def rows(path):
    with path.open(newline='') as f:return list(csv.DictReader(f))
def numeric_rows(path):return [{k:v if k=='region' else float(v) for k,v in r.items()} for r in rows(path)]
def load(path):return json.loads(path.read_text())
def write_json(name,value):(OUT/name).write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')
def write_csv(name,data):
    fieldnames=[]
    for row in data:
        for key in row:
            if key not in fieldnames:fieldnames.append(key)
    with (OUT/name).open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=fieldnames);w.writeheader();w.writerows(data)
def pairs(a,b,time_key='time_s'):
    i=j=0;result=[]
    while i<len(a) and j<len(b):
        x,y=float(a[i][time_key]),float(b[j][time_key])
        if same_time(x,y):result.append((a[i],b[j],abs(x-y)));i+=1;j+=1
        elif x<y:i+=1
        else:j+=1
    return result
def rel(a,b,floor=1e-30):return abs(a-b)/max(abs(a),abs(b),floor)
def pair_transport(grid):
    budget={role:numeric_rows(RUNS[grid][role]/'transport_budget.csv') for role in ['M0','M1']};meta=load(RUNS[grid]['M1']/'transport_metadata.json');out=[]
    for region in [x['name'] for x in meta['cuts'][:7]]:
        result=compare_region(*[[r for r in budget[role] if r['region']==region] for role in ['M0','M1']],common_accepted_only=True)
        for r in result:r['region']=region
        out.extend(result)
    write_csv(f'paired_transport_{grid}.csv',out);return out
def final_mode_energy(grid,role):
    if grid=='R24':path=ROOT/'tmp/wp6_4T_exploratory_20261103'/('spatial_M0' if role=='M0' else 'spatial')/'selected_fields.npz'
    else:path=ROOT/'tmp/wp7_shortwave_validation_20261103'/f'spatial_{role}_selected_fields.npz'
    d=np.load(path);keys=[k for k in d.files if k.startswith('mode_energy_')];key=max(keys,key=lambda x:int(x.rsplit('_',1)[1]));return np.asarray(d[key],float)
def main():
    acceptance={};weighted={};spatial={};stops={};runtime=[]
    for grid in RUNS:
        acceptance[grid]={};weighted[grid]={};spatial[grid]={};stops[grid]={}
        for role in ['M0','M1']:
            run=RUNS[grid][role];finished=load(run/'finished.json');stop=load(run/'stop_record.json');check=load(run/'post/independent_checks.json')
            budget=rows(run/'transport_budget.csv');max_s=max(abs(float(r['salt_residual_m'])) for r in budget);max_w=max(abs(float(r['water_residual_m'])) for r in budget)
            acceptance[grid][role]=dict(returncode=finished['solver']['returncode'],accepted_steps=stop['accepted_steps'],stop=stop,independent_checks=check['passed'],max_salt_residual_m=max_s,max_water_residual_m=max_w,passed=finished['solver']['returncode']==0 and check['passed'] and max_s<1e-14 and max_w<1e-11)
            weighted[grid][role]=rows(run/'post/weighted_diagnostics.csv');spatial[grid][role]=rows(SPATIAL[grid][role]/'spatial_metrics.csv');stops[grid][role]=stop
            runtime.append(dict(grid=grid,role=role,solver_wall_s=finished['solver']['wall_s'],cpu_s=finished['solver']['cpu_s'],post_archive_wall_s=finished.get('post_archive_wall_s',finished.get('post_wall_s')),accepted_steps=stop['accepted_steps'],archive_bytes=finished['archive']['bytes']))
    write_csv('runtime_cost_all.csv',runtime)
    grid_pairs={role:pairs(weighted['R24'][role],weighted['R48'][role]) for role in ['M0','M1']}
    focus={role:[p for p in grid_pairs[role] if float(p[0]['time_s'])>=3*T-1e-7] for role in ['M0','M1']}
    qnames=['weighted_top_sigma','weighted_fourier','weighted_projection'];role_qnames={'M0':['weighted_top_sigma','weighted_fourier'],'M1':qnames}
    growth={}
    for role in ['M0','M1']:
        growth[role]=dict(shared_states=len(grid_pairs[role]),focus_states=len(focus[role]),max_time_representation_error_s=max(p[2] for p in grid_pairs[role]))
        for q in role_qnames[role]:
            values=[rel(float(a[q]),float(b[q])) for a,b,_ in focus[role]]
            growth[role][q]=dict(max_relative_difference=max(values),median_relative_difference=float(np.median(values)),final_relative_difference=values[-1],gate_20pct=max(values)<=.2)
    checkpoint=[]
    for role in ['M0','M1']:
        for factor in [3,4]:
            cand=min(grid_pairs[role],key=lambda p:abs(float(p[0]['time_s'])-factor*T));a,b,e=cand
            row=dict(role=role,nominal_T=factor,actual_R24_s=float(a['time_s']),actual_R48_s=float(b['time_s']),representation_error_s=e)
            for q in role_qnames[role]+['weighted_top_mean','top_max','domain_max','max_mirror_error','max_endpoint_difference','odd_neumann_energy_fraction']:
                row['R24_'+q]=float(a[q]);row['R48_'+q]=float(b[q]);row[q+'_relative_difference']=rel(float(a[q]),float(b[q]))
            checkpoint.append(row)
    write_csv('R24_R48_checkpoints.csv',checkpoint)
    spectral={};mode_rows=[]
    for role in ['M0','M1']:
        a=spatial['R24'][role][-1];b=spatial['R48'][role][-1];ea=final_mode_energy('R24',role);eb=final_mode_energy('R48',role)
        topa=[int(i+1) for i in np.argsort(ea)[::-1][:3]];topb=[int(i+1) for i in np.argsort(eb)[::-1][:3]]
        ma=int(a['dominant_neumann_mode']);mb=int(b['dominant_neumann_mode']);wa=2*(.6/28)/ma;wb=2*(.6/28)/mb
        spectral[role]=dict(R24_dominant_mode=ma,R48_dominant_mode=mb,R24_wavelength_m=wa,R48_wavelength_m=wb,wavelength_relative_difference=rel(wa,wb),R24_cells_per_wavelength=float(a['dominant_cells_per_wavelength']),R48_cells_per_wavelength=float(b['dominant_cells_per_wavelength']),R24_top_three_modes=topa,R48_top_three_modes=topb,top_three_order_stable=topa==topb,
          R24_odd_fraction=float(a['odd_neumann_energy_fraction']),R48_odd_fraction=float(b['odd_neumann_energy_fraction']),R24_depth90_mm=float(a['depth90_mm']),R48_depth90_mm=float(b['depth90_mm']),R24_depth99_mm=float(a['depth99_mm']),R48_depth99_mm=float(b['depth99_mm']),R24_max_relative_anomaly=float(a['max_relative_anomaly']),R48_max_relative_anomaly=float(b['max_relative_anomaly']))
        for m in range(1,max(len(ea),len(eb))+1):mode_rows.append(dict(role=role,mode=m,wavelength_m=2*(.6/28)/m,R24_energy_fraction=ea[m-1]/ea.sum() if m<=len(ea) else None,R48_energy_fraction=eb[m-1]/eb.sum() if m<=len(eb) else None))
    write_csv('final_mode_energy_R24_R48.csv',mode_rows)
    transport={grid:pair_transport(grid) for grid in RUNS};transport_compare=[];transport_summary={}
    for region in sorted({r['region'] for r in transport['R24']}):
        a=[r for r in transport['R24'] if r['region']==region];b=[r for r in transport['R48'] if r['region']==region];pp=pairs(a,b)
        focusp=[p for p in pp if float(p[0]['time_s'])>=3*T-1e-7];last=pp[-1];x,y,e=last
        row=dict(region=region,actual_R24_s=x['time_s'],actual_R48_s=y['time_s'],representation_error_s=e,R24_delta_C_m=x['delta_C_m'],R48_delta_C_m=y['delta_C_m'],delta_C_relative_difference=rel(x['delta_C_m'],y['delta_C_m'],1e-13),R24_M1_total_down_m_s=x['M1_total_down_m_s'],R48_M1_total_down_m_s=y['M1_total_down_m_s'],R24_M0_total_down_m_s=x['M0_total_down_m_s'],R48_M0_total_down_m_s=y['M0_total_down_m_s'],R24_budget_envelope_m=x['budget_envelope_m'],R48_budget_envelope_m=y['budget_envelope_m'])
        transport_compare.append(row)
        rvals=[rel(p[0]['delta_C_m'],p[1]['delta_C_m'],1e-13) for p in focusp];directions=[math.copysign(1,p[0]['delta_C_m'])==math.copysign(1,p[1]['delta_C_m']) for p in focusp if abs(p[0]['delta_C_m'])>1e-13 or abs(p[1]['delta_C_m'])>1e-13]
        transport_summary[region]=dict(shared_states=len(pp),focus_states=len(focusp),max_delta_C_relative_difference=max(rvals),final_delta_C_relative_difference=rvals[-1],direction_stable=all(directions),R24_resolved=abs(x['delta_C_m'])>max(1e-13,x['budget_envelope_m']),R48_resolved=abs(y['delta_C_m'])>max(1e-13,y['budget_envelope_m']))
    write_csv('transport_R24_R48_final.csv',transport_compare)
    events={}
    for grid in RUNS:
        legs={role:validate_leg(rows(RUNS[grid][role]/'surface_events.csv'),stops[grid][role],acceptance[grid][role]['returncode']) for role in ['M0','M1']};paired=paired_times(legs['M0'],legs['M1'])
        events[grid]=dict(M0={k:legs['M0'][k] for k in ['stop_reason','surface','any_local','interior_first']},M1={k:legs['M1'][k] for k in ['stop_reason','surface','any_local','interior_first']},pair={k:paired[k] for k in ['comparison_status','surface_time_difference','any_local_time_difference','interior_time_difference','different_first_event_types','max_time_representation_error_s']},shared_states=len(paired['indices']))
    m1=spectral['M1'];growth_gate=all(growth['M1'][q]['gate_20pct'] for q in qnames);spectral_gate=m1['wavelength_relative_difference']<=.1 and m1['R48_cells_per_wavelength']>=8 and m1['top_three_order_stable'];complete_4T=all(stops[g][r]['reason']=='requested_horizon' and abs(float(stops[g][r]['last_accepted_time_s'])-4*T)<1e-6 for g in RUNS for r in ['M0','M1'])
    near_limit=m1['R48_cells_per_wavelength']<12
    quantitative_transport=all(x['direction_stable'] and x['max_delta_C_relative_difference']<=.2 and x['R24_resolved'] and x['R48_resolved'] for x in transport_summary.values())
    decision=dict(scope='command/11 completed R48 pair',complete_4T=complete_4T,numerical_acceptance=all(acceptance[g][r]['passed'] for g in acceptance for r in acceptance[g]),spectral_convergence=spectral_gate,growth_convergence_20pct=growth_gate,near_resolution_limit=near_limit,quantitative_transport_convergence=quantitative_transport,
      R96='PROPOSE_ONLY_NO_GO' if (not spectral_gate or near_limit) else 'NOT_REQUIRED_BY_CURRENT_GRID_GATE',multiwave='NO-GO pending converged background/symmetry and local mechanism',later_8T='NO-GO pending R48/R96 decision and threshold interpretation',t_eff=None,crystallization_delay_label=None,training_labels_created=False)
    result=dict(acceptance=acceptance,growth=growth,spectral=spectral,transport=transport_summary,events=events,stops=stops,decision=decision)
    write_json('r48_comparison.json',result);write_json('decision.json',decision);print(json.dumps(result,indent=2))
if __name__=='__main__':main()
