"""Frozen background gates. Always records failure; never launches a solver."""
import json,re,sys
import numpy as np
from wp8_strict_campaign import ROOT,OUT,CACHE,sha,save,T
from analyze_wp8 import rows,pairs
from single_wavelength_post import widths,params
from wp6_paired_stopping import validate_leg,same_time

def accepted_prefix_audit(ev,stop,logs):
    """Audit saved accepted states without declaring the interrupted run complete."""
    tt=[float(r['time_s']) for r in ev]
    valid=bool(ev) and tt[0]==0 and stop['accepted_steps']==len(ev)-1 and same_time(float(stop['last_accepted_time_s']),tt[-1])
    for i,r in enumerate(ev):
        valid=valid and int(r['step'])==i and all(np.isfinite(float(r[k])) for k in ['time_s','dt_s','x_top_max','x_domain_max'])
        if i:valid=valid and tt[i]>tt[i-1] and same_time(tt[i],tt[i-1]+float(r['dt_s']))
        valid=valid and float(r['x_domain_max'])>=float(r['x_top_max'])
        for col,flag in [('x_top_max','surface_threshold_reached'),('x_domain_max','global_threshold_reached')]:
            crossed=float(r[col])>=.1;valid=valid and bool(int(r[flag]))==crossed and (not crossed or i==len(ev)-1)
    successes=list(re.finditer(r'WP8_CONSTRAINT iterations=',logs));steps=re.findall(r'Time step (\d+) done',logs)
    valid=valid and len(successes)==len(ev)-1 and steps==[str(i) for i in range(1,len(ev))]
    if 'Paper 2 run failed' in logs:
        tail=logs[successes[-1].end():] if successes else logs
        valid=valid and stop['reason']=='exception_or_incomplete' and 'WP8_ITER 20 ' in tail and 'WP8_CONSTRAINT' not in tail
    return bool(valid)

def validate(stage):
    d=OUT/stage;plan=json.loads((d/'plan.json').read_text());g=plan['gates'];finished=json.loads((d/'finished.json').read_text());stop=json.loads((d/'stop_record.json').read_text());ev=rows(d/'surface_events.csv');bud=rows(d/'transport_budget.csv');ff=json.loads((d/'post/field_manifest.json').read_text());data=dict(np.load(CACHE/stage/'post/profiles.npz'));par=params(d/'params.input');phi=float(par['Problem.Porosity']);checks={};metrics={};issues=[]
    checks['successful_exit']=finished['returncode']==0 and finished['trigger'] is None
    try:events=validate_leg(ev,stop,finished['returncode']);checks['accepted_times_and_stopping']=True
    except ValueError as e:events=dict(status='failed',reason=str(e));checks['accepted_times_and_stopping']=False;issues.append(str(e))
    checks['all_accepted_fields']=len(ev)==len(ff)==len(data['x_mean'])
    metrics['max_horizontal_x_range']=max(v['x_range'] for v in ff);metrics['max_horizontal_p_range']=max(v['p_range'] for v in ff)
    checks['frozen_uniformity']=metrics['max_horizontal_x_range']<=g['lateral_range_max'];checks['exact_equality']=metrics['max_horizontal_x_range']==0 and metrics['max_horizontal_p_range']==0
    metrics['salt_residual_m']=max(abs(float(r['salt_residual_m'])) for r in bud);metrics['water_residual_m']=max(abs(float(r['water_residual_m'])) for r in bud)
    checks['salt_closure']=metrics['salt_residual_m']<=g['salt_residual_m'];checks['water_closure']=metrics['water_residual_m']<=g['water_residual_m']
    checks['input_hash']=sha(d/'params.input')==plan['input_sha256']==finished['input_sha256'];checks['elf_hash']=finished['elf_sha256']==plan['elf_sha256'];checks['worker_hash']=finished['worker_sha256']==sha(ROOT/'scripts/wp8_strict_worker.py')
    receipt=json.loads((d/'archive_receipt.json').read_text());checks['archive_verified']=receipt['verified_frames']==len(ff) and receipt['archive_sha256']==finished['archive']['sha256'] and receipt['profile_sha256']==sha(CACHE/stage/'post/profiles.npz')
    checks['wall_budget']=finished['solver_and_live_audit_wall_s']<=plan['wall_cap_s'];checks['archive_budget']=finished['archive_budget_pass'];checks['memory_budget']=finished['peak_child_rss_kib']<=plan['memory_gib']*1024**2
    logs=(d/'solver.log').read_text();accepted=re.findall(r'WP8_CONSTRAINT iterations=(\d+) residual=([^ ]+) shift=([^\s]+)',logs);checks['one_Newton_success_per_step']=len(accepted)==len(ev)-1;checks['no_failed_Newton_attempt']='Paper 2 run failed' not in logs and not re.search(r'WP8_ITER 20 ',logs)
    checks['accepted_prefix_time_event_Newton_integrity']=accepted_prefix_audit(ev,stop,logs)
    metrics['max_Newton_iterations']=max((int(v[0]) for v in accepted),default=None)
    dt=float(par['TimeLoop.MaxTimeStepSize']);checks['fixed_accepted_dt']=all(abs(float(r['dt_s'])-dt)<1e-7 for r in ev[1:])
    checks['field_event_maxima']=all(abs(f['max_x']-float(r['x_domain_max']))<=1e-14 for f,r in zip(ff,ev))
    zs=data['z_m'];wz=widths(zs);storage=[];waterstorage=[];reintegration=[]
    for region in sorted({r['region'] for r in bud}):
        rr=[r for r in bud if r['region']==region];cut=float(rr[0]['actual_z_over_H'])*.2;mask=zs>=cut;cc=0.
        for i,r in enumerate(rr):
            storage.append(abs(phi*(data['sx_mean'][i,mask]@wz[mask])-float(r['inventory_salt_m'])));waterstorage.append(abs(phi*(data['sw_mean'][i,mask]@wz[mask])-float(r['inventory_water_m'])))
            if i:cc+=float(r['dt_s'])*float(r['cut_salt_total_down_m_s'])
            reintegration.append(abs(cc-float(r['cumulative_cut_salt_m'])))
    metrics.update(storage_reconstruction_salt_m=max(storage),storage_reconstruction_water_m=max(waterstorage),flux_reintegration_m=max(reintegration))
    checks['storage_reconstruction']=max(storage)<=g['storage_reconstruction_salt_m'];checks['flux_reintegration']=max(reintegration)<=g['flux_reintegration_abs_m']
    x0=float(par['Problem.SaltMoleFraction']);mass=x0*float(par['Problem.MolarMassNaCl'])/((1-x0)*float(par['Problem.MolarMassH2O'])+x0*float(par['Problem.MolarMassNaCl']));rho=float(par['Problem.PureWaterDensity'])+1000*mass*(float(par['Problem.DensityLinearCoefficient'])+mass*float(par['Problem.DensityQuadraticCoefficient']));p0=1e5-(float(par['Problem.WaterLevel'])+zs)*rho*float(par['Problem.GravitationalAcceleration'])
    metrics['initial_x_error']=float(abs(data['x_mean'][0]-x0).max());metrics['initial_p_error_Pa']=float(abs(data['pw_mean'][0]-p0).max());checks['initial_x']=metrics['initial_x_error']<=g['initial_x_abs'];checks['initial_p']=metrics['initial_p_error_Pa']<=g['initial_pressure_Pa']
    if stage=='T':
        ref=ROOT/'analysis/wp7_shortwave_validation_20261103/runs/WP7-1D-background-T';refdata=dict(np.load(ROOT/'tmp/wp8_zero_mode_20261103/B/profiles_spectra.npz'));reft=refdata['time_s'];refx=refdata['x_mean'];label='legacy_background'
    else:
        ref=OUT/'T';refdata=dict(np.load(CACHE/'T/post/profiles.npz'));reft=[float(r['time_s']) for r in rows(ref/'surface_events.csv')];refx=refdata['x_mean'];label='accepted_T_run'
    tt=[float(r['time_s']) for r in ev];pp=pairs(reft,tt);checks['common_states_present']=len(pp)>0;checks['same_vertical_nodes']=np.array_equal(refdata['z_m'],zs);metrics['shared_reference_states']=len(pp);metrics['reference_comparison_last_s']=float(tt[pp[-1][1]])
    metrics['profile_bridge_error']=max(float(abs(refx[i]-data['x_mean'][j]).max()) for i,j in pp);checks['profile_bridge']=metrics['profile_bridge_error']<=g['profile_T_abs']
    old=rows(ref/'transport_budget.csv');bridge={}
    for key in ['inventory_salt_m','inventory_water_m','cut_salt_adv_down_m_s','cut_salt_diff_down_m_s','cut_salt_total_down_m_s','cumulative_cut_salt_m']:
        errors=[]
        for region in {r['region'] for r in bud}:
            aa=[r for r in old if r['region']==region];bb=[r for r in bud if r['region']==region];pairs_=pairs([float(r['time_s']) for r in aa],[float(r['time_s']) for r in bb])
            for i,j in pairs_:
                assert aa[i]['actual_z_over_H']==bb[j]['actual_z_over_H'];errors.append(abs(float(aa[i][key])-float(bb[j][key])))
        bridge[key]=max(errors)
    checks['inventory_bridge']=bridge['inventory_salt_m']<=g['inventory_match_abs_m'];checks['total_flux_bridge']=bridge['cut_salt_total_down_m_s']<=g['flux_total_abs_m_s']
    oldinitial={r['region']:float(r['inventory_salt_m']) for r in old if r['step']=='0'};metrics['initial_inventory_relative_error']=max(abs(float(r['inventory_salt_m'])-oldinitial[r['region']])/oldinitial[r['region']] for r in bud if r['step']=='0');checks['initial_inventory']=metrics['initial_inventory_relative_error']<=g['initial_storage_relative']
    if stage=='4T':
        refmanifest=json.loads((ref/'post/field_manifest.json').read_text());metrics['first_T_equal_field_hashes']=sum(refmanifest[i]['sha256']==ff[j]['sha256'] for i,j in pp);checks['first_T_fields_byte_identical']=metrics['first_T_equal_field_hashes']==len(pp);checks['first_T_budgets_identical']=all(v==0 for v in bridge.values())
    oldpar=params(ROOT/'analysis/wp7_shortwave_validation_20261103/runs/WP7-1D-background-T/params.input');changes=[k for k in par if par[k]!=oldpar.get(k)];checks['frozen_physics']=set(changes)<= {'Problem.Name','TimeLoop.TEnd'}
    checks={k:bool(v) for k,v in checks.items()};prefix_exclusions=['successful_exit','accepted_times_and_stopping','no_failed_Newton_attempt'];prefix_ok=all(v for k,v in checks.items() if k not in prefix_exclusions)
    r=dict(passed=all(checks.values()),accepted_prefix_passed=prefix_ok,accepted_prefix_last_time_s=float(ev[-1]['time_s']),prefix_completion_checks_excluded=prefix_exclusions,stage=stage,checks=checks,metrics=metrics,bridge_reference=label,bridge_absolute_errors=bridge,events=events,issues=issues,cost=finished,changed_input_keys=changes,failures=[k for k,v in checks.items() if not v],no_automatic_retry=True)
    target=d/'validation.json'
    if stage=='T' and (OUT/'4T/plan.json').exists():
        approved=json.loads((OUT/'4T/plan.json').read_text())['previous_T_gate_sha256'];assert sha(target)==approved, 'approved T gate was modified'
        target=d/'validation_recheck.json'  # Do not rewrite the gate that admitted 4T.
    save(target,r);print(json.dumps({k:v for k,v in r.items() if k not in ['cost','events']},indent=2));return r
if __name__=='__main__':sys.exit(0 if validate(sys.argv[1])['passed'] else 1)
