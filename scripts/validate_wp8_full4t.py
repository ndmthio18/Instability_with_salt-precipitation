"""Command 15 acceptance gates for the single continuous corrected strict background.

Audits every accepted step independently of the solver and never launches anything.
"""
import json,re,sys
import numpy as np
import wp8_full4t_campaign as f
from analyze_wp8 import rows,pairs
from single_wavelength_post import widths,params
from wp6_paired_stopping import validate_leg,same_time

ROOT=f.ROOT;OUT=f.OUT;CACHE=f.CACHE;STAGE=f.STAGE;sha=f.sha;save=f.save
HIST=ROOT/'analysis/wp8_strict_background_20261103'
WP7BASE=ROOT/'analysis/wp7_shortwave_validation_20261103/runs/WP7-1D-background-T/params.input'


def load(stage_dir,cache_dir):
    d=dict(ev=rows(stage_dir/'surface_events.csv'),bud=rows(stage_dir/'transport_budget.csv'),
           ff=json.loads((stage_dir/'post/field_manifest.json').read_text()),
           par=params(stage_dir/'params.input'),prof=dict(np.load(cache_dir/'post/profiles.npz')))
    return d


def analytic_initial(par,zs):
    x0=float(par['Problem.SaltMoleFraction'])
    mass=x0*float(par['Problem.MolarMassNaCl'])/((1-x0)*float(par['Problem.MolarMassH2O'])+x0*float(par['Problem.MolarMassNaCl']))
    rho=float(par['Problem.PureWaterDensity'])+1000*mass*(float(par['Problem.DensityLinearCoefficient'])+mass*float(par['Problem.DensityQuadraticCoefficient']))
    p0=1e5-(float(par['Problem.WaterLevel'])+zs)*rho*float(par['Problem.GravitationalAcceleration'])
    return x0,p0


def validate():
    d=OUT/STAGE;plan=json.loads((d/'plan.json').read_text());g=plan['gates'];crit=plan['criterion']
    fin=json.loads((d/'finished.json').read_text());stop=json.loads((d/'stop_record.json').read_text())
    cur=load(d,CACHE/STAGE);ev=cur['ev'];bud=cur['bud'];ff=cur['ff'];par=cur['par'];a=cur['prof']
    di=rows(d/'convergence_v2.csv');accepted=[r for r in di if int(r['accepted'])]
    checks={};met={};issues=[]

    # --- completion and stopping semantics -------------------------------
    checks['successful_exit']=fin['returncode']==0 and fin['trigger'] is None
    try:
        events=validate_leg(ev,stop,fin['returncode']);checks['accepted_times_and_stopping']=True
    except ValueError as e:
        events=dict(status='failed',reason=str(e));checks['accepted_times_and_stopping']=False;issues.append(str(e))
    checks['reached_requested_horizon']=stop['reason']=='requested_horizon' and stop['accepted_steps']==plan['steps_cap']
    checks['bounded_steps']=stop['accepted_steps']<=plan['steps_cap']
    checks['frames_steps']=len(ev)==len(ff)==len(a['x_mean'])==stop['accepted_steps']+1
    met['accepted_steps']=stop['accepted_steps'];met['last_accepted_time_s']=float(ev[-1]['time_s'])
    met['horizon_shortfall_s']=plan['horizon_s']-float(ev[-1]['time_s'])
    # The solver accumulates t += dt in double precision, so the endpoint cannot equal
    # the decimal TEnd to a few ULP after 4800 additions. Test the accumulation itself
    # (exact, reproducible) and bound the residual drift by that accumulation budget.
    acc=0.;drift=0.
    for i,r in enumerate(ev):
        if i:acc+=plan['dt_s']
        drift=max(drift,abs(float(r['time_s'])-acc))
    met['fixed_dt_accumulation_drift_s']=drift
    met['horizon_accumulation_budget_s']=plan['steps_cap']*np.spacing(plan['horizon_s'])
    checks['exact_fixed_dt_accumulation']=drift==0.
    checks['horizon_within_accumulation']=abs(met['horizon_shortfall_s'])<=met['horizon_accumulation_budget_s'] and float(ev[-1]['time_s'])<=plan['horizon_s']+1e-7

    # --- horizontal uniformity at every accepted step --------------------
    met['max_horizontal_x_range']=max(v['x_range'] for v in ff)
    met['max_horizontal_p_range']=max(v['p_range'] for v in ff)
    checks['frozen_uniformity']=met['max_horizontal_x_range']<=g['lateral_range_max']
    checks['exact_horizontal_equality']=met['max_horizontal_x_range']==g['expected_lateral_range'] and met['max_horizontal_p_range']==g['expected_lateral_range']

    # --- fixed timestep and threshold status ------------------------------
    dt=float(par['TimeLoop.MaxTimeStepSize'])
    checks['fixed_accepted_dt']=all(abs(float(r['dt_s'])-dt)<1e-7 for r in ev[1:]) and abs(dt-plan['dt_s'])<=0
    checks['threshold_flags']=all((float(r['x_top_max'])>=plan['threshold'])==bool(int(r['surface_threshold_reached']))
        and (float(r['x_domain_max'])>=plan['threshold'])==bool(int(r['global_threshold_reached']))
        and float(r['x_domain_max'])>=float(r['x_top_max']) for r in ev)
    met['max_surface_x']=max(float(r['x_top_max']) for r in ev)
    met['max_domain_x']=max(float(r['x_domain_max']) for r in ev)
    met['threshold_observed']=bool(met['max_domain_x']>=plan['threshold'])
    checks['no_accepted_state_after_crossing']=(not met['threshold_observed']) or bool(int(ev[-1]['global_threshold_reached']) or int(ev[-1]['surface_threshold_reached']))
    checks['field_event_maxima']=all(abs(x['max_x']-float(r['x_domain_max']))<=1e-14 for x,r in zip(ff,ev))

    # --- Newton acceptance -------------------------------------------------
    log=(d/'solver.log').read_text()
    succ=re.findall(r'WP8_CONSTRAINT iterations=(\d+) residual=([^ ]+) shift=([^\s]+)',log)
    checks['one_Newton_success_per_step']=len(succ)==stop['accepted_steps']
    checks['no_failed_Newton_attempt']='Paper 2 run failed' not in log and not re.search(r'WP8_ITER 20 ',log)
    checks['accepted_records_match_steps']=len(accepted)==stop['accepted_steps'] and [int(r['step']) for r in accepted]==list(range(1,len(ev)))
    checks['no_accepted_record_beyond_last_iteration']=all(int(r['iteration'])<plan['criterion']['max_Newton_iteration'] for r in accepted)
    it=np.array([int(r['iteration'])+1 for r in accepted])
    met['Newton_iterations_min']=int(it.min());met['Newton_iterations_max']=int(it.max());met['Newton_iterations_mean']=float(it.mean())
    met['Newton_iterations_histogram']={str(int(k)):int(v) for k,v in zip(*np.unique(it,return_counts=True))}
    checks['accepted_scaled_contract']=all(
        float(r['relative_shift'])<=crit['relative_shift'] and
        float(r['scaled_nonlinear'])<=crit['scaled_residual'] and
        float(r['linear_eta_component'])<=crit['componentwise_linear_backward_error'] and
        float(r['regional_salt_error_m'])<=crit['regional_step_salt_m'] and
        float(r['regional_water_error_m'])<=crit['regional_step_water_m'] for r in accepted)
    for k in ['water_max','salt_max','scaled_nonlinear','linear_eta_component','relative_shift','regional_salt_error_m','regional_water_error_m','min_pivot','pivot_growth']:
        met['accepted_max_'+k]=max(float(r[k]) for r in accepted)
    met['accepted_min_min_pivot']=min(float(r['min_pivot']) for r in accepted)

    # --- conservative diagnostics, reconstructed independently -------------
    met['max_salt_closure_m']=max(abs(float(r['salt_residual_m'])) for r in bud)
    met['max_water_closure_m']=max(abs(float(r['water_residual_m'])) for r in bud)
    checks['salt_closure']=met['max_salt_closure_m']<=g['salt_residual_m']
    checks['water_closure']=met['max_water_closure_m']<=g['water_residual_m']
    phi=float(par['Problem.Porosity']);z=a['z_m'];w=widths(z)
    reg=sorted({r['region'] for r in bud})
    checks['seven_planes_plus_free_domain']=len(reg)==8 and 'free_domain' in reg and sum(x.startswith('upper_') for x in reg)==7
    storage=[];water=[];rein=[];ident=[];wrein=[];wident=[]
    for region in reg:
        rr=[r for r in bud if r['region']==region]
        checks['regional_rows_'+region]=len(rr)==len(ev) and all(same_time(float(r['time_s']),float(e['time_s'])) for r,e in zip(rr,ev))
        mask=z>=float(rr[0]['actual_z_over_H'])*.2;C=0.;CW=0.;I0=float(rr[0]['inventory_salt_m']);W0=float(rr[0]['inventory_water_m'])
        for i,r in enumerate(rr):
            storage.append(abs(phi*a['sx_mean'][i,mask]@w[mask]-float(r['inventory_salt_m'])))
            water.append(abs(phi*a['sw_mean'][i,mask]@w[mask]-float(r['inventory_water_m'])))
            if i:
                C+=float(r['dt_s'])*float(r['cut_salt_total_down_m_s'])
                CW+=float(r['dt_s'])*(float(r['cut_water_down_m_s'])+float(r['other_water_out_m_s']))
            rein.append(abs(C-float(r['cumulative_cut_salt_m'])))
            ident.append(abs((float(r['inventory_salt_m'])-I0+C)-float(r['salt_residual_m'])))
            wrein.append(abs(CW-float(r['cumulative_net_water_out_m'])))
            wident.append(abs(float(r['inventory_water_m'])-W0+CW-float(r['water_residual_m'])))
    met.update(storage_salt_error_m=max(storage),storage_water_error_m=max(water),
               flux_reintegration_error_m=max(rein),budget_identity_error_m=max(ident),
               water_flux_reintegration_error_m=max(wrein),water_budget_identity_error_m=max(wident))
    checks['independent_storage']=max(storage)<=g['storage_reconstruction_salt_m']
    checks['independent_flux']=max(rein)<=g['flux_reintegration_abs_m']
    checks['independent_budget']=max(ident)<=1e-14
    checks['independent_water_storage']=max(water)<=g['water_residual_m']
    checks['independent_water_flux']=max(wrein)<=g['water_residual_m']
    checks['independent_water_budget']=max(wident)<=g['water_residual_m']

    # --- original initial state and frozen physics -------------------------
    x0,p0=analytic_initial(par,z)
    met['initial_x_error']=float(abs(a['x_mean'][0]-x0).max());met['initial_p_error_Pa']=float(abs(a['pw_mean'][0]-p0).max())
    checks['initial_x']=met['initial_x_error']<=g['initial_x_abs'];checks['initial_p']=met['initial_p_error_Pa']<=g['initial_pressure_Pa']
    base=params(WP7BASE);diff=[k for k in par if par[k]!=base.get(k)]
    checks['frozen_physics']=all(k in ['Problem.Name','TimeLoop.TEnd'] or k.startswith('WP8.') for k in diff) and all(k in par for k in base)
    checks['identical_TEnd_value']=float(par['TimeLoop.TEnd'])==plan['horizon_s']
    checks['no_replay_or_restart']=par.get('WP8.AllowDiagnosticReplay')=='false' and not plan['restart'] and not plan['replay'] and not (d/'replay.dat').exists()

    # --- provenance and budgets --------------------------------------------
    rec=json.loads((d/'archive_receipt.json').read_text())
    checks['hashes']=sha(d/'params.input')==fin['input_sha256']==plan['input_sha256'] and fin['elf_sha256']==plan['elf_sha256'] and fin['worker_sha256']==sha(ROOT/'scripts/wp8_strict_worker.py')
    checks['authorization_hash']=sha(f.AUTH)==plan['authorization_sha256']
    checks['archive_verified']=rec['verified_frames']==len(ff) and rec['archive_sha256']==fin['archive']['sha256'] and rec['profile_sha256']==sha(CACHE/STAGE/'post/profiles.npz')
    checks['wall_budget']=fin['solver_and_live_audit_wall_s']<=plan['wall_cap_s']
    checks['archive_budget']=fin['archive']['bytes']<=plan['archive_cap_mib']*1024**2
    checks['memory_budget']=fin['peak_child_rss_kib']<=plan['memory_gib']*1024**2
    checks['single_core']=fin['cpu'] is not None and plan['cpus']==1

    # --- did the corrected solver pass the old failure window? -------------
    hstop=json.loads((HIST/'4T/stop_record.json').read_text())
    met['historical_last_accepted_step']=hstop['accepted_steps']
    met['historical_last_accepted_time_s']=float(hstop['last_accepted_time_s'])
    met['historical_failed_attempt_step']=hstop['accepted_steps']+1
    checks['passed_old_failure_window']=stop['accepted_steps']>hstop['accepted_steps'] and fin['returncode']==0
    win=[r for r in accepted if hstop['accepted_steps']-9<=int(r['step'])<=hstop['accepted_steps']+10]
    met['failure_window_steps']=[int(r['step']) for r in win]
    met['failure_window_iterations']=[int(r['iteration'])+1 for r in win]
    met['failure_window_max_scaled_nonlinear']=max(float(r['scaled_nonlinear']) for r in win)
    met['failure_window_max_salt_residual']=max(float(r['salt_max']) for r in win)

    checks={k:bool(v) for k,v in checks.items()}
    met={k:(float(v) if isinstance(v,np.floating) else v) for k,v in met.items()}
    out=dict(passed=all(checks.values()),stage=STAGE,command='command/15.md',
             checks=checks,metrics=met,events=events,issues=issues,
             failures=[k for k,v in checks.items() if not v],
             censored_by_threshold=bool(met['threshold_observed']),
             changed_input_keys=diff,no_automatic_retry=True,cost=fin)
    save(OUT/'validation.json',out)
    print(json.dumps({k:v for k,v in out.items() if k not in ['cost','events']},indent=2,default=str))
    return out


if __name__=='__main__':
    sys.exit(0 if validate()['passed'] else 1)
