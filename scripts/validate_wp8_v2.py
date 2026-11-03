"""Independent accepted-field, budget, event and equation audit; no launches."""
import json,re,csv
import numpy as np
import wp8_v2_regression as v
from analyze_wp8 import rows,pairs
from single_wavelength_post import widths,params
from wp6_paired_stopping import validate_leg,same_time
c=v.c

def audit(name):
    d=c.OUT/name;p=json.loads((d/'plan.json').read_text());f=json.loads((d/'finished.json').read_text());s=json.loads((d/'stop_record.json').read_text());ev=rows(d/'surface_events.csv');bud=rows(d/'transport_budget.csv');di=rows(d/'convergence_v2.csv');ff=json.loads((d/'post/field_manifest.json').read_text());ar=json.loads((d/'archive_receipt.json').read_text());a=dict(np.load(c.CACHE/name/'post/profiles.npz'));par=params(d/'params.input');g=p['gates'];checks={};met={}
    accepted=[r for r in di if int(r['accepted'])]; failed=f['returncode']!=0
    expected=(name=='forced_failure')
    checks['exit_semantics']=(failed if expected else not failed) if name!='late_legacy_window' else True
    checks['bounded_steps']=s['accepted_steps']<=p['steps_cap']
    checks['frames_steps']=len(ev)==len(ff)==len(a['x_mean'])==s['accepted_steps']+1
    checks['one_success_per_accepted_step']=len(accepted)==s['accepted_steps'] and [int(r['step']) for r in accepted]==list(range(1,len(ev)))
    log=(d/'solver.log').read_text();checks['solver_success_records']=log.count('WP8_CONSTRAINT iterations=')==s['accepted_steps']
    checks['failed_path_no_acceptance']=not failed or (s['reason']=='exception_or_incomplete' and int(di[-1]['accepted'])==0 and int(di[-1]['step'])==len(ev))
    checks['forced_failure_no_step']=not expected or s['accepted_steps']==0
    checks['time_axis']=all(int(r['step'])==i and (float(r['time_s'])==0 if not i else same_time(float(r['time_s']),float(ev[i-1]['time_s'])+float(r['dt_s']))) for i,r in enumerate(ev)) and same_time(float(s['last_accepted_time_s']),float(ev[-1]['time_s']))
    checks['threshold_flags']=all((float(r['x_top_max'])>=.1)==bool(int(r['surface_threshold_reached'])) and (float(r['x_domain_max'])>=.1)==bool(int(r['global_threshold_reached'])) and float(r['x_domain_max'])>=float(r['x_top_max']) for r in ev)
    events={'status':'failed_expected' if expected else 'failed','accepted_prefix_only':True}
    if not failed:events=validate_leg(ev,s,f['returncode'])
    checks['maxima_field_event']=all(abs(x['max_x']-float(r['x_domain_max']))<=1e-14 for x,r in zip(ff,ev))
    checks['exact_horizontal_equality']=all(x['x_range']==0 and x['p_range']==0 for x in ff)
    checks['fixed_dt']=all(abs(float(r['dt_s'])-113.68588647695695)<1e-7 for r in ev[1:])
    checks['hashes']=c.sha(d/'params.input')==f['input_sha256']==p['input_sha256'] and f['elf_sha256']==p['elf_sha256'] and f['worker_sha256']==c.sha(c.ROOT/'scripts/wp8_strict_worker.py')
    checks['archive']=ar['verified_frames']==len(ff) and ar['archive_sha256']==f['archive']['sha256'] and ar['profile_sha256']==c.sha(c.CACHE/name/'post/profiles.npz')
    checks['cost_caps']=f['solver_and_live_audit_wall_s']<=p['wall_cap_s'] and f['archive']['bytes']<=p['archive_cap_mib']*1024**2 and f['peak_child_rss_kib']<=2*1024**2
    base=params(c.ROOT/'analysis/wp7_shortwave_validation_20261103/runs/WP7-1D-background-T/params.input');diff=[k for k in par if par[k]!=base.get(k)]
    checks['frozen_physics']=all(k in ['Problem.Name','TimeLoop.TEnd'] or k.startswith('WP8.') for k in diff) and all(k in par for k in base)
    met['max_salt_closure_m']=max(abs(float(r['salt_residual_m'])) for r in bud);met['max_water_closure_m']=max(abs(float(r['water_residual_m'])) for r in bud)
    checks['salt_closure']=met['max_salt_closure_m']<=g['salt_residual_m'];checks['water_closure']=met['max_water_closure_m']<=g['water_residual_m']
    phi=float(par['Problem.Porosity']);z=a['z_m'];w=widths(z);storage=[];water=[];rein=[];residual=[];waterrein=[];waterresidual=[];reg=sorted({r['region'] for r in bud});checks['seven_planes_plus_free_domain']=len(reg)==8 and 'free_domain' in reg and sum(x.startswith('upper_') for x in reg)==7
    for region in reg:
        rr=[r for r in bud if r['region']==region];checks['regional_rows_'+region]=len(rr)==len(ev) and all(same_time(float(r['time_s']),float(e['time_s'])) for r,e in zip(rr,ev));mask=z>=float(rr[0]['actual_z_over_H'])*.2;C=0;CW=0
        for i,r in enumerate(rr):
            storage.append(abs(phi*a['sx_mean'][i,mask]@w[mask]-float(r['inventory_salt_m'])));water.append(abs(phi*a['sw_mean'][i,mask]@w[mask]-float(r['inventory_water_m'])))
            if i:C+=float(r['dt_s'])*float(r['cut_salt_total_down_m_s'])
            if i:CW+=float(r['dt_s'])*(float(r['cut_water_down_m_s'])+float(r['other_water_out_m_s']))
            rein.append(abs(C-float(r['cumulative_cut_salt_m'])))
            residual.append(abs((float(r['inventory_salt_m'])-float(rr[0]['inventory_salt_m'])+C)-float(r['salt_residual_m'])))
            waterrein.append(abs(CW-float(r['cumulative_net_water_out_m'])))
            waterresidual.append(abs(float(r['inventory_water_m'])-float(rr[0]['inventory_water_m'])+CW-float(r['water_residual_m'])))
    met.update(storage_salt_error_m=max(storage),storage_water_error_m=max(water),flux_reintegration_error_m=max(rein),budget_identity_error_m=max(residual))
    checks['independent_storage']=max(storage)<=g['storage_reconstruction_salt_m'];checks['independent_flux']=max(rein)<=g['flux_reintegration_abs_m'];checks['independent_budget']=max(residual)<=1e-14
    met.update(water_flux_reintegration_error_m=max(waterrein),water_budget_identity_error_m=max(waterresidual))
    checks['independent_water_storage']=max(water)<=g['water_residual_m'];checks['independent_water_flux']=max(waterrein)<=g['water_residual_m'];checks['independent_water_budget']=max(waterresidual)<=g['water_residual_m']
    checks['accepted_scaled_contract']=True
    if p['criterion']=='scaled':
        checks['accepted_scaled_contract']=all(float(r['relative_shift'])<=1e-10 and float(r['scaled_nonlinear'])<=p['roundoff_factor']*np.finfo(float).eps and float(r['linear_eta_component'])<=1e-12 and float(r['regional_salt_error_m'])<=1e-12/4800 and float(r['regional_water_error_m'])<=1e-10/4800 for r in accepted)
    if accepted:
        for key in ['water_max','salt_max','scaled_nonlinear','linear_eta_component','relative_shift','regional_salt_error_m','regional_water_error_m']:met['accepted_max_'+key]=max(float(r[key]) for r in accepted)
    # ASCII diagnostic replay is compared as a short window, never called a restart.
    old,index=v.CASES[name][:2];refstage=old or 'T';ref=c.ROOT/f'analysis/wp8_strict_background_20261103/{refstage}';ra=dict(np.load(c.ROOT/f'tmp/wp8_strict_background_20261103/{refstage}/post/profiles.npz'));re=rows(ref/'surface_events.csv');rb=rows(ref/'transport_budget.csv')
    pp=pairs([float(r['time_s']) for r in re],[float(r['time_s'])+p['physical_time_offset_s'] for r in ev]);checks['reference_common_initial']=bool(pp) and pp[0]==(index,0);checks['reference_nodes']=np.array_equal(z,ra['z_m']);met['reference_states']=len(pp)
    met['profile_bridge']=max(float(abs(ra['x_mean'][i]-a['x_mean'][j]).max()) for i,j in pp);met['pressure_bridge_Pa']=max(float(abs(ra['pw_mean'][i]-a['pw_mean'][j]).max()) for i,j in pp)
    checks['profile_bridge']=met['profile_bridge']<=g['profile_T_abs'];checks['initial_pressure']=float(abs(ra['pw_mean'][index]-a['pw_mean'][0]).max())<=g['initial_pressure_Pa'];checks['initial_x']=float(abs(ra['x_mean'][index]-a['x_mean'][0]).max())<=g['initial_x_abs']
    for key,tol in [('inventory_salt_m',g['inventory_match_abs_m']),('cut_salt_total_down_m_s',g['flux_total_abs_m_s'])]:
        errors=[]
        for region in reg:
            aa=[r for r in rb if r['region']==region];bb=[r for r in bud if r['region']==region]
            errors.extend(abs(float(aa[i][key])-float(bb[j][key])) for i,j in pp)
        met['bridge_'+key]=max(errors);checks['bridge_'+key]=max(errors)<=tol
    initial=[]
    for region in reg:
        aa=[r for r in rb if r['region']==region][index];bb=[r for r in bud if r['region']==region][0];initial.append(abs(float(aa['inventory_salt_m'])-float(bb['inventory_salt_m']))/abs(float(aa['inventory_salt_m'])))
    met['initial_inventory_relative']=max(initial);checks['initial_inventory']=max(initial)<=g['initial_storage_relative']
    checks={k:bool(b) for k,b in checks.items()};met={k:float(x) if isinstance(x,np.floating) else x for k,x in met.items()}
    result=dict(passed=all(checks.values()),checks=checks,metrics=met,failed_solver=failed,events=events,physical_time_offset_s=p['physical_time_offset_s'],diagnostic_replay=bool(old),restart_equivalence=False,changed_keys=diff,cost=f)
    c.save(d/'validation.json',result);return result

def equations():
    d=c.OUT/'late_equation_probe';a=np.loadtxt(d/'jacobian.csv',delimiter=',');u,r,du=np.loadtxt(d/'state.csv',delimiter=',').T
    scale=np.maximum(abs(u),np.tile([1.,.0036],len(u)//2));A=a*scale[None,:];row=np.sum(abs(A),axis=1);scaled=A/row[:,None]
    e=r.astype(np.longdouble)-a.astype(np.longdouble)@du.astype(np.longdouble);den=abs(r.astype(np.longdouble))+abs(a.astype(np.longdouble))@abs(du.astype(np.longdouble))
    fd=rows(d/'fd_columns.csv')+rows(c.OUT/'late_equation_probe_fine/fd_columns.csv');errors=[]
    for col in sorted({int(x['column']) for x in fd}):
        for h in sorted({float(x['h']) for x in fd if int(x['column'])==col}):
            rr=[x for x in fd if int(x['column'])==col and float(x['h'])==h];f=np.array([float(x['fd']) for x in rr]);j=np.array([float(x['assembled']) for x in rr]);errors.append(dict(column=col,h=h,normwise_relative=float(np.max(abs(f-j))/np.max(abs(j)))))
    ulp=rows(d/'ulp_probes.csv');noise=[]
    for comp in [0,1]:
        for sign in [-1,1]:
            rr=[x for x in ulp if int(x['component'])==comp and int(x['direction'])==sign];delta=np.array([float(x['residual'])-float(x['baseline']) for x in rr]);noise.append(dict(component=comp,direction=sign,water_max_delta=float(max(abs(delta[::2]))),salt_max_delta=float(max(abs(delta[1::2])))))
    sv=np.linalg.svd(scaled,compute_uv=False)
    out=dict(probe_is_nearby_ASCII_state_not_original_failed_iterate=True,dimension=len(u),raw_condition_2=float(np.linalg.cond(a)),row_column_scaled_condition_2=float(sv[0]/sv[-1]),scaled_min_singular=float(sv[-1]),linear_eta_component_recomputed=float(max(abs(e)/np.maximum(den,1e-300))),linear_absolute_residual_recomputed=float(max(abs(e))),numpy_longdouble_mantissa_bits=int(np.finfo(np.longdouble).nmant),water_row_scale_minmax=[float(min(row[2::2])),float(max(row[2::2]))],salt_row_scale_minmax=[float(min(row[3::2])),float(max(row[3::2]))],scaled_nonlinear=float(max(abs(r)/row)),fd_columns=errors,ulp_probe=noise,epsilon=float(np.finfo(float).eps),fd_scope='8 columns x 6 central differences; every reduced row; not exhaustive Jacobian validation')
    out['best_fd_relative_by_column']={str(col):min(x['normwise_relative'] for x in errors if x['column']==col) for col in sorted({x['column'] for x in errors})}
    out['FD_scale_study_pass']=max(out['best_fd_relative_by_column'].values())<1e-5
    c.save(c.OUT/'equation_diagnostics.json',out);return out

def main():
    results={n:audit(n) for n in v.CASES if (c.OUT/n/'archive_receipt.json').exists()};eq=equations()
    for n,r in results.items():print(n,r['passed'],[k for k,b in r['checks'].items() if not b])
    c.save(c.OUT/'regression_summary.json',dict(all_regression_checks_passed=all(r['passed'] for r in results.values()),runs={n:dict(passed=r['passed'],failed_solver=r['failed_solver'],metrics=r['metrics']) for n,r in results.items()}))
    print(json.dumps(eq,indent=2))
if __name__=='__main__':main()
