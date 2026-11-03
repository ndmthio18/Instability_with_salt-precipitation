"""Offline WP6 registration only. No SSH, compiler, scheduler, or solver calls."""
import csv,hashlib,json,math
from pathlib import Path
import numpy as np
from single_wavelength_post import params
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'analysis/wp6_bounded_nonlinear_plan_20261103'
WP5=ROOT/'analysis/single_wavelength_validation_20261103'
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def writecsv(name,rows):
    with (OUT/name).open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
def main():
    OUT.mkdir(exist_ok=True,parents=True)
    for p,h in json.loads((WP5/'reference_hashes.json').read_text()).items():assert sha(ROOT/p)==h
    p=params(WP5/'runs/R24-M1-full/params.input');T=float(p['TimeLoop.TEnd']);dt=float(p['TimeLoop.MaxTimeStepSize'])
    cache=ROOT/'tmp/single_wavelength_validation_20261103/reference/M1/post/profiles.npz'
    z=np.load(cache)['z_m'];mid=(z[:-1]+z[1:])/2;H=z[-1]-z[0]
    # Odd subdivision preserves every original BOX face midpoint exactly.
    fine=np.r_[np.concatenate([np.linspace(a,b,4)[:-1] for a,b in zip(z[:-1],z[1:])]),z[-1]]
    fm=(fine[:-1]+fine[1:])/2
    cuts=[]
    specs=[('legacy_'+str(q),H*(1-q),'primary') for q in [.9,.75,.5]]+[(f'shallow_{d:g}mm',d/1000,'secondary_preregistered') for d in [1,2,5,10]]
    for label,depth,role in specs:
        j=int(np.argmin(abs(mid-(z[-1]-depth))));zz=float(mid[j]);jf=int(np.argmin(abs(fm-zz)));assert abs(fm[jf]-zz)<1e-13
        cuts.append(dict(id=label,role=role,requested_depth_mm=depth*1000,actual_z_m=zz,actual_z_over_H=zz/H,actual_depth_mm=(z[-1]-zz)*1000,base_lower_vertex=j,refined_lower_vertex=jf,nested_face_mismatch_m=float(fm[jf]-zz)))
    assert len(set(r['actual_z_m'] for r in cuts))==7
    writecsv('observation_planes.csv',cuts)
    writecsv('vertical_grid_200_vertices.csv',[dict(index=i,z_m=float(v)) for i,v in enumerate(z)])
    writecsv('vertical_grid_600_nested_vertices.csv',[dict(index=i,z_m=float(v)) for i,v in enumerate(fine)])
    stages=[dict(id='bridge',multiplier=1,pairs=1,physical_end_s=T,physical_days=T/86400,steps_per_leg=600,entry='isolated diagnostic build and static tests pass')]
    for k in [2,4,8]:stages.append(dict(id=f'N28_{k}T',multiplier=k,pairs=1,physical_end_s=k*T,physical_days=k*T/86400,steps_per_leg=600*k,entry='previous stage accepted; no threshold; no unresolved validation trigger; budget remaining'))
    writecsv('stage_schedule.csv',stages)
    # Worst of measured R24 M0/M1 full wall times, rounded upward.
    perf=list(csv.DictReader((WP5/'performance_all_runs.csv').open()))
    base=max(float(r['wall_s']) for r in perf if r['run'] in ['R24-M0-full','R24-M1-full']);assert base<82
    cost=[]
    for s in stages:
        sec=2*s['multiplier']*82
        cost.append(dict(task=s['id'],runs=2,reference_T_units=s['multiplier'],nominal_core_h=sec/3600,planning_core_h=5*sec/3600,formula='2 legs * multiplier * 82 s; planning factor 5',status='not_started'))
    for name,factor in [('trigger_dt_half',2),('trigger_vertical_nested600',6),('trigger_R48',2.4)]:
        sec=2*8*82*factor;cost.append(dict(task=name,runs=2,reference_T_units=8,nominal_core_h=sec/3600,planning_core_h=5*sec/3600,formula=f'worst trigger at 8T; per-leg cost factor {factor}; planning factor 5',status='conditional_not_started'))
    cost.append(dict(task='full_domain_audit',runs=2,reference_T_units='signal_window <= 8T',nominal_core_h='',planning_core_h=24,formula='12 h per leg hard cap; admission requires forecast <= 9.6 h/leg',status='conditional_budget_admission'))
    cost.append(dict(task='nx60_short_package',runs='4 + at most 2 full-domain 10-step timing controls',reference_T_units=.5,nominal_core_h='',planning_core_h=4,formula='unmeasured mode-specific throughput; reserved cap, not ETA',status='prepare_after_canonical_review'))
    writecsv('cost_budget.csv',cost)
    legacy=ROOT/'analysis/dense_history_20261102/evidence/runs/r06_c36877_n60';old=params(next(legacy.glob('*native.input')))
    changes={k:dict(nx28=p[k],nx60=old.get(k)) for k in p if p[k]!=old.get(k)}
    (OUT/'nx60_input_differences.json').write_text(json.dumps(changes,indent=2)+'\n')
    reg=dict(schema_version=1,status='plan_only_no_submission',wp5_reference_commit='b626428e832fa4f23fc08ff962edb1935646a63f',T_s=T,dt_s=dt,canonical_nx=28,width_m=float(p['Grid.Positions0'].split()[-1]),amplitude=1e-6,concentration_threshold=.1,max_window_T=8,start_from_zero=True,restart_qualified=False,per_run_wall_limit_h=12,admission_forecast_h=9.6,max_simultaneous_ranks=2,total_core_h_cap=56,reduced_bucket_core_h=28,full_domain_bucket_core_h=24,nx60_bucket_core_h=4,max_successful_solver_runs=22,max_successful_solver_runs_scope='2 bridge + 6 extension + 6 sensitivity + 2 full-domain + 4 nx60 short + 2 timing; each conditional, no bulk submission',observation_planes=cuts,thresholds=dict(salt_closure_absolute_m=1e-14,water_closure_absolute_m=1e-11,initial_pair_inventory_absolute_m=1e-15,candidate_pair_C_floor_m=1e-13,candidate_pair_C_to_residual=10,candidate_persistence_dt_units=10,candidate_sign_fraction=.8,pair_C_required_to_empirical_error=5,pair_dC_relative_stability=.2,mean_inventory_profile_relative=.01,mode_relative=.05,odd_neumann_energy_fraction=1e-6,endpoint_difference_to_current_perturbation=.001,threshold_bracket_relative=.005,near_threshold_concentration=.09),unknown_errors=['U_dt','U_z','U_x','U_domain','event_localization'],t_eff=None,preregistration_change_policy='amend before inspecting dependent new results; preserve original registration; never replace negative deep-cut results with positive shallow cuts')
    (OUT/'preregistration.json').write_text(json.dumps(reg,indent=2)+'\n')
    reg['full_domain_control']=dict(width_m=.6000000000000001,cells_x=672,cells_z=200,nx=28,same_dx_as_reduced=True,forecast_s_per_step=6.4653*672/600,forecast_safety_factor=1.5,must_cover_signal_window=True)
    reg['nx60_preparation']=dict(width_m=.6000000000000001/60,cells_x=[20,40],cells_z=200,physical_nx=60,short_end_s=T/2,dt_s=dt,profile='historical A60 with source hash',canonical_replacement=False,full_reference_gate='pending unless same-ELF common-window validation is feasible',scope='separate initialization-protocol sensitivity, not pure wavenumber causality')
    reg['thresholds']['nonlinear_structure_review_harmonic_fraction']=.01
    reg['thresholds']['nonlinear_structure_review_relative_sigma']=.01
    reg['thresholds']['endpoint_absolute_floor']=1e-12
    (OUT/'preregistration.json').write_text(json.dumps(reg,indent=2)+'\n')
    sources=[ROOT/'command/6.md',WP5/'preregistration.json',WP5/'full_validation.json',WP5/'performance_all_runs.csv',WP5/'runs/R24-M1-full/params.input',WP5/'runs/R24-M1-full/eigenprofile.dat',WP5/'reference_hashes.json',cache,ROOT/'model/transport_v1/source/appl/simpleRichardsnc/transportdiagnostics.hh',next(legacy.glob('*native.input')),next(legacy.glob('*eigenprofile.dat')),legacy/'audit.json']
    (OUT/'source_hashes.json').write_text(json.dumps({str(q.relative_to(ROOT)).replace('\\','/'):sha(q) for q in sources},indent=2)+'\n')
    checks=dict(reference_source_hashes_unchanged=True,distinct_planes=len(cuts),nested_grid_vertices=len(fine),all_cut_coordinates_preserved=True,amplitude_unchanged=True,nx28_physical_k=2*math.pi/float(p['Grid.Positions0'].split()[-1]),nx60_physical_k=2*math.pi*60/.6000000000000001,solver_started=False,source_modified=False)
    (OUT/'plan_checks.json').write_text(json.dumps(checks,indent=2)+'\n')
    manifest={str(q.relative_to(ROOT)).replace('\\','/'):dict(sha256=sha(q),bytes=q.stat().st_size) for q in sorted([q for q in OUT.iterdir() if q.is_file() and q.name!='artifact_manifest.json']+[Path(__file__)])}
    (OUT/'artifact_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(dict(checks=checks,planes=cuts,costs=cost),indent=2))
if __name__=='__main__':main()
