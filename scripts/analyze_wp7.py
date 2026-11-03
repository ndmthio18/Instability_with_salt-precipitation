"""Assemble the bounded WP7 evidence without inventing a missing R48 pair."""
import csv,json,math
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'analysis/wp7_shortwave_validation_20261103'
WP6=ROOT/'analysis/wp6_4T_exploratory_20261103/runs'
WP6TMP=ROOT/'tmp/wp6_4T_exploratory_20261103/runs'
def load(path):return json.loads(path.read_text())
def rows(path):
    with path.open() as f:return list(csv.DictReader(f))
def keyed(data,maximum_step=None):
    return {(int(r['step']),r['region']):r for r in data if maximum_step is None or int(r['step'])<=maximum_step}
def main():
    one=OUT/'runs/WP7-1D-background-T';r48=OUT/'runs/WP7-R48-M1-bench10'
    p1=np.load(one/'post/profiles.npz');p2=np.load(WP6TMP/'WP6-M0-explore4T/post/profiles.npz')
    exact=np.where(p1['time_s']==p2['time_s'][:len(p1['time_s'])])[0];assert np.array_equal(exact,np.arange(1200))
    profile_delta=p1['x_mean'][:1200]-p2['x_mean'][:1200]
    profile_scale=float(np.max(np.abs(p2['x_mean'][:1200])))
    a=keyed(rows(one/'transport_budget.csv'),1199);b=keyed(rows(WP6/'WP6-M0-explore4T/transport_budget.csv'),1199);assert a.keys()==b.keys()
    fields=['inventory_salt_m','inventory_water_m','cut_salt_adv_down_m_s','cut_salt_diff_down_m_s','cut_salt_total_down_m_s','cut_water_down_m_s','cumulative_cut_salt_m','cumulative_net_water_out_m','salt_residual_m','water_residual_m']
    comparisons={}
    for field in fields:
        delta=np.array([float(a[k][field])-float(b[k][field]) for k in a]);scale=max(abs(float(b[k][field])) for k in b)
        comparisons[field]=dict(max_absolute_difference=float(np.max(np.abs(delta))),reference_max_absolute=scale,global_normalized_difference=float(np.max(np.abs(delta))/max(scale,1e-30)))
    sigma1=float(np.max(np.abs(p1['x_sigma'])));range1=2*sigma1
    bg_checks=load(one/'post/independent_checks.json');profile_abs=float(np.max(np.abs(profile_delta)))
    background=dict(
      interpretation='one-cell horizontal extrusion of the same 2D equations, retained separately from 2D M0',
      exact_common_frames=1200,last_exact_common_time_s=float(p1['time_s'][1199]),terminal_time_difference_s=float(p1['time_s'][-1]-p2['time_s'][1200]),
      profile_max_absolute_difference=profile_abs,profile_global_normalized_difference=profile_abs/profile_scale,
      one_cell_max_lateral_sigma=sigma1,one_cell_max_lateral_range=range1,lateral_range_method='two equal-width endpoint vertices: range=2*sigma',
      lateral_uniformity_gate=range1<=1e-14,profile_gate=profile_abs/profile_scale<=1e-6,
      budget_comparisons=comparisons,independent_storage_and_budget_checks=bg_checks,
      closure_gate=bg_checks['passed'],overall_preregistered_gate=range1<=1e-14 and profile_abs/profile_scale<=1e-6 and bg_checks['passed'])
    r48_budget=rows(r48/'transport_budget.csv');old_budget=rows(WP6/'WP6-M1-explore4T/transport_budget.csv')
    r0={x['region']:x for x in r48_budget if x['step']=='0'};o0={x['region']:x for x in old_budget if x['step']=='0'};assert r0.keys()==o0.keys()
    initial={}
    for field in ['inventory_salt_m','inventory_water_m']:
        diffs=[abs(float(r0[k][field])-float(o0[k][field])) for k in r0];scale=max(abs(float(o0[k][field])) for k in o0)
        initial[field]=dict(max_absolute_difference=max(diffs),global_normalized_difference=max(diffs)/scale)
    admission=load(OUT/'benchmark_admission.json');audit=load(r48/'post/field_audit.json');checks=load(r48/'post/independent_checks.json')
    W=0.6/28
    r48_summary=dict(status='benchmark_only; 4T pair not launched',benchmark=admission,initial_inventory_vs_R24=initial,field_audit=audit,independent_checks=checks,
      width_m=W,dx_m=W/48,canonical_cells_per_wavelength=48,initial_physical_wavenumber_rad_m=2*math.pi/W,
      R24_late_m10_wavelength_m=2*W/10,R48_cells_per_R24_late_m10_wavelength=2*48/10,
      R48_cells_per_m12_wavelength=2*48/12,three_to_four_T_comparison='not available: long pair blocked by durable-storage gate')
    local=load(OUT/'local_flux_smoke/bridge_result.json');build=load(OUT/'build/attempt_01/result.json')
    local_summary=dict(static_checks=load(OUT/'local_flux_static_checks.json'),build=dict(configure_exit=build['configure_exit'],build_exit=build['build_exit'],executable=build['executable']),bridge=local,
      enabled_metrics=['local up/down water and advective salt flux','down-flow area fraction','area-weighted concentration-flow covariance','horizontal flux-concentration association','depth-plane penetration'],
      evidence_boundary='The one-step bridge validates conservative identities only. No late-time local-flow mechanism was measured.')
    decision=dict(scope='command/11 only',R48_4T_pair='NO-GO: cost admitted but registered durable-storage gate failed',R96='NO-GO',multiwave='NO-GO',later_8T='NO-GO',
      background_model='partially verified: profiles, budgets and closure pass; strict 1e-14 lateral-range gate fails',
      R24_R48_late_wavelength='unresolved because R48 3T-4T data do not exist',local_scvf='implementation compiled and conservative bridge passed; late-time mechanism pending',
      t_eff=None,crystallization_delay_label=None,training_labels_created=False,next_gate='provide >=200 GiB durable writable storage with >=20% filesystem free, then run the frozen R48 M0/M1 4T pair and archive each leg immediately')
    (OUT/'background_validation.json').write_text(json.dumps(background,indent=2)+'\n')
    (OUT/'r48_preflight_summary.json').write_text(json.dumps(r48_summary,indent=2)+'\n')
    (OUT/'local_scvf_summary.json').write_text(json.dumps(local_summary,indent=2)+'\n')
    (OUT/'decision.json').write_text(json.dumps(decision,indent=2)+'\n')
    with (OUT/'background_budget_comparison.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['quantity','max_absolute_difference','reference_max_absolute','global_normalized_difference'])
        for k,v in comparisons.items():w.writerow([k,v['max_absolute_difference'],v['reference_max_absolute'],v['global_normalized_difference']])
    with (OUT/'runtime_cost.csv').open('w',newline='') as f:
        w=csv.writer(f);w.writerow(['run','solver_wall_s','cpu_s','post_archive_wall_s','accepted_steps','archive_bytes'])
        for d in [r48,one]:
            x=load(d/'finished.json');w.writerow([d.name,x['solver']['wall_s'],x['solver']['cpu_s'],x['post_archive_wall_s'],x['solver']['stop']['accepted_steps'],x['archive']['bytes']])
    print(json.dumps(decision,indent=2))
if __name__=='__main__':main()
