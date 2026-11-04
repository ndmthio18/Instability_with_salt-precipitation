"""WP9: physical acceptance of an R96 leg against its frozen protocol.

Usage: python3 wp9_validate_r96.py [M1|M0]

Every gate is the one registered in the leg's protocol before the run. No gate is
loosened here; a failing gate is reported, not adjusted. No solver is invoked.
"""
import collections,json,math,sys
from pathlib import Path
import numpy as np
import wp9_common as c

W7=c.ROOT/'analysis/wp7_shortwave_validation_20261103/runs'
LEGS={
 'M1':dict(name='WP9-R96-M1-4T',host='sal',base=c.OUT,run=c.OUT/'run',
           protocol=c.OUT/'r96_protocol.json',auth=c.OUT/'user_authorization.json',
           archive=c.OUT/'archive_receipt.json',reduction=c.OUT/'reduction_receipt.json',
           pairing=c.OUT/'frame_pairing_receipt.json',
           spatial=c.OUT/'fields/spatial_metrics.csv',
           npz=c.CACHE/'R96_M1_profiles_spectra.npz',reference=W7/'WP7-R48-M1-4T',
           out=c.OUT/'r96_validation.json'),
 'M0':dict(name='WP9-R96-M0-4T',host='hpcc',base=c.OUT/'r96_m0',run=c.OUT/'r96_m0/run',
           protocol=c.OUT/'r96_m0/protocol.json',auth=c.OUT/'r96_m0/user_authorization.json',
           archive=c.OUT/'r96_m0/archive_receipt.json',reduction=c.OUT/'r96_m0/reduction_receipt.json',
           pairing=c.OUT/'r96_m0/frame_pairing_receipt.json',
           spatial=c.OUT/'r96_m0/fields/spatial_metrics.csv',
           npz=c.CACHE/'R96_M0_profiles_spectra.npz',reference=W7/'WP7-R48-M0-4T',
           out=c.OUT/'r96_m0/validation.json'),
}
SALT_RESIDUAL=1e-14
WATER_RESIDUAL=1e-11
OK_STOPS=('requested_horizon','surface_solubility_threshold','interior_solubility_threshold')


def regions(path):
    out=collections.defaultdict(list)
    for r in c.rows(path):out[r['region']].append(r)
    return out


def main(leg='M1'):
    L=LEGS[leg];RUN=L['run'];NAME=L['name'];R48=L['reference']
    proto=json.loads(L['protocol'].read_text(encoding='utf-8'))
    plan=json.loads((RUN/'plan.json').read_text())
    fin=json.loads((RUN/'finished.json').read_text())
    rec=json.loads(L['archive'].read_text())
    red=json.loads(L['reduction'].read_text())
    ev=c.rows(RUN/'surface_events.csv')
    bud=regions(RUN/'transport_budget.csv')
    dt=proto['numerics']['dt_s'];steps=proto['numerics']['steps']
    solver=fin['solver'];stop=solver['stop']
    checks={};met={}

    # 1 completion
    checks['exit_zero']=solver['returncode']==0
    checks['stop_reason_admitted']=stop['reason'] in OK_STOPS
    met['stop_reason']=stop['reason'];met['accepted_steps']=stop['accepted_steps']
    threshold_stop=stop['reason']!='requested_horizon'
    checks['accepted_steps']=(stop['accepted_steps']==steps) if not threshold_stop else stop['accepted_steps']<steps
    met['last_accepted_time_s']=stop['last_accepted_time_s']
    met['last_accepted_over_T']=stop['last_accepted_time_s']/c.T_S

    # 2 fixed timestep, the exact accumulation, and the end-of-window truncation
    #
    # The time loop never steps past the requested end. At the campaign step the
    # accumulated axis lands just below the decimal TEnd, so R24, R48 and the strict
    # background take 4800 full steps and stop short of it. At dt/2 the accumulation
    # overshoots TEnd, so the solver shortens the last step to land on it exactly. That
    # truncation is a property of the run and is measured here, not hidden; WP6 recorded
    # the same behaviour for its own coarse leg.
    t=np.array([float(r['time_s']) for r in ev])
    d=np.array([float(r['dt_s']) for r in ev])
    interior=d[1:-1] if len(d)>2 else d[1:]
    met['distinct_accepted_dt']=sorted({float(x) for x in d[1:]})[:4]
    met['final_step_dt_s']=float(d[-1])
    met['final_step_truncation_s']=float(dt-d[-1])
    checks['fixed_dt_on_every_step_but_the_last']=bool(np.all(interior==dt))
    checks['final_step_not_longer_than_dt']=bool(d[-1]<=dt)
    budget=steps*math.ulp(c.HORIZON_S)
    met['time_axis_budget_s']=budget
    checks['final_step_truncation_within_budget']=bool(0.0<=dt-d[-1]<=budget)
    # the accumulation is exact up to the step before the last one
    acc=0.0;drift=0.0
    for i,v in enumerate(t[:-1]):
        if i:acc+=dt
        drift=max(drift,abs(float(v)-acc))
    met['fixed_dt_accumulation_drift_s']=drift
    checks['exact_fixed_dt_accumulation']=drift==0.0
    met['endpoint_error_s']=float(abs(t[-1]-c.HORIZON_S))
    checks['endpoint_on_the_requested_horizon']=bool(
        abs(t[-1]-c.HORIZON_S)<=budget) if not threshold_stop else True
    checks['no_timestep_reduction']=not [s for s in solver['failure_lines']
                                         if 'reduc' in s.lower() and 'RetryTimeStepReductionFactor' not in s]

    # 3 frames and archive
    met['frames']=fin['frames'] if 'frames' in fin else rec['frames']
    expected_frames=proto['numerics']['expected_frames']
    checks['frame_count']=rec['frames']==expected_frames or threshold_stop
    checks['archive_hash_matches']=rec['archive_matches_finished']
    checks['all_frame_hashes_verified']=bool(red['reduced_audit']['all_field_hashes_verified']
                                             and red['spatial_audit']['all_field_hashes_verified'])
    # the two reviewed field passes ran against the symlink view, whose manifest carries the
    # same sha256 values with retargeted paths; both ends of that mapping are checked here
    pair=json.loads(L['pairing'].read_text(encoding='utf-8'))
    checks['field_manifest_consistent']=(red['reduced_audit']['field_manifest_sha256']
                                         ==red['spatial_audit']['field_manifest_sha256']
                                         ==pair['view_manifest_sha256'])
    checks['frame_view_preserves_hashes']=bool(pair['sha256_values_unchanged']
                                               and pair['source_manifest_sha256']==rec['field_manifest_sha256'])
    met['frame_view']=dict(stride=pair['stride'],frames=pair['frames'],
                           recorded_steps=pair['recorded_steps'],
                           first_kept_step=pair['first_kept_step'],last_kept_step=pair['last_kept_step'])
    checks['frame_view_covers_every_even_step']=(pair['recorded_steps']==steps+1
                                                 and pair['frames']==steps//pair['stride']+1) or threshold_stop

    # 4 conservative closure, every region and every accepted step
    worst_s=worst_w=0.0;identity=0.0
    for region,rows in bud.items():
        s=np.array([float(r['salt_residual_m']) for r in rows])
        w=np.array([float(r['water_residual_m']) for r in rows])
        worst_s=max(worst_s,float(np.max(np.abs(s))));worst_w=max(worst_w,float(np.max(np.abs(w))))
        I=np.array([float(r['inventory_salt_m']) for r in rows])
        cum=np.array([float(r['cumulative_cut_salt_m']) for r in rows])
        identity=max(identity,float(np.max(np.abs((I-I[0])+(cum-cum[0])-s))))
    met['max_salt_residual_m']=worst_s;met['max_water_residual_m']=worst_w
    met['max_budget_identity_error_m']=identity
    checks['salt_closure']=worst_s<=SALT_RESIDUAL
    checks['water_closure']=worst_w<=WATER_RESIDUAL
    checks['budget_identity']=identity<=SALT_RESIDUAL

    # 5 independent reconstruction from the hash-verified fields
    ind=rec['independent_checks']
    met['independent_vtu']=dict(frames=ind['frames'],recorded_steps=ind['recorded_steps'],
        frame_stride=ind['frame_stride'],
        max_inventory_salt_difference_m=ind['max_inventory_salt_difference_m'],
        max_inventory_water_difference_m=ind['max_inventory_water_difference_m'],
        max_online_absolute_error=max(ind['online_absolute_errors'].values()),
        max_salt_closure_m=max(v['salt_closure_m'] for v in ind['independent_budgets'].values()),
        max_water_closure_m=max(v['water_closure_m'] for v in ind['independent_budgets'].values()),
        max_cumulative_salt_reintegration_error=max(v['cumulative_salt_reintegration_error']
                                                    for v in ind['independent_budgets'].values()))
    checks['independent_vtu_quadrature']=bool(ind['passed'])
    with np.load(L['npz']) as f:
        zs=np.asarray(f['z_m']);salt=np.asarray(f['salt_mean'])
    wz=np.diff(np.r_[zs[0],(zs[1:]+zs[:-1])/2,zs[-1]])
    inv=salt@wz
    whole=np.array([float(r['whole_salt_inventory_m']) for r in bud['free_domain']])
    stride=ind['frame_stride']
    met['independent_inventory_max_difference_m']=float(np.max(np.abs(inv-whole[::stride][:len(inv)])))
    checks['independent_storage']=met['independent_inventory_max_difference_m']<=SALT_RESIDUAL

    # 6 initial state and finiteness
    met['initial_surface_mean']=float(ev[0]['x_top_area_mean'])
    ev48=c.rows(R48/'surface_events.csv')
    met['initial_surface_mean_difference_to_R48']=abs(float(ev[0]['x_top_area_mean'])
                                                      -float(ev48[0]['x_top_area_mean']))
    checks['initial_state_matches_R48']=met['initial_surface_mean_difference_to_R48']<=1e-12
    finite=all(np.isfinite([float(r[k]) for r in ev for k in
                            ['x_top_area_mean','x_top_max','x_domain_max','sigma_top_node']]))
    checks['all_finite']=bool(finite)

    # 7 threshold status
    mx=max(float(r['x_domain_max']) for r in ev)
    met['max_local_x']=mx;met['fraction_of_solubility']=mx/c.SOLUBILITY
    met['threshold_reached']=any(r['surface_threshold_reached']=='1' or r['global_threshold_reached']=='1' for r in ev)
    checks['threshold_status_consistent']=(met['threshold_reached']==threshold_stop)

    # 8 shared time axis with R48
    t48=np.array([float(r['time_s']) for r in ev48])
    n=min(len(t48),len(t[::stride]))
    dd=np.abs(t[::stride][:n]-t48[:n])
    met['paired_states']=int(n);met['max_paired_time_difference_s']=float(dd.max())
    met['paired_time_budget_s']=proto['numerics']['time_axis']['tolerance_s']
    checks['paired_time_axis_within_budget']=bool(dd.max()<=met['paired_time_budget_s'])

    # 9 lateral symmetry screen, retained descriptively
    sp=c.rows(L['spatial'])
    met['max_mirror_difference']=max(float(r['max_mirror_difference']) for r in sp)
    met['max_endpoint_difference']=max(float(r['max_endpoint_difference']) for r in sp)
    met['max_relative_anomaly_final']=float(sp[-1]['max_relative_anomaly'])
    checks['symmetry_screen_recorded']=True

    # 10 resources
    met['solver_wall_s']=solver['wall_s'];met['peak_rss_kib']=solver['peak_rss_kib']
    met['newton_iterations']=solver['newton_iterations']
    met['newton_per_step']=solver['newton_iterations']/max(stop['accepted_steps'],1)
    met['archive_bytes']=rec['bytes']
    wall_cap=(proto.get('resources') or proto['execution'])['hard_wall_s']
    met['hard_wall_s']=wall_cap
    checks['within_hard_wall']=solver['wall_s']<=wall_cap
    checks['one_newton_pass_per_step']=solver['newton_iterations']>=stop['accepted_steps']

    # 11 provenance
    checks['frozen_input_used']=plan['input_sha256']==c.sha(RUN/'params.input')
    checks['elf_is_the_R24_R48_executable']=plan['elf_sha256']==proto['numerics']['executable']['sha256']
    checks['authorized']=c.sha(L['auth'])==plan['user_authorization_sha256']
    checks['protocol_unchanged']=c.sha(L['protocol'])==plan['preregistration_sha256']

    met['host']=L['host'];met['threads']=solver.get('threads',1)
    out=dict(leg=leg,run=str(RUN.relative_to(c.ROOT)),name=NAME,host=L['host'],
             passed=all(checks.values()),checks=checks,metrics=met,
             failures=[k for k,v in checks.items() if not v],
             censored_by_threshold=bool(met['threshold_reached']),
             gates=dict(salt_residual_m=SALT_RESIDUAL,water_residual_m=WATER_RESIDUAL,
                        source='the WP7 thresholds, registered unchanged in the protocol'))
    c.save(L['out'],out)
    print(json.dumps(dict(passed=out['passed'],failures=out['failures'],
                          metrics={k:met[k] for k in ['stop_reason','accepted_steps','frames',
                              'max_salt_residual_m','max_water_residual_m','max_local_x',
                              'solver_wall_s','newton_per_step','max_paired_time_difference_s']
                              if k in met}),indent=2,default=str))
    assert out['passed'],out['failures']


if __name__=='__main__':main(sys.argv[1] if len(sys.argv)>1 else 'M1')
