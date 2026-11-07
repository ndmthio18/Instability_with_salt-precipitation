import csv
import json
import pathlib
import sys

root=pathlib.Path.cwd()
stage=sys.argv[1]
folder=root/('stage_'+stage)
def read(path):
    with path.open() as f:
        return [{k:float(v) for k,v in row.items()} for row in csv.DictReader(f)]

summary=[]
for path in sorted(folder.iterdir()):
    if not path.is_dir(): continue
    history=path/'coupled_history.csv'
    if not history.exists():
        summary.append(dict(case=path.name,status='no_accepted_history'))
        continue
    rows=read(history); first,last=rows[0],rows[-1]
    assert first['time_s']==0 and first['active']==0 and first['top_solid_kg_m2']==0
    for a,b in zip(rows,rows[1:]):
        assert b['time_s']>a['time_s']
        assert abs(b['time_s']-a['time_s']-b['dt_s'])<1e-9
        assert b['active']>=a['active']
    salt=max(abs(r['salt_residual_kg']) for r in rows)/first['soil_salt_kg']
    water=max(abs(r['water_residual_kg']) for r in rows)/first['soil_water_kg']
    solid=max(abs(r['solid_ledger_residual_kg']) for r in rows)/first['soil_salt_kg']
    assert max(salt,water,solid)<1e-8,(path.name,salt,water,solid)
    active=[r for r in rows if r['active']]
    activation=active[0] if active else None
    exchange_salt=exchange_water=None
    first_precip=next((r['time_s'] for r in rows if r['top_solid_kg_m2']>1e-15),None)
    if activation:
        assert activation['pore_salt_transfer_kg']>0 and activation['pore_water_transfer_kg']>0
        assert activation['top_soil_sw']>=.99
        assert all(r['pore_salt_transfer_kg']==activation['pore_salt_transfer_kg'] and
                   r['pore_water_transfer_kg']==activation['pore_water_transfer_kg'] for r in active)
        assert abs(activation['crust_salt_kg']-activation['pore_salt_transfer_kg'])<1e-12*first['soil_salt_kg']
        assert abs(activation['crust_water_kg']-activation['pore_water_transfer_kg'])<1e-12*first['soil_water_kg']
        if 'interface_liquid_kg_m2' in activation:
            exchange_salt=max(abs(r['soil_salt_kg']-activation['soil_salt_kg']+.001*r['interface_salt_kg_m2']) for r in active)
            exchange_water=max(abs(r['soil_water_kg']-activation['soil_water_kg']+
                .001*(r['interface_liquid_kg_m2']-r['interface_salt_kg_m2'])) for r in active)
            assert exchange_salt<1e-8*first['soil_salt_kg'] and exchange_water<1e-8*first['soil_water_kg']
    code=int((path/'exit_status.txt').read_text())
    events=(path/'coupling_events.csv').read_text()
    event_rows=list(csv.DictReader(events.splitlines()))
    attempts=[r for r in event_rows if r['event'].startswith('activation_')]
    blocked=[r for r in event_rows if r['event']=='activation_blocked_saturation']
    item=dict(case=path.name,exit_code=code,accepted_time_s=last['time_s'],
              status=('activation_complete' if stage=='b' and activation else
                      'horizon_complete' if code==0 and last['time_s']==2500 else 'stopped_prefix'),
              precipitation_s=first_precip,activation_s=activation['time_s'] if activation else None,
              activation_height_m=activation['h_m'] if activation else None,
              first_activation_attempt_s=float(attempts[0]['time_s']) if attempts else None,
              first_refused_donor_sw=float(blocked[0]['value1']) if blocked else None,
              subcrystal_activation=bool(activation and activation['h_m']<3.5e-6),
              activation_blocked_saturation='activation_blocked_saturation' in events,
              salt_residual_relative=salt,water_residual_relative=water,solid_ledger_relative=solid,
              interface_salt_residual_relative=exchange_salt/first['soil_salt_kg'] if exchange_salt is not None else None,
              interface_water_residual_relative=exchange_water/first['soil_water_kg'] if exchange_water is not None else None,
              final_height_m=last['h_m'],final_h_mass_m=last['h_mass_m'],
              final_top_solid_kg_m2=last['top_solid_kg_m2'],final_internal_solid_kg_m2=last['internal_solid_kg_m2'],
              final_epsilon_min=last['epsilon_min'],final_evaporation=last['evaporation_kg_m2_s'],
              final_soil_top_mass_fraction=last['top_soil_mass_fraction'],
              final_soil_top_sw=last['top_soil_sw'],coupled_trials=last['coupled_trials'],rejected_trials=last['rejected_trials'],
              rejected_trial_injected='injected rollback test' in events)
    # Actual common-time comparison, without extrapolating the 1000-s P2-G run.
    if stage=='c':
        old=read(root.parent/'p2_g_surface_precipitation_20261107/runs/n80/surface_precipitation.csv')
        old_by_time={r['time_s']:r for r in old}
        prefix=[r for r in rows if not r['active'] and r['time_s'] in old_by_time]
        prefix_difference=max(abs(r['soil_salt_kg']-old_by_time[r['time_s']]['soil_salt_kg']) for r in prefix)
        item['precrust_soil_salt_difference_kg']=prefix_difference
        if path.name=='nominal': assert prefix_difference<1e-12*first['soil_salt_kg']
        common=[r for r in rows if r['time_s'] in old_by_time]
        if common:
            r=common[-1]; g=old_by_time[r['time_s']]
            item['p2g_comparison']=dict(time_s=r['time_s'],ale_height_m=r['h_m'],ale_h_mass_m=r['h_mass_m'],
                p2g_h_mass_m=g['h_mass_m'],ale_prat_prediction_m=r['h_prat_m'],p2g_prat_prediction_m=g['h_prat_m'],
                soil_concentration_difference=r['top_soil_mass_fraction']-g['top_mass_fraction'],
                evaporation_relative_difference=r['evaporation_kg_m2_s']/g['evaporation_kg_m2_s']-1)
    summary.append(item)

result=dict(stage=stage,cases=summary)
if stage=='b':
    by={r['case']:r for r in summary}
    assert by['activation_1um']['status']=='activation_complete'
    assert by['activation_rollback']['status']=='activation_complete'
    assert by['activation_rollback']['rejected_trial_injected']
    assert by['activation_rollback']['rejected_trials']>=1
    assert all(r.get('exit_code')==0 for r in summary)
    result['implementation_gate']='PASS'
else:
    by={r['case']:r for r in summary}
    if by['nominal']['status']=='horizon_complete':
        assert by['coupled_rollback']['status']=='horizon_complete'
        assert by['coupled_rollback']['rejected_trial_injected'] and by['coupled_rollback']['rejected_trials']>=1
        result['sensitivity']={}
        for name,tolerance in [('soil160',.02),('crust128',.02),('dt_half',.01)]:
            r=by[name]; b=by['nominal']
            if r['activation_s'] is None or r['accepted_time_s']!=b['accepted_time_s']:
                result['sensitivity'][name]=dict(pass_gate=False,reason='different activation or accepted time')
            else:
                metrics={k:abs(r[k]-b[k])/max(abs(r[k]),1e-30) for k in
                         ('final_height_m','final_h_mass_m','final_evaporation','final_epsilon_min',
                          'final_top_solid_kg_m2','final_internal_solid_kg_m2','final_soil_top_mass_fraction')}
                result['sensitivity'][name]=dict(pass_gate=max(metrics.values())<tolerance,
                                                  relative_to=name,metrics=metrics,tolerance=tolerance)
    result['physical_validation']=False
(folder/'verification.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result,indent=2))
