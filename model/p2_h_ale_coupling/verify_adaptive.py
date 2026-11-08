"""Audit full accepted histories and allocation profiles, retaining failed gates."""
import csv,hashlib,json,pathlib,statistics,sys
root=pathlib.Path.cwd()
out=root/'audit'; out.mkdir(exist_ok=True)
def read(p):
    if not p.exists(): return []
    with p.open() as f: return [{k:float(v) for k,v in r.items()} for r in csv.DictReader(f)]
def events(p):
    with p.open() as f: return list(csv.DictReader(f))
summary=[]; controlled=[]; errors=[]
frozen_root=root.parent/'p2_h_ale_coupling_20261107/repo'
unchanged={}
frozen_paths=['model/p2_h_ale_coupling/ale.hh','model/p2_h_ale_coupling/bdf.hh',
    'model/p2_g_surface_precipitation/surface.hh','model/p2_d_coupling/cpp/pratcrust.hh',
    'model/p2_f_analytical_growth/integration/p2f_1d.input']
frozen_paths += [str(p.relative_to(root/'repo')) for p in (root/'repo/model/p1_r2_conservative_richards/source/appl').rglob('*') if p.is_file()]
for name in frozen_paths:
    actual=hashlib.sha256((root/'repo'/name).read_bytes()).hexdigest()
    assert actual==hashlib.sha256((frozen_root/name).read_bytes()).hexdigest(),name
    unchanged[name]=actual
for stage,target in [('controlled',controlled),('smoke',summary),('matrix',summary),('refusal',summary)]:
    if not (root/stage/'cases.json').exists(): continue
    for spec in json.loads((root/stage/'cases.json').read_text()):
        path=root/stage/spec['name']; rows=read(path/'coupled_history.csv')
        item=dict(spec,stage=stage,exit_code=int((path/'exit_status.txt').read_text()))
        if not rows: item['status']='no_accepted_history'; errors.append(item); target.append(item); continue
        first,last=rows[0],rows[-1]; active=[r for r in rows if r['active']]; activation=active[0] if active else None
        assert first['active']==0 and first['h_m']==0 and first['top_solid_kg_m2']==0
        assert all(b['time_s']>a['time_s'] and b['active']>=a['active'] for a,b in zip(rows,rows[1:]))
        salt=max(abs(r['salt_residual_kg']) for r in rows)/first['soil_salt_kg']
        water=max(abs(r['water_residual_kg']) for r in rows)/first['soil_water_kg']
        assert max(salt,water)<=1e-8,(spec['name'],salt,water)
        e=events(path/'coupling_events.csv'); attempts=[r for r in e if r['event'].startswith('activation_')]
        accepted=[r for r in e if r['event']=='accepted_activation_attempt']
        assert len(accepted)==int(bool(active))
        item.update(status=('controlled_only' if stage=='controlled' else
            'coupled_horizon' if active and item['exit_code']==0 and last['time_s']==spec.get('tend',2500) else
            'precrust_horizon' if not active and item['exit_code']==0 and last['time_s']==spec.get('tend',2500) else 'stopped_prefix'),
            accepted_time_s=last['time_s'],precipitation_s=next((r['time_s'] for r in rows if r['top_solid_kg_m2']>1e-15),None),
            first_attempt_s=float(attempts[0]['time_s']) if attempts else None,
            activation_s=activation['time_s'] if activation else None,
            salt_residual_relative=salt,water_residual_relative=water,
            rejected_trials=last['rejected_trials'],refusal_reasons=sorted(set(r['event'] for r in e if 'blocked' in r['event'] or 'rejected' in r['event'])),
            final={k:last[k] for k in ('h_m','h_mass_m','top_solid_kg_m2','internal_solid_kg_m2','epsilon_min','evaporation_kg_m2_s','top_soil_sw','top_soil_mass_fraction')})
        steps=[r['dt_s'] for r in rows[1:]]
        item['accepted_dt_statistics']=dict(count=len(steps),minimum=min(steps) if steps else None,
            maximum=max(steps) if steps else None,median=statistics.median(steps) if steps else None,
            newton_reduced_steps=sum(r['retried'] for r in read(path/'accepted_steps.csv')))
        if activation:
            assert all(r['pore_salt_transfer_kg']==activation['pore_salt_transfer_kg'] and r['pore_water_transfer_kg']==activation['pore_water_transfer_kg'] for r in active)
            es=max(abs(r['soil_salt_kg']-activation['soil_salt_kg']+.001*r['interface_salt_kg_m2']) for r in active)/first['soil_salt_kg']
            ew=max(abs(r['soil_water_kg']-activation['soil_water_kg']+.001*(r['interface_liquid_kg_m2']-r['interface_salt_kg_m2'])) for r in active)/first['soil_water_kg']
            assert max(es,ew)<=1e-8
            item.update(interface_salt_relative=es,interface_water_relative=ew)
        da=read(path/'donor_attempts.csv'); profiles=read(path/'donor_profiles.csv')
        max_sr=max_wr=0
        for attempt in da:
            sr=abs(attempt['salt_transfer_residual_kg'])/first['soil_salt_kg']; wr=abs(attempt['water_transfer_residual_kg'])/first['soil_water_kg']
            assert max(sr,wr)<=1e-12
            max_sr=max(max_sr,sr); max_wr=max(max_wr,wr)
        for aid in sorted(set(p['attempt'] for p in profiles)):
            ps=[p for p in profiles if p['attempt']==aid]; a=next(a for a in da if a['attempt']==aid)
            assert len(set(p['dof'] for p in ps))==len(ps)
            assert abs(sum(p['volume_m3'] for p in ps)-.0002)<1e-14
            assert all(p['extracted_kg']==0 or p['box_deepest_m']<=spec['depth_mm']*.001+1e-14 for p in ps)
            for p in ps:
                assert p['extracted_kg']<=p['available_kg']+1e-18
                old=.41*p['volume_m3']*p['sw_before']*p['rho_before']
                new=.41*p['volume_m3']*p['sw_after']*p['rho_after']
                assert abs(old*p['c_before']-new*p['c_after']-p['extracted_kg']*p['c_before'])<=1e-12*first['soil_salt_kg']
                assert abs(old*(1-p['c_before'])-new*(1-p['c_after'])-p['extracted_kg']*(1-p['c_before']))<=1e-12*first['soil_water_kg']
                if p['extracted_kg']>0: assert p['sw_after']>=.99
                else: assert p['p_before_pa']==p['p_after_pa'] and p['c_before']==p['c_after']
            if a['sufficient']:
                assert abs(sum(p['extracted_kg'] for p in ps)-a['required_kg'])<=1e-12*first['soil_water_kg']
                assert abs(sum(p['extracted_kg']*p['c_before'] for p in ps)-a['salt_transfer_kg'])<=1e-12*first['soil_salt_kg']
        if da:
            a=next((a for a in da if accepted and a['attempt']==float(accepted[0]['value1'])),da[0])
            ps=[p for p in profiles if p['attempt']==a['attempt']]
            depths=sorted(set(p['depth_m'] for p in ps)); spacings=[b-a for a,b in zip(depths,depths[1:])]
            item['allocation']=dict(a,pressure_change_max_pa=max((abs(p['p_after_pa']-p['p_before_pa']) for p in ps),default=0),
                saturation_change_max=max((abs(p['sw_after']-p['sw_before']) for p in ps),default=0),
                concentration_change_max=max((abs(p['c_after']-p['c_before']) for p in ps),default=0),
                transfer_salt_max_relative=max_sr,transfer_water_max_relative=max_wr,
                min_soil_spacing_m=min(spacings) if spacings else None,max_soil_spacing_m=max(spacings) if spacings else None,
                top_box_depth_m=max((p['box_deepest_m'] for p in ps if p['depth_m']==0),default=None))
            # Raw attempts retain the greedy proposal's traversal depth even
            # when total capacity is insufficient. Actual transferred depth is
            # zero for a refused allocation; no partial withdrawal was applied.
            item['allocation']['proposed_depth_m']=a['used_depth_m']
            if not a['sufficient']: item['allocation']['used_depth_m']=0
        if spec.get('precrust'):
            item['thresholds_s']={str(h):next((r['time_s'] for r in rows if r['h_mass_m']>=h*1e-6),None) for h in (1,2,5,10)}
        target.append(item)

by={r['name']:r for r in summary if r['stage']=='matrix'}
def compare(a,b):
    ra=read(root/'matrix'/a/'coupled_history.csv'); rb=read(root/'matrix'/b/'coupled_history.csv')
    ta={r['time_s']:r for r in ra if r['active']}; tb={r['time_s']:r for r in rb if r['active']}
    times=sorted(ta.keys() & tb.keys()); result=dict(coarse=a,fine=b,gate=False)
    if not times: return dict(result,reason='no common accepted active state')
    t=times[-1]; keys=('h_m','h_mass_m','evaporation_kg_m2_s','top_solid_kg_m2','internal_solid_kg_m2','epsilon_min')
    differences={k:abs(ta[t][k]-tb[t][k])/max(abs(tb[t][k]),1e-30) for k in keys}
    gate=t==2500 and max(differences[k] for k in keys[:3])<=.02
    return dict(result,time_s=t,relative_differences=differences,gate=gate,reason='final gate' if t==2500 else 'accepted prefix only')
comparisons=[]
if by:
    for n in (40,80,160,320):
        reference={r['time_s']:r for r in read(root/'matrix'/f'pre_n{n}'/'coupled_history.csv')}
        for h in (1,2,5,10):
            for d in (2,5,10):
                name=f'new_n{n}_h{h}_d{d}'; rows=read(root/'matrix'/name/'coupled_history.csv')
                prefix=[r for r in rows if not r['active'] and r['time_s'] in reference]
                delta=max(abs(r[k]-reference[r['time_s']][k]) for r in prefix for k in
                    ('soil_salt_kg','soil_water_kg','top_solid_kg_m2','top_soil_mass_fraction','top_soil_sw','vapour_kg'))
                assert delta==0,(name,'precrust mismatch',delta)
                by[name]['precrust_prefix_max_absolute_difference']=delta
                by[name]['precrust_prefix_common_end_s']=prefix[-1]['time_s']
                by[name]['precrust_prefix_common_rows']=len(prefix)
    for h in (1,2,5,10):
        for d in (2,5,10):
            for a,b in ((40,80),(80,160),(160,320)):
                comparisons.append(compare(f'new_n{a}_h{h}_d{d}',f'new_n{b}_h{h}_d{d}'))
    for n in (160,320):
        c=compare(f'new_n{n}_h1_d5',f'dt_half_n{n}')
        c['kind']='timestep'; c['gate']=c.get('time_s')==2500 and max(c.get('relative_differences',{'missing':1}).values())<=.01
        comparisons.append(c)
    # Old mode must reproduce the previous run, including refused activation.
    regressions=[]
    for name,reference in [('old_n80','nominal'),('old_n160','soil160')]:
        current=read(root/'matrix'/name/'coupled_history.csv')
        frozen=read(root.parent/'p2_h_ale_coupling_20261107/stage_c'/reference/'coupled_history.csv')
        assert len(current)==len(frozen)
        delta=max(abs(a[k]-b[k]) for a,b in zip(current,frozen) for k in a)
        assert delta<1e-12,(name,delta)
        regressions.append(dict(case=name,maximum_absolute_column_difference=delta,pass_gate=True))
    for name in ('rollback_activation','rollback_coupled'):
        assert by[name]['rejected_trials']>=1 and any('injected' in r for r in by[name]['refusal_reasons'])
else: regressions=[]
precondition320=None
if by:
    names=('old_n320','pre_n320','new_n320_h1_d2')
    hashes={name:hashlib.sha256((root/'matrix'/name/'coupled_history.csv').read_bytes()).hexdigest() for name in names}
    precondition320=dict(history_sha256=hashes,identical_histories=len(set(hashes.values()))==1,
        accepted_times_s={name:by[name]['accepted_time_s'] for name in names},
        activation_times_s={name:by[name]['activation_s'] for name in names})
    dt_hash=hashlib.sha256((root/'matrix/dt_half_n320/coupled_history.csv').read_bytes()).hexdigest()
    precondition320['dt_half_history_sha256']=dt_hash
    precondition320['dt_half_identical_history']=dt_hash==hashes['old_n320']
refusal_check=None
if (root/'refusal/depth_refused/exit_status.txt').exists():
    a=read(root/'refusal/depth_refused/coupled_history.csv');b=read(root/'refusal/precrust_reference/coupled_history.csv')
    assert len(a)==len(b) and a[-1]['time_s']==400
    assert all(x==y for x,y in zip(a,b)), 'refused donation mutated accepted transport state'
    refusal_check=dict(equal_entire_accepted_history=True,rows=len(a),final_time_s=400)
depth_sensitivity=[]; precrust_comparison=[]; threshold_sensitivity=[]; age_diagnostic=[]
if by:
    for n in (40,80,160,320):
        for h in (1,2,5,10):
            for a,b in ((2,5),(5,10)):
                c=compare(f'new_n{n}_h{h}_d{a}',f'new_n{n}_h{h}_d{b}');c['kind']='donor_depth';depth_sensitivity.append(c)
    for n in (40,80,160,320):
        for a,b in ((1,2),(2,5),(5,10)):
            c=compare(f'new_n{n}_h{a}_d5',f'new_n{n}_h{b}_d5');c.pop('gate',None)
            c['kind']='activation_threshold_descriptive_only';threshold_sensitivity.append(c)
        r=by[f'pre_n{n}']; history=read(root/'matrix'/r['name']/'coupled_history.csv')
        selected=[x for x in history if x['time_s'] in (0,100,500,1000,2500)]
        precrust_comparison.append(dict(n=n,first_precipitation_s=r['precipitation_s'],thresholds_s=r['thresholds_s'],
            samples=[{k:x[k] for k in ('time_s','top_solid_kg_m2','h_mass_m','top_soil_mass_fraction','top_soil_sw')} for x in selected]))
    # Secondary phase-lag diagnostic only: NOT the common-physical-time gate.
    for h in (1,2,5,10):
        trajectories={}
        for n in (40,80,160,320):
            rows=read(root/'matrix'/f'new_n{n}_h{h}_d5'/'coupled_history.csv')
            trajectories[n]={r['time_s']-r['activation_s']:r for r in rows if r['active']}
        available=[set(v) for v in trajectories.values() if v]
        ages=set.intersection(*available) if available else set()
        if not ages: continue
        age=max(ages)
        for a,b in ((40,80),(80,160),(160,320)):
            if age not in trajectories[a] or age not in trajectories[b]:continue
            ra=trajectories[a][age];rb=trajectories[b][age]
            age_diagnostic.append(dict(h_um=h,coarse=a,fine=b,age_s=age,
                coarse_physical_time_s=ra['time_s'],fine_physical_time_s=rb['time_s'],
                height_relative=abs(ra['h_m']-rb['h_m'])/rb['h_m'],
                solid_relative=abs(ra['h_mass_m']-rb['h_mass_m'])/rb['h_mass_m']))
result=dict(cases=summary,controlled=controlled,comparisons=comparisons,depth_sensitivity=depth_sensitivity,
    threshold_sensitivity=threshold_sensitivity,
    activation_age_diagnostic=age_diagnostic,
    precrust_comparison=precrust_comparison,old_mode_regressions=regressions,refusal_check=refusal_check,
    preactivation_320_control=precondition320,
    unchanged_physics_sources=unchanged,integrity_errors=errors,physical_validation=False)
(out/'verification.json').write_text(json.dumps(result,indent=2))
with (out/'convergence.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=['coarse','fine','time_s','gate','reason','h_relative','solid_relative','evaporation_relative']);w.writeheader()
    for c in comparisons:
        d=c.get('relative_differences',{}); row={k:c.get(k) for k in ['coarse','fine','time_s','gate','reason']}
        row.update(h_relative=d.get('h_m'),solid_relative=d.get('h_mass_m'),evaporation_relative=d.get('evaporation_kg_m2_s'));w.writerow(row)
print(json.dumps(dict(cases=len(summary),controlled=len(controlled),comparisons=len(comparisons),integrity_errors=errors),indent=2))
for name,items in [('case_summary',summary),('controlled_allocation',controlled)]:
    fields=['name','n','h_um','depth_mm','mode','status','exit_code','accepted_time_s','precipitation_s','first_attempt_s','activation_s',
        'salt_residual_relative','water_residual_relative','required_kg','available_kg','sufficient','donor_count','used_depth_m',
        'pressure_change_max_pa','saturation_change_max','concentration_change_max','transfer_salt_max_relative','transfer_water_max_relative']
    with (out/(name+'.csv')).open('w',newline='') as f:
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader()
        for item in items:
            data=dict(item);data.update(item.get('allocation',{}));writer.writerow({k:data.get(k) for k in fields})
