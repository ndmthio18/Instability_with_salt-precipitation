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
        assert all(first[k]==0 for k in ('h_mass_m','internal_solid_kg_m2','crust_water_kg','crust_salt_kg'))
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


# Independently validate the actual runtime BOX boundary faces, not a reconstructed
# nominal grid. Geometric normals are independent of the production predicate.
for item in controlled+summary:
    path=root/item['stage']/item['name']; faces=read(path/'boundary_faces.csv')
    cfg=read(path/'activation_config.csv')[0]
    assert cfg['threshold_m']==item['h_um']*1e-6 and cfg['birth_porosity']==.43
    assert cfg['crust_cells']==64 and cfg['initial_height_m']==0 and cfg['initial_solid_kg_m2']==0
    assert faces and len(set(f['id'] for f in faces))==len(faces)
    for f in faces:
        expected=abs(f['ny']-1)<1e-12 and abs(f['y_m']-.2)<1e-12
        assert bool(f['selected_top'])==expected,(item['name'],f)
    actual=sum(f['area_m2'] for f in faces if f['selected_top'])
    legacy_side=sum(f['area_m2'] for f in faces if f['legacy_top'] and abs(f['ny'])<.5)
    assert abs(actual-.001)<1e-15
    assert (legacy_side>0)==(item['n']==320)
    item['boundary']=dict(top_area_m2=actual,legacy_side_area_m2=legacy_side,selected_side_area_m2=0)
    item['activation_config']=cfg
    if item['name'].startswith('crystal_') or item['stage']=='smoke':
        args=json.loads((path/'command.json').read_text());assert '-P2H.ActivationHeight' not in args
    if item.get('allocation',{}).get('sufficient'):
        a=item['allocation']; assert a['h_m']>=cfg['threshold_m']-1e-19
        mixed=a['salt_transfer_kg']/(a['salt_transfer_kg']+a['water_transfer_kg'])
        ceq=.1*.05844/(.1*.05844+.9*.01801528)
        ps=[r for r in read(path/'donor_profiles.csv') if r['attempt']==a['attempt'] and r['extracted_kg']>0]
        a['initial_mixed_concentration']=mixed
        a['donor_concentration_min']=min(r['c_before'] for r in ps)
        a['donor_concentration_max']=max(r['c_before'] for r in ps)
        a['initial_dissolved_excess_kg_m2']=max(0,a['salt_transfer_kg']-ceq/(1-ceq)*a['water_transfer_kg'])/.001
        fields=read(path/'crust_fields.csv')
        item['stored_crust_times_s']=sorted(set(f['time_s'] for f in fields))
        born=[f for f in fields if f['time_s']==item['activation_s']]
        if born:
            assert len(born)==64 and all(abs(f['epsilon']-.43)<1e-14 for f in born)
            assert abs(max(f['z_top_m'] for f in born)-a['h_m'])<1e-14
            mixed=a['salt_transfer_kg']/(a['salt_transfer_kg']+a['water_transfer_kg'])
            assert all(abs(f['concentration']-mixed)<1e-14 for f in born)
            history=read(path/'coupled_history.csv')
            birth=next(r for r in history if r['time_s']==item['activation_s'])
            assert abs(birth['h_m']-birth['top_solid_kg_m2']/(2165*.57))<1e-14
            assert birth['internal_solid_kg_m2']==0
            last=history[-1];active=[r for r in history if r['active']]
            item['post_activation_budget']=dict(duration_s=last['time_s']-birth['time_s'],
                interface_salt_kg_m2=last['interface_salt_kg_m2'],
                interface_water_kg_m2=last['interface_liquid_kg_m2']-last['interface_salt_kg_m2'],
                dissolved_salt_change_kg_m2=(last['crust_salt_kg']-birth['crust_salt_kg'])/.001,
                top_solid_change_kg_m2=last['top_solid_kg_m2']-birth['top_solid_kg_m2'],
                internal_solid_kg_m2=last['internal_solid_kg_m2'],
                max_top_soil_concentration=max(r['top_soil_mass_fraction'] for r in active))
    if item.get('precrust'):
        item['threshold_3p5_s']=next((r['time_s'] for r in read(path/'coupled_history.csv') if r['h_mass_m']>=3.5e-6),None)
by={r['name']:r for r in summary if r['stage']=='matrix'}
comparisons=[]; regressions=[]
def compare(a,b,kind):
    ra=read(root/'matrix'/a/'coupled_history.csv');rb=read(root/'matrix'/b/'coupled_history.csv')
    ta={r['time_s']:r for r in ra if r['active']};tb={r['time_s']:r for r in rb if r['active']}
    common=sorted(ta.keys() & tb.keys()); result=dict(coarse=a,fine=b,kind=kind,gate=False)
    if not common:return dict(result,reason='no exactly matched accepted active timestamp; final state unavailable')
    t=common[-1];keys=('h_m','h_mass_m','evaporation_kg_m2_s','top_solid_kg_m2','internal_solid_kg_m2','epsilon_min')
    delta={k:abs(ta[t][k]-tb[t][k])/max(abs(tb[t][k]),1e-30) for k in keys}
    return dict(result,time_s=t,relative_differences=delta,gate=t==2500 and max(delta[k] for k in keys[:3])<=(.01 if kind=='timestep' else .02))
if by:
    for n in (40,80,160,320):
        c=by[f'crystal_n{n}'];ref={r['time_s']:r for r in read(root/'matrix'/f'pre_n{n}'/'coupled_history.csv')}
        rows=read(root/'matrix'/f'crystal_n{n}'/'coupled_history.csv')
        common=[r for r in rows if not r['active'] and r['time_s'] in ref]
        diff=max(abs(r[k]-ref[r['time_s']][k]) for r in common for k in ('soil_salt_kg','soil_water_kg','top_solid_kg_m2','top_soil_sw','top_soil_mass_fraction'))
        assert diff==0;c['precrust_prefix_max_difference']=diff;c['precrust_common_end_s']=common[-1]['time_s']
    for a,b in ((40,80),(80,160),(160,320)): comparisons.append(compare(f'crystal_n{a}',f'crystal_n{b}','grid'))
    for n in (160,320):comparisons.append(compare(f'crystal_n{n}',f'dt_half_n{n}','timestep'))
    for n in (80,160):
        a=root/'matrix'/f'old_n{n}'/'coupled_history.csv';b=root.parent/'p2_h_adaptive_donor_20261108/matrix'/f'old_n{n}'/'coupled_history.csv'
        equal=a.read_bytes()==b.read_bytes();assert equal
        regressions.append(dict(case=f'old_n{n}',identical_history=equal,sha256=hashlib.sha256(a.read_bytes()).hexdigest()))
result=dict(cases=summary,controlled=controlled,comparisons=comparisons,regressions=regressions,unchanged_physics_sources=unchanged,integrity_errors=errors,physical_validation=False)
(out/'verification.json').write_text(json.dumps(result,indent=2))
for name,items in [('case_summary',summary),('controlled_allocation',controlled),('convergence',comparisons)]:
    flat=[]
    for c in items:
        row={k:v for k,v in c.items() if not isinstance(v,(dict,list))}
        for k,v in c.items():
            if isinstance(v,dict):row.update({k+'.'+a:b for a,b in v.items() if not isinstance(b,(list,dict))})
        flat.append(row)
    if flat:
        with (out/(name+'.csv')).open('w',newline='') as f:
            w=csv.DictWriter(f,fieldnames=sorted(set().union(*(set(c) for c in flat))));w.writeheader();w.writerows(flat)
print(json.dumps(dict(cases=len(summary),controlled=len(controlled),integrity_errors=errors,comparisons=comparisons)),flush=True)
if len(sys.argv)>1 and sys.argv[1]=='smoke':
    assert len(controlled)==4 and all(c['allocation']['sufficient'] for c in controlled)
    assert len(summary)==1 and summary[0]['accepted_time_s']==120 and summary[0]['exit_code']==0
assert not errors
