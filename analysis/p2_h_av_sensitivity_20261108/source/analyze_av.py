"""Conservative post-processing for P2-H reactive-area sensitivity runs."""
import csv,json,math,pathlib,sys
from collections import defaultdict

root=pathlib.Path.cwd()
output=root/'analysis'
results=output/'results'
for path in (results,results/'time_series',results/'reaction_profiles',results/'porosity_profiles'):
    path.mkdir(parents=True,exist_ok=True)

def read(path):
    with path.open(newline='') as handle:
        return [{key:float(value) if key not in ('status',) else value for key,value in row.items()}
                for row in csv.DictReader(handle)]

def write(path,rows,fields=None):
    rows=list(rows)
    if not rows:
        return
    fields=fields or list(rows[0])
    with path.open('w',newline='') as handle:
        writer=csv.DictWriter(handle,fields,extrasaction='ignore')
        writer.writeheader(); writer.writerows(rows)

def inverse_mass(cells,target):
    if target<=0:
        return cells[0]['z_bottom_m']
    accumulated=0
    for cell in cells:
        density=max(cell['reaction_kg_m3_s'],0.)
        mass=density*(cell['z_top_m']-cell['z_bottom_m'])
        if accumulated+mass>=target and density>0:
            return cell['z_bottom_m']+(target-accumulated)/density
        accumulated+=mass
    return cells[-1]['z_top_m']

def reaction_metric(cells):
    cells=sorted(cells,key=lambda row:row['z_bottom_m'])
    masses=[max(row['reaction_kg_m3_s'],0.)*(row['z_top_m']-row['z_bottom_m']) for row in cells]
    total=sum(masses)
    result=dict(time_s=cells[0]['time_s'],status=cells[0]['status'],h_m=cells[0]['h_m'])
    if total<=1e-18:
        result.update(positive_rate_kg_m2_s=total,l95_bottom_m=math.nan,l95_shortest_m=math.nan,
                      l95_bottom_over_h=math.nan,l95_shortest_over_h=math.nan,reaction_centroid_m=math.nan,
                      maximum_position_m=math.nan,maximum_reaction_kg_m3_s=0,reacting_volume_m=0)
        return result
    target=.95*total
    bottom=inverse_mass(cells,target)-cells[0]['z_bottom_m']
    edges=[0.]
    for mass in masses:
        edges.append(edges[-1]+mass)
    starts={0.,.05*total}
    for edge in edges:
        if 0<=edge<=.05*total: starts.add(edge)
        shifted=edge-target
        if 0<=shifted<=.05*total: starts.add(shifted)
    shortest=min(inverse_mass(cells,start+target)-inverse_mass(cells,start) for start in starts)
    centroid=sum(.5*(row['z_bottom_m']+row['z_top_m'])*mass for row,mass in zip(cells,masses))/total
    maximum=max(cells,key=lambda row:row['reaction_kg_m3_s'])
    peak=maximum['reaction_kg_m3_s']
    reacting=sum(row['z_top_m']-row['z_bottom_m'] for row in cells if row['reaction_kg_m3_s']>max(peak*1e-12,1e-18))
    h=result['h_m']
    result.update(positive_rate_kg_m2_s=total,l95_bottom_m=bottom,l95_shortest_m=shortest,
                  l95_bottom_over_h=bottom/h,l95_shortest_over_h=shortest/h,
                  reaction_centroid_m=centroid,
                  maximum_position_m=.5*(maximum['z_bottom_m']+maximum['z_top_m']),
                  maximum_reaction_kg_m3_s=peak,reacting_volume_m=reacting)
    return result

specs=[]
for stage in ('primary_exact','dt_half_exact','ale_refine_exact','ale256_exact'):
    casefile=root/stage/'cases.json'
    if casefile.exists():
        for item in json.loads(casefile.read_text()):
            item=dict(item,stage=stage,path=root/stage/item['name'])
            specs.append(item)

summary=[]; all_metrics={}; histories={}; failures=[]
for spec in specs:
    path=spec['path']
    history_file=path/'coupled_history.csv'
    if not history_file.exists():
        failures.append(dict(case=spec['name'],stage=spec['stage'],reason='missing coupled_history.csv'))
        continue
    history=read(history_file); histories[spec['name']]=history
    accepted=defaultdict(list); rejected=defaultdict(list)
    diagnostics=path/'ale_diagnostics.csv'
    if diagnostics.exists():
        for row in read(diagnostics):
            (accepted if row['status']=='accepted' else rejected)[row['time_s']].append(row)
    metrics=[reaction_metric(accepted[time]) for time in sorted(accepted)]
    all_metrics[spec['name']]=metrics
    write(results/'time_series'/f"{spec['name']}.csv",
          [row for index,row in enumerate(history) if index==0 or index==len(history)-1 or row['active']==1 and
           (abs(row['time_s']-row['activation_s'])<1e-10 or abs(row['time_s']/10-round(row['time_s']/10))<1e-10)])
    write(results/'time_series'/f"{spec['name']}_reaction_metrics.csv",metrics)
    selected=set()
    active=[row for row in history if row['active']==1]
    if active:
        activation=active[0]['time_s']; available=sorted(accepted)
        for target in (activation,activation+1,activation+5,activation+20,activation+100,activation+500,history[-1]['time_s']):
            if available and target<=available[-1]+1e-10:
                selected.add(min(available,key=lambda time:abs(time-target)))
    profile_rows=[row for time in sorted(selected) for row in accepted[time]]
    write(results/'reaction_profiles'/f"{spec['name']}.csv",profile_rows)
    write(results/'porosity_profiles'/f"{spec['name']}.csv",profile_rows)
    first=history[0]; last=history[-1]
    activation=next((row for row in history if row['active']==1),None)
    metric=min(metrics,key=lambda row:abs(row['time_s']-last['time_s'])) if metrics else {}
    exit_code=int((path/'exit_status.txt').read_text()) if (path/'exit_status.txt').exists() else -999
    completed=exit_code==0 and last['time_s']>=2499.999 if spec['stage']!='smoke' else exit_code==0
    first_rejected_time=min(rejected) if rejected else math.nan
    rejected_min=min((row for rows in rejected.values() for row in rows),key=lambda row:row['epsilon'],default={})
    salt_scale=first['soil_salt_kg']+max(first['top_solid_kg_m2'],1e-30)*.001
    water_scale=first['soil_water_kg']
    threshold_times={}
    for threshold in (.20,.10,.02):
        crossing=next((row['time_s'] for row in history if row['active']==1 and row['epsilon_min']<=threshold),math.nan)
        threshold_times[f't_epsilon_{str(threshold).replace(".","p")}_s']=crossing
    summary.append(dict(case=spec['name'],stage=spec['stage'],av_multiplier=spec['av'],reaction_coefficient_m_s=1e-3,
        soil_cells=spec['n'],ale_cells=spec['ale'],dt_s=spec['dt'],activation_s=activation['time_s'] if activation else math.nan,
        initial_ale_concentration=next(iter(accepted.values()))[0]['c'] if accepted else math.nan,
        actual_end_s=last['time_s'],final_height_m=last['h_m'],final_top_growth_velocity_m_s=last['top_growth_velocity_m_s'],
        cumulative_positive_internal_kg_m2=last['internal_positive_kg_m2'],
        cumulative_internal_dissolution_kg_m2=last['internal_dissolution_kg_m2'],
        net_internal_solid_kg_m2=last['internal_solid_kg_m2'],cumulative_new_top_solid_kg_m2=last['top_new_solid_kg_m2'],
        internal_fraction_net=last['internal_fraction_net'],internal_fraction_positive=last['internal_fraction_positive'],
        l95_bottom_m=metric.get('l95_bottom_m',math.nan),l95_shortest_m=metric.get('l95_shortest_m',math.nan),
        l95_bottom_over_h=metric.get('l95_bottom_over_h',math.nan),l95_shortest_over_h=metric.get('l95_shortest_over_h',math.nan),
        epsilon_min=last['epsilon_min'],salt_relative_residual=max(abs(row['salt_residual_kg']) for row in history)/salt_scale,
        water_relative_residual=max(abs(row['water_residual_kg']) for row in history)/water_scale,
        quadrature_net_relative=max(abs(row['internal_quadrature_residual_kg_m2']) for row in history)/max(abs(last['internal_solid_kg_m2']),1e-30),
        first_rejected_time_s=first_rejected_time,rejected_epsilon_min=rejected_min.get('epsilon',math.nan),
        rejected_cell=int(rejected_min.get('i',-1)),rejected_cell_z_bottom_m=rejected_min.get('z_bottom_m',math.nan),
        rejected_c=rejected_min.get('c',math.nan),rejected_av_m_inv=rejected_min.get('av_m_inv',math.nan),
        rejected_reaction_kg_m3_s=rejected_min.get('reaction_kg_m3_s',math.nan),rejected_deff_m2_s=rejected_min.get('deff_m2_s',math.nan),
        rejected_bottom_liquid_kg_m2_s=rejected_min.get('bottom_liquid_kg_m2_s',math.nan),
        rejected_bottom_salt_kg_m2_s=rejected_min.get('bottom_salt_kg_m2_s',math.nan),
        exit_code=exit_code,status='COMPLETED' if completed else 'STOPPED',**threshold_times))
    if rejected:
        first_time=min(rejected)
        write(results/'reaction_profiles'/f"{spec['name']}_first_rejected.csv",rejected[first_time])
    if not completed and spec['stage']!='smoke':
        failures.append(dict(case=spec['name'],stage=spec['stage'],reason='early stop',last_accepted_s=last['time_s'],
                             first_rejected_s=first_rejected_time,rejected_epsilon_min=rejected_min.get('epsilon',math.nan)))

write(results/'summary.csv',summary)
write(results/'failures.csv',failures)

# Strict preactivation identity: physical columns only; alpha is intentionally different.
identity=[]
for n in (80,160):
    names=[f'n{n}_av{av:g}' for av in (.1,1.,10.)]
    if not all(name in histories for name in names): continue
    ref=histories[names[0]]
    stop=min(next(row['time_s'] for row in histories[name] if row['active']==1) for name in names)
    keys=('time_s','h_m','h_mass_m','top_solid_kg_m2','soil_salt_kg','soil_water_kg','top_soil_mass_fraction','top_soil_sw')
    for name in names[1:]:
        other=histories[name]
        pairs=[(a,b) for a,b in zip(ref,other) if a['time_s']<=stop and b['time_s']<=stop]
        identity.append(dict(soil_cells=n,reference=names[0],case=name,through_s=stop,
            row_count_reference=sum(row['time_s']<=stop for row in ref),row_count_case=sum(row['time_s']<=stop for row in other),
            maximum_absolute_difference=max(abs(a[key]-b[key]) for a,b in pairs for key in keys)))
write(results/'preactivation_identity.csv',identity)

activation_identity=[]
for n in (80,160):
    names=[f'n{n}_av{av:g}' for av in (.1,1.,10.)]
    if not all(name in histories for name in names): continue
    def initial_state(name):
        case=root/'primary_exact'/name
        diagnostic=[row for row in read(case/'ale_diagnostics.csv') if row['status']=='accepted']
        activation=min(row['time_s'] for row in diagnostic)
        ale=[row for row in diagnostic if row['time_s']==activation]
        soil=[row for row in read(case/'soil_fields.csv') if row['time_s']==activation]
        donor=read(case/'donor_profiles.csv')
        return activation,ale,soil,donor
    reference=initial_state(names[0])
    for name in names[1:]:
        current=initial_state(name)
        ale_keys=('time_s','z_bottom_m','z_top_m','epsilon','c','h_m','speed_m_s','bottom_c',
                  'pore_salt_kg_m2','pore_water_kg_m2','solid_kg_m2')
        soil_keys=('time_s','y_m','mass_fraction','saturation','density_kg_m3')
        donor_keys=tuple(reference[3][0]) if reference[3] else ()
        activation_identity.append(dict(soil_cells=n,reference=names[0],case=name,
            activation_difference_s=abs(reference[0]-current[0]),
            initial_ale_max_abs=max(abs(a[key]-b[key]) for a,b in zip(reference[1],current[1]) for key in ale_keys),
            post_transfer_soil_max_abs=max(abs(a[key]-b[key]) for a,b in zip(reference[2],current[2]) for key in soil_keys),
            donor_profile_max_abs=max((abs(a[key]-b[key]) for a,b in zip(reference[3],current[3]) for key in donor_keys),default=0),
            ale_rows_reference=len(reference[1]),ale_rows_case=len(current[1]),soil_rows_reference=len(reference[2]),soil_rows_case=len(current[2]),
            donor_rows_reference=len(reference[3]),donor_rows_case=len(current[3])))
write(results/'activation_identity.csv',activation_identity)

def interpolate(rows,key,target):
    rows=sorted(rows,key=lambda row:row[key])
    if target<rows[0][key] or target>rows[-1][key]: return None
    for left,right in zip(rows,rows[1:]):
        if left[key]<=target<=right[key]:
            fraction=0 if right[key]==left[key] else (target-left[key])/(right[key]-left[key])
            out={name:left[name]+fraction*(right[name]-left[name]) for name in left if isinstance(left[name],float)}
            out[key]=target
            return out
    return rows[-1]

comparison_keys=('h_m','top_growth_velocity_m_s','internal_positive_kg_m2','internal_dissolution_kg_m2',
                 'internal_solid_kg_m2','top_new_solid_kg_m2','internal_fraction_net','epsilon_min',
                 'top_soil_mass_fraction','interface_salt_kg_m2')
physical=[]; elapsed=[]; by_height=[]
for n in (80,160):
    available=[f'n{n}_av{av:g}' for av in (.1,1.,10.) if f'n{n}_av{av:g}' in histories]
    if len(available)!=3: continue
    common_start=max(next(row['activation_s'] for row in histories[name] if row['active']) for name in available)
    common_end=min(histories[name][-1]['time_s'] for name in available)
    for time in (676.,680.,695.,725.,740.,1000.,1500.,2000.,2500.):
        if common_start<=time<=common_end:
            for name in available:
                row=interpolate(histories[name],'time_s',time)
                physical.append(dict(case=name,soil_cells=n,av_multiplier=float(name.rsplit('av',1)[1]),time_s=time,
                                     **{key:row[key] for key in comparison_keys}))
    for duration in (1.,5.,20.,100.,500.,1000.):
        if common_start+duration<=common_end:
            for name in available:
                activation=next(row['activation_s'] for row in histories[name] if row['active'])
                row=interpolate(histories[name],'time_s',activation+duration)
                elapsed.append(dict(case=name,soil_cells=n,av_multiplier=float(name.rsplit('av',1)[1]),elapsed_after_activation_s=duration,
                                    physical_time_s=activation+duration,**{key:row[key] for key in comparison_keys}))
    active={name:[row for row in histories[name] if row['active']] for name in available}
    minimum_height=max(rows[0]['h_m'] for rows in active.values()); maximum_height=min(rows[-1]['h_m'] for rows in active.values())
    for height_um in (3.52,3.56,3.60,3.64,4.,6.,8.,10.):
        target=height_um*1e-6
        if minimum_height<=target<=maximum_height:
            for name in available:
                row=interpolate(active[name],'h_m',target)
                by_height.append(dict(case=name,soil_cells=n,av_multiplier=float(name.rsplit('av',1)[1]),height_um=height_um,
                                      physical_time_s=row['time_s'],**{key:row[key] for key in comparison_keys if key!='h_m'}))
write(results/'comparison_physical_time.csv',physical)
write(results/'comparison_elapsed_time.csv',elapsed)
write(results/'comparison_height.csv',by_height)

# A3 is the only case whose 64-cell reaction zone failed the exploratory ALE
# resolution gate.  Compare the 64/128/256 trajectories only where every
# trajectory has an accepted state; never compare their unequal stop endpoints.
resolution=[]
resolution_names=('n80_av10','n80_av10_ale128','n80_av10_ale256')
if all(name in histories for name in resolution_names):
    resolution_fields=('h_m','internal_positive_kg_m2','internal_dissolution_kg_m2',
                       'internal_solid_kg_m2','top_new_solid_kg_m2','internal_fraction_net','epsilon_min')
    for time in (700.,725.,730.):
        if all(histories[name][0]['time_s']<=time<=histories[name][-1]['time_s'] for name in resolution_names):
            for name in resolution_names:
                row=interpolate(histories[name],'time_s',time)
                resolution.append(dict(comparison='physical_time',target=time,case=name,
                    ale_cells=64 if name=='n80_av10' else int(name.rsplit('ale',1)[1]),
                    physical_time_s=time,**{key:row[key] for key in resolution_fields}))
    active_resolution={name:[row for row in histories[name] if row['active']] for name in resolution_names}
    for height_um in (3.52,3.56,3.60):
        target=height_um*1e-6
        if all(rows[0]['h_m']<=target<=rows[-1]['h_m'] for rows in active_resolution.values()):
            for name,rows in active_resolution.items():
                row=interpolate(rows,'h_m',target)
                resolution.append(dict(comparison='height_um',target=height_um,case=name,
                    ale_cells=64 if name=='n80_av10' else int(name.rsplit('ale',1)[1]),
                    physical_time_s=row['time_s'],**{key:row[key] for key in resolution_fields}))
write(results/'ale_resolution_common_state.csv',resolution)

verification=dict(
    case_count=len(summary),
    completed_count=sum(row['status']=='COMPLETED' for row in summary),
    preactivation_identity=identity,
    activation_identity=activation_identity,
    all_zero_initial_solid=all(histories[row['case']][0]['top_solid_kg_m2']==0 and histories[row['case']][0]['h_mass_m']==0 for row in summary),
    maximum_salt_relative_residual=max((row['salt_relative_residual'] for row in summary),default=math.nan),
    maximum_water_relative_residual=max((row['water_relative_residual'] for row in summary),default=math.nan),
    stopped_cases=[row['case'] for row in summary if row['status']!='COMPLETED'])
(output/'verification.json').write_text(json.dumps(verification,indent=2,allow_nan=True))
print(json.dumps(verification,indent=2,allow_nan=True))
