import pathlib,csv,json,math
root=pathlib.Path.cwd()
out=root/'results'
summary=[]
ceq=.1*.05844/(.9*.01801528+.1*.05844)
for path in sorted((root/'runs').iterdir()):
    with (path/'surface_precipitation.csv').open() as f:
        rows=[{k:float(v) for k,v in r.items()} for r in csv.DictReader(f)]
    first,last=rows[0],rows[-1]
    onset=next((r['time_s'] for r in rows if r['precipitated_step_kg_m2']>1e-15),None)
    switch=next((r['time_s'] for r in rows if r['prat_active']),None)
    compare=[r for r in rows if r['prat_demand_kg_m2']>0]
    actual=sum(r['precipitated_step_kg_m2'] for r in compare)
    demand=sum(r['prat_demand_kg_m2'] for r in compare)
    sr=max(abs(r['salt_residual_kg']) for r in rows)/first['soil_salt_kg']
    wr=max(abs(r['water_residual_kg']) for r in rows)/first['soil_water_kg']
    item=dict(case=path.name,time_s=last['time_s'],onset_s=onset,continuum_s=switch,
              solid_kg_m2=last['solid_kg_m2'],h_mass_m=last['h_mass_m'],h_prat_m=last['h_prat_m'],
              actual_comparison_kg_m2=actual,prat_demand_kg_m2=demand,
              demand_over_actual=demand/actual if actual else None,
              salt_residual_relative=sr,water_residual_relative=wr,
              max_top_mass_fraction=max(r['top_mass_fraction'] for r in rows))
    assert sr<1e-8 and wr<1e-8,item
    if path.name!='no_precip_control':
        assert item['max_top_mass_fraction']<=ceq+1e-12,item
        assert onset is not None and switch is not None,item
    # Compact exact records: include onset/switch and final state plus regular samples.
    chosen=[r for i,r in enumerate(rows) if i%20==0 or i==len(rows)-1 or r['time_s'] in (onset,switch)]
    with (out/(path.name+'.csv')).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(chosen)
    summary.append(item)
with (out/'summary.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=summary[0]);w.writeheader();w.writerows(summary)
(out/'summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2))
