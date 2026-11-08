"""Independent fields, ledgers, regressions and predeclared scientific gates."""
import hashlib,json,os,pathlib
import numpy as np
import pandas as pd
assert os.environ.get('SLURM_JOB_ID')
root=pathlib.Path.cwd();out=root/'audit';out.mkdir(exist_ok=True)
def read(p):return pd.read_csv(p)
def save(rows,name):pd.DataFrame(rows).to_csv(out/name,index=False)
checks=[];convergence=[];activation=[];budget=[];failures=[];independent_series=[]
for stage in ('diagnostic','ale_diagnostic','correction/fixed_smoke','correction/fixed_tests','correction/fixed_pre','correction/fixed_coupled'):
    for d in sorted((root/stage).glob('*')):
        if not d.is_dir() or not (d/'exit_status.txt').exists():continue
        h=read(d/'coupled_history.csv');a,b=h.iloc[0],h.iloc[-1];active=h[h.active==1]
        assert a.h_m==0 and a.top_solid_kg_m2==0 and a.internal_solid_kg_m2==0
        # An independent reconstruction from DOF volumes, rho, Sw and C.
        f=read(d/'grid_profiles.csv');f['salt']=.41*f.box_volume_m3*f.sw*f.rho*f.c
        f['water']=.41*f.box_volume_m3*f.sw*f.rho*(1-f.c)
        m=f.groupby('time_s')[['salt','water']].sum();joined=h.set_index('time_s').join(m,how='inner')
        independent=joined.copy()
        es=(joined.soil_salt_kg-joined.salt).abs().max()/a.soil_salt_kg
        ew=(joined.soil_water_kg-joined.water).abs().max()/a.soil_water_kg
        assert max(es,ew)<1e-12
        q=read(d/'grid_budgets.csv');q=q[q.active==0]
        implicit='fixed' in stage and 'projection' not in d.name and not d.name.startswith('undersaturated_projection')
        if implicit and len(q):
            assert q.projection_salt_kg.abs().max()<1e-20,'implicit mode projected salt'
            top=h[h.active==0].top_soil_mass_fraction
            assert top.max()<=.2649408697676898+5.1e-13
            physical=q[q.depth_m==q.depth_m.max()]
            precip=(physical.dt_s*physical.boundary_precip_kg_m2_s).sum()
            lastpre=h[h.active==0].iloc[-1]
            assert abs(precip-lastpre.top_solid_kg_m2)<1e-12
        item=dict(case=d.name,stage=stage,accepted_time_s=b.time_s,exit_code=int((d/'exit_status.txt').read_text()),
            global_salt_relative=h.salt_residual_kg.abs().max()/a.soil_salt_kg,
            global_water_relative=h.water_residual_kg.abs().max()/a.soil_water_kg,
            independent_soil_salt_relative=es,independent_soil_water_relative=ew)
        assert max(item['global_salt_relative'],item['global_water_relative'])<=1e-8
        if len(active):
            first=active.iloc[0];ad=read(d/'ale_diagnostics.csv');ad=ad[ad.status=='accepted'].copy()
            rho=(first.crust_salt_kg+first.crust_water_kg)/(.001*.43*first.h_m)
            ad['dz']=ad.z_top_m-ad.z_bottom_m
            ad['salt']=rho*ad.epsilon*ad.c*ad.dz*.001
            ad['water']=rho*ad.epsilon*(1-ad.c)*ad.dz*.001
            ad['solid']=2165*(1-ad.epsilon)*ad.dz*.001
            fields=ad.groupby('time_s')[['salt','water','solid']].sum()
            independent=independent.join(fields.rename(columns={'salt':'ale_salt','water':'ale_water','solid':'ale_solid'}))
            joined=active.set_index('time_s').join(fields,how='inner')
            assert len(joined)==len(active)
            cs=(joined.crust_salt_kg-joined.salt).abs().max()/a.soil_salt_kg
            cw=(joined.crust_water_kg-joined.water).abs().max()/a.soil_water_kg
            ss=(joined.solid-.001*(joined.top_solid_kg_m2+joined.internal_solid_kg_m2)).abs().max()/a.soil_salt_kg
            assert max(cs,cw,ss)<1e-12
            si=(active.soil_salt_kg-first.soil_salt_kg+.001*active.interface_salt_kg_m2).abs().max()/a.soil_salt_kg
            wi=(active.soil_water_kg-first.soil_water_kg+.001*(active.interface_liquid_kg_m2-active.interface_salt_kg_m2)).abs().max()/a.soil_water_kg
            assert max(si,wi)<1e-8
            da=read(d/'donor_attempts.csv');good=da[da.sufficient==1].iloc[-1]
            assert abs(good.salt_transfer_residual_kg)/a.soil_salt_kg<=1e-12
            assert abs(good.water_transfer_residual_kg)/a.soil_water_kg<=1e-12
            item.update(independent_crust_salt_relative=cs,independent_crust_water_relative=cw,
                        independent_solid_relative=ss,interface_salt_relative=si,interface_water_relative=wi)
            activation.append(dict(case=d.name,stage=stage,time_s=first.time_s,h_um=first.h_m*1e6,
                c_initial=good.salt_transfer_kg/(good.salt_transfer_kg+good.water_transfer_kg),donor_depth_mm=good.used_depth_m*1000,
                donor_count=good.donor_count,transfer_salt_relative=abs(good.salt_transfer_residual_kg)/a.soil_salt_kg,
                transfer_water_relative=abs(good.water_transfer_residual_kg)/a.soil_water_kg))
        for t,r in independent.iterrows():
            crust_salt=(r['ale_salt']+r['ale_solid']) if r.active else r.top_solid_kg_m2*.001
            crust_water=r['ale_water'] if r.active else 0
            sr=(r.salt+crust_salt-independent.salt.iloc[0])/independent.salt.iloc[0]
            wr=(r.water+crust_water+r.vapour_kg-independent.water.iloc[0])/independent.water.iloc[0]
            assert max(abs(sr),abs(wr))<1e-8
            independent_series.append(dict(case=d.name,stage=stage,time_s=t,salt_relative=sr,water_relative=wr))
        budget.append(item)
        if b.time_s<2500 and stage not in ('correction/fixed_smoke','correction/fixed_tests') and d.name!='historical160_ale':
            failures.append(item)

old=root.parent/'p2_h_crystal_boundary_20261108/matrix'
for n in (40,80):
    d=old/f'crystal_n{n}';h=read(d/'coupled_history.csv');a=h[h.active==1].iloc[0]
    da=read(d/'donor_attempts.csv');g=da[da.sufficient==1].iloc[-1]
    activation.append(dict(case=f'crystal_n{n}',stage='historical_reference',time_s=a.time_s,h_um=a.h_m*1e6,
        c_initial=g.salt_transfer_kg/(g.salt_transfer_kg+g.water_transfer_kg),donor_depth_mm=g.used_depth_m*1000,
        donor_count=g.donor_count,transfer_salt_relative=abs(g.salt_transfer_residual_kg)/g.initial_salt_kg,
        transfer_water_relative=abs(g.water_transfer_residual_kg)/g.initial_water_kg))
for n in (40,80,160,320):
    p=root/f'diagnostic/historical{n}/coupled_history.csv';ref=old/f'pre_n{n}/coupled_history.csv'
    if p.exists() and (p.parent/'exit_status.txt').exists():
        same=p.read_bytes()==ref.read_bytes();checks.append(dict(test=f'observer_regression_{n}',status='PASS' if same else 'FAIL'));assert same
p=root/'correction/fixed_tests'
if (p/'undersaturated_projection/exit_status.txt').exists():
    a=read(p/'undersaturated_projection/coupled_history.csv');b=read(p/'undersaturated_implicit_surface/coupled_history.csv')
    same=np.max(np.abs(a.to_numpy()-b.to_numpy()))==0
    assert same;checks.append(dict(test='inactive_boundary_exact_regression',status='PASS'))
if (p/'projection40_regression/exit_status.txt').exists():
    assert (p/'projection40_regression/coupled_history.csv').read_bytes()==(old/'pre_n40/coupled_history.csv').read_bytes()
    checks.append(dict(test='retained_projection_exact_regression',status='PASS'))

p=root/'correction/fixed_coupled/fixed80_rollback'
if (p/'exit_status.txt').exists():
    events=read(p/'coupling_events.csv');h=read(p/'coupled_history.csv')
    assert events.event.str.startswith('rejected_injected rollback test').sum()==1
    assert (events.event=='accepted_activation_attempt').sum()==1
    assert (h.active.diff().fillna(0)>=0).all()
    active=h[h.active==1]
    assert active.pore_salt_transfer_kg.nunique()==1 and active.pore_water_transfer_kg.nunique()==1
    checks.append(dict(test='forced_rejection_accepted_inventory_rollback',status='PASS'))

for label,left,right,limit in [('40to80','fixed40','fixed80',.02),('80to160','fixed80','fixed160',.02),('160to320','fixed160','fixed320',.02),('dt160','fixed160','fixed_half160',.01),('dt320','fixed320','fixed_half320',.01)]:
    paths=[root/'correction/fixed_coupled'/n/'coupled_history.csv' for n in (left,right)]
    if not all(p.exists() and (p.parent/'exit_status.txt').exists() for p in paths):continue
    hs=[read(p) for p in paths];a,b=[h.iloc[-1] for h in hs]
    complete=all(h.time_s==2500 and h.active==1 for h in (a,b))
    for quantity in ('h_m','h_mass_m','top_solid_kg_m2','internal_solid_kg_m2','evaporation_kg_m2_s','epsilon_min'):
        difference=abs(a[quantity]-b[quantity])/abs(b[quantity]) if complete else np.nan
        primary=quantity in ('h_m','h_mass_m','evaporation_kg_m2_s')
        convergence.append(dict(comparison=label,quantity=quantity,relative_difference=difference,limit=limit,
            primary_gate=primary,status=('PASS' if difference<=limit else 'FAIL') if complete else 'FAIL_INCOMPLETE'))
    if complete:
        fs=[read(p.parent/'ale_diagnostics.csv') for p in paths]
        fs=[f[(f.status=='accepted')&(f.time_s==2500)] for f in fs]
        z=np.linspace(0,1,501);eps=[]
        for f in fs:eps.append(np.interp(z,(f.z_bottom_m+f.z_top_m)/(2*f.h_m),f.epsilon))
        convergence.append(dict(comparison=label,quantity='porosity_profile_linf_absolute',relative_difference=np.max(np.abs(eps[0]-eps[1])),limit=np.nan,primary_gate=False,status='DIAGNOSTIC'))

save(checks,'verification_checks.csv');save(convergence,'convergence.csv');save(activation,'activation_comparison.csv');save(budget,'independent_mass_budgets.csv');save(failures,'failed_cases.csv');save(independent_series,'independent_residual_series.csv')
print(json.dumps(dict(checks=checks,completed_budget_audits=len(budget),failed_cases=failures),indent=2))
