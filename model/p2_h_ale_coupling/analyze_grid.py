"""Compute-node diagnostic audit; failures are data, never filtered out."""
import csv, hashlib, json, pathlib, sys, os
assert os.environ.get("SLURM_JOB_ID")
import numpy as np
import pandas as pd
root=pathlib.Path.cwd(); out=root/'audit'; out.mkdir(exist_ok=True)
def save(rows,name): pd.DataFrame(rows).to_csv(out/name,index=False)
def historical():
    campaign=root.parent if root.name=='correction' else root
    old=campaign.parent/'p2_h_crystal_boundary_20261108/matrix'; rows=[]
    for n in (40,80,160,320):
        p=old/f'pre_n{n}'; h=pd.read_csv(p/'coupled_history.csv'); a,b=h.iloc[0],h.iloc[-1]
        f=pd.read_csv(p/'soil_fields.csv'); z=f[f.time_s==0].groupby('y_m').mean(numeric_only=True).sort_index()
        y=z.index.to_numpy(); bounds=np.r_[0,(y[:-1]+y[1:])/2,.2]; dz=np.diff(bounds)
        item=dict(n=n,time_s=b.time_s,surface_solid_g_m2=b.top_solid_kg_m2*1000,
                  soil_dissolved_change_g_m2=(b.soil_salt_kg-a.soil_salt_kg)/.001*1000,
                  initial_soil_salt_kg=a.soil_salt_kg,final_soil_salt_kg=b.soil_salt_kg,
                  budget_kg=b.salt_residual_kg,top_cell_m=y[-1]-y[-2])
        for t in (0,2500):
            z=f[f.time_s==t].groupby('y_m').mean(numeric_only=True).sort_index()
            c=z.mass_fraction.to_numpy(); m=.41*z.saturation.to_numpy()*z.density_kg_m3.to_numpy()
            for depth in (.0001,.0005,.001,.002,.005,.01):
                overlap=np.maximum(0,bounds[1:]-np.maximum(bounds[:-1],.2-depth))
                item[f'salt_{depth:g}m_{t}s_g_m2']=np.sum(m*c*overlap)*1000
                item[f'excess_{depth:g}m_{t}s_g_m2']=np.sum(m*np.maximum(0,(c-.2649408697676898)/(1-.2649408697676898))*overlap)*1000
        rows.append(item)
    save(rows,'historical_mass_location.csv')
    print(pd.DataFrame(rows)[['n','surface_solid_g_m2','soil_dissolved_change_g_m2','excess_0.01m_2500s_g_m2']].to_string(index=False))

def audit(stage):
    if not (root/stage).exists():
        print('NOT TESTED: stage directory absent:',stage)
        return
    summaries=[]; cuts=[]; regression=[]
    for p in sorted((root/stage).iterdir()):
        if not p.is_dir() or not (p/'exit_status.txt').exists(): continue
        h=pd.read_csv(p/'coupled_history.csv'); a,b=h.iloc[0],h.iloc[-1]; active=h[h.active==1]
        grid=pd.read_csv(p/'grid_profiles.csv'); init=grid[grid.time_s==0]
        d=np.sort(init.depth_m.unique()); item=dict(case=p.name,stage=stage,exit_code=int((p/'exit_status.txt').read_text()),
            layers=len(d)-1,top_cell_m=d[1],time_s=b.time_s,
            first_precip_s=h.loc[h.top_solid_kg_m2>1e-15,'time_s'].min(),
            threshold_s=h.loc[h.h_mass_m>=3.5e-6,'time_s'].min(),
            activation_s=active.time_s.min(),surface_g_m2=b.top_solid_kg_m2*1000,
            height_um=b.h_m*1e6,hmass_um=b.h_mass_m*1e6,internal_g_m2=b.internal_solid_kg_m2*1000,
            epsilon_min=b.epsilon_min,evaporation=b.evaporation_kg_m2_s,
            salt_relative=h.salt_residual_kg.abs().max()/a.soil_salt_kg,
            water_relative=h.water_residual_kg.abs().max()/a.soil_water_kg)
        assert max(item['salt_relative'],item['water_relative'])<=1e-8
        assert abs(init.box_volume_m3.sum()-.0002)<1e-14
        q=pd.read_csv(p/'grid_budgets.csv'); q=q[q.active==0]
        if len(q):
            boundary=q.get("boundary_precip_kg_m2_s",0)*.001
            reconstruct=q.get("boundary_reconstruction_salt_kg_s",0)
            net=q.adv_in_kg_s+q.diff_in_kg_s+reconstruct-boundary
            rs=q.salt_transport_kg-q.salt_before_kg-q.dt_s*net
            rw=q.water_transport_kg-q.water_before_kg-q.dt_s*(q.liquid_in_kg_s-q.adv_in_kg_s-q.diff_in_kg_s-reconstruct-.001*q.evap_kg_m2_s)
            item.update(cut_salt_step_relative=rs.abs().max()/a.soil_salt_kg,cut_water_step_relative=rw.abs().max()/a.soil_water_kg)
            assert max(item['cut_salt_step_relative'],item['cut_water_step_relative'])<1e-8,(p.name,item)
            for depth,g in q.groupby('depth_m'):
                first,last=g.iloc[0],g.iloc[-1]
                cuts.append(dict(case=p.name,depth_m=depth,time_s=last.time_s,
                    dissolved_change_g_m2=(last.salt_after_kg-first.salt_before_kg)*1e6,
                    surface_g_m2=last.surface_solid_kg*1e6,
                    regional_projection_g_m2=g.projection_salt_kg.sum()*1e6,
                    adv_integral_g_m2=(g.dt_s*g.adv_in_kg_s).sum()*1e6,
                    diff_integral_g_m2=(g.dt_s*g.diff_in_kg_s).sum()*1e6,
                    excess_g_m2=last.supersaturation_excess_kg*1e6,
                    closure_kg=last.salt_after_kg-first.salt_before_kg+g.projection_salt_kg.sum()-(g.dt_s*(g.adv_in_kg_s+g.diff_in_kg_s+g.get("boundary_reconstruction_salt_kg_s",0)-.001*g.get("boundary_precip_kg_m2_s",0))).sum(),
                    boundary_precip_integral_g_m2=(g.dt_s*g.get("boundary_precip_kg_m2_s",0)).sum()*1000,
                    reconstructed_cut_integral_g_m2=(g.dt_s*g.get("boundary_reconstruction_salt_kg_s",0)).sum()*1e6))
        if p.name in ('historical40','historical80','historical160','historical320'):
            old=root.parent/f'p2_h_crystal_boundary_20261108/matrix/pre_n{p.name[10:]}/coupled_history.csv'
            same=(p/'coupled_history.csv').read_bytes()==old.read_bytes()
            regression.append(dict(case=p.name,whole_history_byte_equal=same))
        if len(active):
            da=pd.read_csv(p/'donor_attempts.csv'); good=da[da.sufficient==1].iloc[-1]
            item.update(c_initial=good.salt_transfer_kg/(good.salt_transfer_kg+good.water_transfer_kg),donor_depth_m=good.used_depth_m,
                        transfer_salt_relative=abs(good.salt_transfer_residual_kg)/a.soil_salt_kg,
                        transfer_water_relative=abs(good.water_transfer_residual_kg)/a.soil_water_kg)
        summaries.append(item)
    save(summaries,stage+'_summary.csv'); save(cuts,stage+'_mass_budgets.csv'); save(regression,stage+'_regression.csv')
    print(pd.DataFrame(summaries).to_string(index=False))

historical()
for stage in sys.argv[1:]: audit(stage)
