"""Full-field failure location and conservative intervention ledgers."""
import json,os,pathlib
import numpy as np
import pandas as pd
assert os.environ.get('SLURM_JOB_ID')
root=pathlib.Path.cwd();out=root/'audit';out.mkdir(exist_ok=True)
def save(rows,name):pd.DataFrame(rows).to_csv(out/name,index=False)
results=[]
for n in (160,320):
    d=root/f'ale_diagnostic/historical{n}_ale';f=pd.read_csv(d/'ale_diagnostics.csv');h=pd.read_csv(d/'coupled_history.csv');h=h[h.active==1]
    a,b=h.iloc[0],h.iloc[-1];accepted=f[f.status=='accepted'];first=accepted.iloc[:64];last=accepted.iloc[-64:]
    tiny=last.loc[last.epsilon.idxmin()]
    rejected=f[f.status=='bound_rejected'];item=dict(n=n,activation_s=a.time_s,last_accepted_s=b.time_s,
        duration_s=b.time_s-a.time_s,min_cell=int(tiny.i),min_cell_z_bottom_um=tiny.z_bottom_m*1e6,
        min_cell_z_top_um=tiny.z_top_m*1e6,epsilon_min=tiny.epsilon,c_at_min=tiny.c,reaction_at_min=tiny.reaction_kg_m3_s,
        initial_dissolved_salt_g_m2=a.crust_salt_kg*1e6,initial_water_g_m2=a.crust_water_kg*1e6,
        initial_excess_g_m2=(a.crust_salt_kg-.2649408697676898/(1-.2649408697676898)*a.crust_water_kg)*1e6,
        salt_input_g_m2=b.interface_salt_kg_m2*1000,
        water_input_g_m2=(b.interface_liquid_kg_m2-b.interface_salt_kg_m2)*1000,
        dissolved_salt_change_g_m2=(b.crust_salt_kg-a.crust_salt_kg)*1e6,
        top_solid_change_g_m2=(b.top_solid_kg_m2-a.top_solid_kg_m2)*1000,
        internal_solid_g_m2=b.internal_solid_kg_m2*1000)
    bottom=accepted[accepted.i==0].set_index('time_s').bottom_c
    adv=(h.interface_liquid_kg_m2.diff().iloc[1:].to_numpy()*bottom.loc[h.time_s.iloc[1:]].to_numpy()).sum()
    item.update(advective_input_g_m2=adv*1000,diffusive_input_g_m2=(b.interface_salt_kg_m2-adv)*1000)
    first.to_csv(out/f'ale{n}_birth.csv',index=False);last.to_csv(out/f'ale{n}_last_accepted.csv',index=False)
    if len(rejected):
        firstrej=rejected.iloc[:64];firstrej.to_csv(out/f'ale{n}_first_rejected.csv',index=False)
        r=firstrej.loc[firstrej.epsilon.idxmin()]
        item.update(first_rejected_time_s=r.time_s,first_rejected_min_porosity=r.epsilon,
                    first_rejected_min_cell=int(r.i),first_rejected_concentration=r.c,
                    rejected_min_c=firstrej.c.min(),rejected_max_c=firstrej.c.max())
    results.append(item)
save(results,'ale_failure_budget.csv')
ablations=[]
for d in sorted((root/'ablations').glob('eq*')):
    if not (d/'exit_status.txt').exists():continue
    h=pd.read_csv(d/'history.csv');a=h[h.status=='accepted'].iloc[0];b=h[h.status=='accepted'].iloc[-1]
    m=pd.read_csv(d/'initial.csv').iloc[0];cmd=json.loads((d/'command.json').read_text())
    item=dict(case=d.name,exit_code=int((d/'exit_status.txt').read_text()),end_s=b.time_s,requested_s=float(cmd[-2]),
        h0_um=m.h_initial*1e6,h_original_um=m.h_original*1e6,epsilon_initial=m.epsilon_initial,c_initial=m.c_initial,
        initialization_conversion_g_m2=m.initial_internal_conversion*1000,
        initial_salt_residual=m.salt_residual,initial_water_residual=m.water_residual,
        h_final_um=b.h_m*1e6,epsilon_min=b.epsilon_min,top_added_g_m2=b.top_added_kg_m2*1000,
        internal_dynamic_g_m2=b.internal_added_kg_m2*1000,
        total_internal_g_m2=(b.internal_added_kg_m2+m.initial_internal_conversion)*1000,
        salt_input_g_m2=b.salt_input_kg_m2*1000,
        water_input_g_m2=(b.liquid_input_kg_m2-b.salt_input_kg_m2)*1000,
        salt_residual_max=h[h.status=='accepted'].salt_residual_kg_m2.abs().max(),
        water_residual_max=h[h.status=='accepted'].water_residual_kg_m2.abs().max())
    assert max(item['salt_residual_max'],item['water_residual_max'])<1e-10
    item['status']='PASS_EXECUTION' if item['exit_code']==0 and abs(item['end_s']-item['requested_s'])<1e-9 else 'FAIL_INCOMPLETE'
    ablations.append(item)
save(ablations,'ablation_comparison.csv')
print(pd.DataFrame(results).to_string(index=False));print(pd.DataFrame(ablations).to_string(index=False))
