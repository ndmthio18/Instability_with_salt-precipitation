"""Machine-readable scientific gates; finishing a run is only one condition."""
import json,os,pathlib
import numpy as np
import pandas as pd
assert os.environ.get('SLURM_JOB_ID')
root=pathlib.Path.cwd();out=root/'audit'
b=pd.read_csv(out/'independent_mass_budgets.csv');c=pd.read_csv(out/'convergence.csv');a=pd.read_csv(out/'activation_comparison.csv')
pre=pd.read_csv(out/'preactivation_convergence.csv');ab=pd.read_csv(out/'ablation_comparison.csv')
final=pd.read_csv(root/'correction/audit/fixed_coupled_summary.csv');main=final[final['case'].isin(['fixed40','fixed80','fixed160','fixed320'])]
result=dict(corrected_all_four_complete=bool((main.time_s==2500).all() and len(main)==4),
    adjacent_primary_gate='PASS' if (c[c.primary_gate].status=='PASS').all() else 'FAIL',
    maximum_combined_residual_all=float(b[['global_salt_relative','global_water_relative']].max().max()),
    maximum_combined_residual_corrected_coupled=float(b[b.stage=='correction/fixed_coupled'][['global_salt_relative','global_water_relative']].max().max()),
    maximum_activation_transfer_relative=float(a[['transfer_salt_relative','transfer_water_relative']].max().max()),
    final_coupled=main.to_dict(orient='records'),
    max_grid_height_relative=float(c[(c.primary_gate)&(c.comparison.str.contains('to'))].relative_difference.max()),
    max_dt_primary_relative=float(c[(c.primary_gate)&(c.comparison.str.startswith('dt'))].relative_difference.max()),
    ablation_temporal_checks=[],near_surface={},pre_surface_refinement={})
for name in ('eq0_n320_native','eq0_n320_replay_total'):
    x=ab[ab['case']==name].iloc[0];y=ab[ab['case']==name+'_half'].iloc[0]
    error=abs(x.internal_dynamic_g_m2-y.internal_dynamic_g_m2)/abs(y.internal_dynamic_g_m2)
    result['ablation_temporal_checks'].append(dict(case=name,internal_solid_relative=error,status='PASS' if error<.01 else 'FAIL'))
q=pd.read_csv(out/'diagnostic_mass_budgets.csv')
for family in pre.family.unique():
    p=pre[(pre.family==family)&(pre.quantity=='surface_g_m2')]
    result['pre_surface_refinement'][family]=float(p.relative_difference.max())
for stage,names in [('diagnostic',['historical160','historical320']),('correction/fixed_pre',['fixed160','fixed320'])]:
    for name in names:
        r=q[(q.stage==stage)&(q['case']==name)&(q.depth_m==.01)].iloc[0]
        result['near_surface'][name]=r.to_dict()
pressure=concentration=0;minimum=1
for n in (40,80,160,320):
    d=pd.read_csv(root/f'correction/fixed_coupled/fixed{n}/donor_profiles.csv');d=d[d.extracted_kg>0]
    pressure=max(pressure,float((d.p_after_pa-d.p_before_pa).abs().max()))
    concentration=max(concentration,float((d.c_after-d.c_before).abs().max()))
    minimum=min(minimum,float(d.sw_after.min()))
result['donor_perturbations']=dict(max_pressure_change_pa=pressure,max_concentration_change=concentration,min_retained_sw=minimum)
result['decision']='PASS numerical exploratory gates for declared 1D equations; NOT TESTED physical validation and asymptotic ALE convergence'
(out/'scientific_gates.json').write_text(json.dumps(result,indent=2,allow_nan=False))
print(json.dumps(result,indent=2,allow_nan=False))
