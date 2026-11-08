"""Physical-region budgets and convergence, including nonlocal projection."""
import os,pathlib
import numpy as np
import pandas as pd
assert os.environ.get('SLURM_JOB_ID')
root=pathlib.Path.cwd();out=root/'audit'
rows=[]
for prefix,stage in [('', 'diagnostic'),('correction/','fixed_pre')]:
    q=pd.read_csv(root/prefix/'audit'/(stage+'_mass_budgets.csv'));q['stage']=prefix+stage
    q['nonlocal_projection_g_m2']=q.surface_g_m2-q.regional_projection_g_m2-q.boundary_precip_integral_g_m2
    q['net_cut_integral_g_m2']=q.adv_integral_g_m2+q.diff_integral_g_m2+q.reconstructed_cut_integral_g_m2
    q['combined_region_surface_residual_g_m2']=q.dissolved_change_g_m2+q.surface_g_m2-q.net_cut_integral_g_m2-q.nonlocal_projection_g_m2
    rows.append(q)
budgets=pd.concat(rows,ignore_index=True);budgets.to_csv(out/'diagnostic_mass_budgets.csv',index=False)
comparisons=[]
groups=[('original projection','diagnostic',[f'historical{n}' for n in (40,80,160,320)]),
        ('corrected original','correction/fixed_pre',[f'fixed{n}' for n in (40,80,160,320)]),
        ('fixed map A projection','diagnostic',[f'familyA{n}' for n in (40,80,160,320,640)]),
        ('fixed map A corrected','correction/fixed_pre',['fixed80','fixedA160','fixedA320','fixedA640']),
        ('fixed surface B projection','diagnostic',[f'familyB{n}_top80' for n in (40,80,160,320)]),
        ('fixed surface B corrected','correction/fixed_pre',[f'fixedB{n}_top80' for n in (40,80,160,320)]),
        ('surface refine B projection','diagnostic',[f'familyB160_top{n}' for n in (40,80,160)]),
        ('surface refine B corrected','correction/fixed_pre',[f'fixedB160_top{n}' for n in (40,80,160)]),
        ('old dt320','diagnostic',['historical320','historical320_half']),
        ('corrected dt160','correction/fixed_pre',['fixed160','fixed_half160']),
        ('corrected dt320','correction/fixed_pre',['fixed320','fixed_half320'])]
for family,stage,names in groups:
    for left,right in zip(names,names[1:]):
        a=budgets[(budgets.stage==stage)&(budgets['case']==left)].set_index('depth_m')
        b=budgets[(budgets.stage==stage)&(budgets['case']==right)].set_index('depth_m')
        if a.empty or b.empty:continue
        assert (a.time_s==2500).all() and (b.time_s==2500).all()
        for depth in b.index:
            for field in ('surface_g_m2','dissolved_change_g_m2','excess_g_m2','net_cut_integral_g_m2'):
                x,y=a.loc[depth,field],b.loc[depth,field];difference=abs(x-y)
                comparisons.append(dict(family=family,left=left,right=right,depth_m=depth,quantity=field,
                    coarse=x,fine=y,absolute_difference_g_m2=difference,
                    relative_difference=difference/abs(y) if abs(y)>1e-10 else np.nan))
pd.DataFrame(comparisons).to_csv(out/'preactivation_convergence.csv',index=False)
print('Maximum cumulative region budget discrepancy (kg/m2):',budgets.combined_region_surface_residual_g_m2.abs().max()/1000)
