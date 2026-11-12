"""Machine-readable acceptance gates; never turn a failed gate into PASS."""
import json
from pathlib import Path
import numpy as np
import pandas as pd
R=Path(__file__).resolve().parent/'results'
s=pd.read_csv(R/'summary.csv');c=pd.read_csv(R/'convergence.csv');f=pd.read_csv(R/'frozen_bvp_terms.csv')
stats={}
for name,g in f.groupby('id'):
    jf=float(s.loc[s.id==name,'jfactor'].iloc[0])
    d=g[g.t>=1/jf].copy()
    if d.empty:continue
    for model in ['A','Bmean','Bupper','Cbvp']:d[model+'_error']=abs(d[model]/d.ALE-1)
    last=d.iloc[-1]
    stats[name]={k:float(last[k]) for k in ['t','ALE','A_error','Bmean_error','Bupper_error','Cbvp_error','concentration_storage_to_reaction']}
budget=max(s[k].abs().max() for k in ['salt_res','water_res','diss_res','local_res'])
out={'baseline':'c6fb6277c4cdf1689d80f4bbdda7a421060d874a','completed_cases':len(s),
 'budget_max_kg_m2':float(budget),'budget_pass':bool(budget<2e-9),
 'all_target_durations_reached':bool((s.status=='end').all()),
 'framework_stops':int(s.status.str.startswith('solid_framework').sum()),
 'pore_filling_stops':int(s.status.str.startswith('pore_filling').sum()),
 'target_completions':int((s.status=='end').sum()),
 'solver_stops':int(s.status.str.startswith('solver_stop').sum()),
 'frozen_bvp_and_effective_properties':stats}
ref=c[(c.u==0)&(c.n==256)&(c.power==2)]
if len(ref):out['u0_n256_vs_finest_max_velocity_difference']=float(ref.velocity_max_relative.iloc[0])
(R/'acceptance.json').write_text(json.dumps(out,indent=2)+'\n',encoding='utf-8')
print(json.dumps(out,indent=2))
