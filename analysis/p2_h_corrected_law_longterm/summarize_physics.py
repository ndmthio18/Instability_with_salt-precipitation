"""Traceable milestone, buffering and same-grid screening tables."""
from pathlib import Path
import numpy as np
import pandas as pd
R=Path(__file__).resolve().parent
s=pd.read_csv(R/'results/summary.csv');rows=[];screen=[]
for row in s.itertuples():
    d=pd.read_csv(R/'raw'/f'{row.id}_history.csv')
    from analyze import law,birth
    v=law(d.h,d.j,d.B,birth,row.beta,row.avmult)[0]
    error=abs(v/d.v-1);ix=np.flatnonzero((d.t>=1/row.jfactor)&(error>.01))
    last=d.iloc[-1];m0=d.solid.iloc[0]
    rows.append(dict(id=row.id,initial_solid=m0,final_solid=last.solid,final_dissolved=last.dissolved,
      final_h=last.h,top_cumulative=last.Qtop,internal_precip_cumulative=last.Qplus,
      dissolution_cumulative=last.Qminus,bottom_salt_cumulative=last.Qb,
      gross_dissolution_over_initial_solid=last.Qminus/m0,
      first_sample_above_1pct_s=d.t.iloc[ix[0]] if len(ix) else np.nan,
      preceding_sample_s=d.t.iloc[max(0,ix[0]-1)] if len(ix) else np.nan))
    if row.n==128 and row.power==1.5 and row.u>=0 and not row.id.startswith('p1'):
        ref=pd.read_csv(R/'raw/u0_n128_p1.5_r1e-08_dt30_history.csv')
        speed0=np.interp(last.t,ref.t,ref.v)
        screen.append(dict(id=row.id,t=last.t,u=row.u,ALE_v=last.v,matched_u0_v=speed0,
             relative_speed_difference=(last.v/speed0-1),diffusion_rate=last.Dminus,
             top_precip_rate=last.G,net_bottom_salt=last.Jb))
pd.DataFrame(rows).to_csv(R/'results/redistribution_summary.csv',index=False)
pd.DataFrame(screen).to_csv(R/'results/screening_comparison.csv',index=False)
parameters=s[['id','u','beta','avmult','jfactor','n','power','rtol']].copy()
parameters['initial_a_v_m_inv']=36*(1-.43)/25e-6*parameters.avmult
parameters['k_s_m_s']=parameters.beta*(1-.43)*1e-3
parameters['constant_j_kg_m2_s']=np.where(parameters.u>=0,1.78e-5*parameters.jfactor,np.nan)
parameters['initial_epsilon']=.43;parameters['D_s_m2_s']=1.3e-9
parameters['initial_D_eff_m2_s']=np.sqrt(.43)*1.3e-9
parameters['k_r_m_s']=.001;parameters['grain_m']=25e-6
parameters.to_csv(R/'results/dimensional_parameters.csv',index=False)
