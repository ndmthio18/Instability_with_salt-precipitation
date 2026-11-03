"""Numerical diagnostic figure; saved data only."""
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from wp8_v2_regression import c
from analyze_wp8 import rows
O=c.OUT
fig,axs=plt.subplots(2,2,figsize=(10,7),layout='constrained')
r=rows(c.ROOT/'analysis/wp8_strict_background_20261103/failed_step_iterations.csv');a=axs[0,0]
for key,limit,label in [('salt_residual',2e-20,'Salt / old salt limit'),('water_residual',2e-18,'Water / old water limit')]:a.semilogy([int(x['iteration']) for x in r],[float(x[key])/limit for x in r],'.-',label=label)
a.axhline(1,color='black',lw=.8,ls='--');a.set(xlabel='Newton iteration',ylabel='Residual / old limit',title='A  Original rejected step 4414');a.legend(fontsize=8)
a=axs[0,1];r=rows(O/'late_legacy_window/convergence_v2.csv');r=[x for x in r if x['step']=='10']
a.semilogy([int(x['iteration']) for x in r],[float(x['salt_max'])/2e-20 for x in r],'.-',label='Legacy salt / old limit')
a.semilogy([int(x['iteration']) for x in r],[float(x['scaled_nonlinear'])/(32*np.finfo(float).eps) for x in r],'.-',label='Scaled residual / new limit')
a.axhline(1,color='black',lw=.8,ls='--');a.set(xlabel='Newton iteration',ylabel='Residual / criterion',title='B  Nearby replay: legacy step 10');a.legend(fontsize=8)
a.text(.40,.45,'Other acceptance gates also apply',transform=a.transAxes,fontsize=7.5)
a=axs[1,0];eq=json.loads((O/'equation_diagnostics.json').read_text())
for col in [2,200,398,400]:
 rr=sorted([x for x in eq['fd_columns'] if x['column']==col],key=lambda x:x['h']);a.loglog([x['h'] for x in rr],[x['normwise_relative'] for x in rr],'.-',label=f'Pressure column {col}')
a.set(xlabel='Central difference step (Pa)',ylabel='Column relative difference',title='C  Jacobian difference-step audit');a.legend(fontsize=8)
a=axs[1,1];names=['early_scaled','nearT_scaled','late_scaled','late_legacy_window'];vals=[json.loads((O/n/'validation.json').read_text())['metrics']['max_salt_closure_m'] for n in names]
a.bar(range(4),vals,color=['#286f9e','#286f9e','#286f9e','#999999']);a.set_yscale('log');a.axhline(1e-12,ls='--',color='black',label='Unchanged salt budget gate');a.set_xticks(range(4),['Early','Near T','Late v2','Late legacy']);a.set(ylabel='Max regional salt closure (m)',title='D  Independent short-window budgets',ylim=(1e-18,3e-12));a.legend(fontsize=8)
for a in axs.flat:a.spines[['top','right']].set_visible(False);a.grid(alpha=.15)
fig.savefig(O/'Figure_Newton_diagnostics.png',dpi=180);fig.savefig(O/'Figure_Newton_diagnostics.pdf');plt.close(fig)
