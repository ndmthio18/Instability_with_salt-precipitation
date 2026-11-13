"""Recompute all tables from accepted histories. No analytical input to the ALE."""
from pathlib import Path
import json
import numpy as np
import pandas as pd
ROOT=Path(__file__).resolve().parent
R=ROOT/'results'; R.mkdir(exist_ok=True)
CEQ=.264940869768; RHO=1202.168577; CR=2165.; EPS=.43; D=1.3e-9; KR=.001
MD=EPS**1.5*D/(1-CEQ); AV=36*(1-EPS)/25e-6; ELL=np.sqrt(MD/(AV*KR))
X=2*(1-EPS)*KR*ELL/MD; PHI=(1-EPS)*(1+X)/X

def enrich(h):
 h=h.copy(); H=h.hc/ELL;sech=2*np.exp(-H)/(1+np.exp(-2*H))
 h['Sref']=CEQ*h.j/(1-CEQ);h['vref']=h.Sref/(CR*PHI)
 h['vt_corr']=h.Sref/(CR*(1-EPS))*X*(1-h.B/(CEQ*h.j)*sech)/(X+np.tanh(H))
 h['P_corr']=(CEQ*h.j-h.B)/(1-CEQ);h['vnet_corr']=h.P_corr/(CR*PHI)
 h['vb_diagnostic']=h.vt_corr-h.vnet_corr
 # Integrate the unmodified analytical speed on the same physical clock.
 # Include all forcing knots so piecewise-linear forcing is integrated exactly.
 if h.j.nunique()>1 or h.B.nunique()>1:
  forcing=pd.read_csv(ROOT/'forcing.csv');ft=forcing.iloc[:,0].to_numpy()
  tt=np.unique(np.r_[h.t.to_numpy(),ft[(ft>=0)&(ft<=h.t.max())]])
  jj=np.interp(tt,ft,forcing.iloc[:,2]);bb=np.interp(tt,ft,forcing.iloc[:,3])
  vv=(CEQ*jj-bb)/((1-CEQ)*CR*PHI)
  integral=np.r_[0,np.cumsum((vv[1:]+vv[:-1])*.5*np.diff(tt))]
  h['hc_corr']=h.hc.iloc[0]+np.interp(h.t,tt,integral)
 else:h['hc_corr']=h.hc.iloc[0]+h.vnet_corr.iloc[0]*h.t
 h['height_error_relative']=(h.hc-h.hc_corr)/h.hc
 h['delta_top']=h.vt-h.vt_corr;h['delta_bottom']=h.vb-h.vb_diagnostic;h['delta_net']=h.vnet-h.vnet_corr
 h['phi_mean']=h.solid/(CR*h.hc);h['solid_rate']=h.G-h.Db+h.Rplus-h.Dminus
 # This algebraic identity is explanatory, not an independent conservation test.
 h['phi_rate_identity']=(h.solid_rate/CR-h.phi_mean*h.vnet)/h.hc
 h['mass_term']=h.solid_rate/(CR*h.phi_mean)
 h['supply_term']=(h.Jsoil-h.P_corr)/(CR*PHI)
 h['storage_term']=-(h.Jsoil-h.solid_rate)/(CR*PHI)
 h['density_term']=h.solid_rate/CR*(1/h.phi_mean-1/PHI)
 h['structure_term']=-h.hc*h.phi_rate_identity/h.phi_mean
 if len(h)>2:
  h['phi_rate_fd']=np.gradient(h.phi_mean,h.t,edge_order=2)
  h['identity_fd_residual']=h.vnet-(h.mass_term-h.hc/h.phi_mean*h.phi_rate_fd)
 return h

def read(case):
 return enrich(pd.read_csv(ROOT/'raw'/f'{case}_history.csv'))

def compare(a,b,kind,keys):
 t=np.intersect1d(a.t,b.t);aa=a.set_index('t').loc[t];bb=b.set_index('t').loc[t]
 assert len(t)>0
 row=dict(kind=kind,end_s=float(t[-1]),n_times=len(t))
 for key in keys:
  d=np.abs(aa[key]-bb[key]);scale=bb.vref if key in ['vt','vb','vnet'] else bb.Sref if key in ['Fplus','Fminus'] else np.maximum(np.abs(bb[key]),1e-9)
  row[key+'_max_abs']=float(d.max());row[key+'_max_normalized']=float((d/scale).max())
 return row

def main():
 summaries=[];allh=[];meta={}
 for config in ['parallel_cases.csv','fast_p01_cases.csv','final_cases.csv','p2_cases.csv','diagnostic_cases.csv','p0_cases.csv','early_cases.csv','retry_cases.csv','production_cases.csv']:
  path=ROOT/config
  if path.exists():
   for c in pd.read_csv(path).to_dict('records'):meta[c['id']]=c
 for p in sorted((ROOT/'raw').glob('*_history.csv')):
  case=p.name[:-12];sfile=ROOT/'raw'/f'{case}_status.txt'
  if not sfile.exists():continue
  status=dict(line.split('=',1) for line in sfile.read_text().splitlines());h=read(case);h['id']=case;allh.append(h)
  row=dict(meta.get(case,{}));row['id']=case;row.update(status);row['samples']=len(h)
  for k in ['vt','vb','vnet','zt','zb','hc','phi_mean','core_phi','phib','vt_corr','vnet_corr','delta_top','delta_bottom','delta_net']:
   row['final_'+k]=float(h[k].iloc[-1])
  for k in ['salt_res','water_res','solid_res','liquid_salt_res','crust_salt_res','liquid_water_res','crust_water_res','water_jump','contact_res']:
   row['max_'+k]=float(h[k].abs().max())
  summaries.append(row)
 summary=pd.DataFrame(summaries);summary.to_csv(R/'summary.csv',index=False)
 (R/'histories').mkdir(exist_ok=True)
 for h in allh:h.to_csv(R/'histories'/(str(h.id.iloc[0])+'.csv.gz'),index=False,compression={'method':'gzip','mtime':0})
 conv=[]
 prefix='fast_' if (ROOT/'raw/fast_p0_ld_x008_r1e-10_status.txt').exists() else 'setup_' if (ROOT/'raw/setup_p0_ld_x008_r1e-10_status.txt').exists() else ''
 for lab in ['u06','u09','x008']:
  a=f'{prefix}p0_ld_{lab}_r1e-09';b=f'{prefix}p0_ld_{lab}_r1e-10'
  conv.append(dict(id=a,reference=b,**compare(read(a),read(b),'P0 tolerance',['vt','vb','vnet','zt','zb','hc','Fplus','Fminus'])))
 pd.DataFrame(conv).to_csv(R/'temporal_convergence.csv',index=False)
 flux=[];windows=[]
 for lab in ['u06','x008']:
  for gap in [1e-10,1e-9,1e-8]:
   for n in [128,256]:
    a=f'{prefix}early_{lab}_n{n}_b{gap:g}';b=f'{prefix}early_{lab}_n{2*n}_b{gap:g}'
    flux.append(dict(id=a,reference=b,gap=gap,**compare(read(a),read(b),'spatial',['Fplus','Fminus','Cb','vt','vb','vnet','hc'])))
  for n in [128,256,512]:
   a=f'{prefix}early_{lab}_n{n}_b1e-09';b=f'{prefix}early_{lab}_n{n}_tight'
   flux.append(dict(id=a,reference=b,gap=1e-9,**compare(read(a),read(b),'temporal',['Fplus','Fminus','Cb','vt','vb','vnet','hc'])))
  names=[f'{prefix}early_{lab}_n{n}_b1e-09' for n in [128,256,512]]+[f'{prefix}early_{lab}_n{n}_tight' for n in [128,256,512]]
  hist=[read(x) for x in names];tt=hist[0].t.to_numpy()
  for h in hist[1:]:tt=np.intersect1d(tt,h.t)
  errors=np.zeros(len(tt));details=[]
  for ai,bi in [(0,1),(1,2),(0,3),(1,4),(2,5)]:
   a=hist[ai].set_index('t').loc[tt];b=hist[bi].set_index('t').loc[tt]
   for key in ['Fplus','Fminus']:errors=np.maximum(errors,np.abs(a[key]-b[key])/b.Sref)
  after=np.maximum.accumulate(errors[::-1])[::-1];passed=np.where(after<.01)[0]
  complete=all(h.t.iloc[-1]>=100 for h in hist)
  windows.append(dict(case=lab,common_end_s=float(tt[-1]),all_complete_100s=complete,earliest_prefix_pass_s=float(tt[passed[0]]) if len(passed) else None,certified_earliest_s=float(tt[passed[0]]) if len(passed) and complete else None,max_normalized_error=float(errors.max())))
  pd.DataFrame({'t':tt,'max_normalized_flux_error':errors}).to_csv(R/f'flux_errors_{lab}.csv',index=False)
 pd.DataFrame(flux).to_csv(R/'spatial_temporal_flux_convergence.csv',index=False)
 (R/'flux_windows.json').write_text(json.dumps(windows,indent=2))
 diagnostics=[]
 for p in (ROOT/'raw').glob('*_diagnostics.csv'):
  try:
   d=pd.read_csv(p)
   if len(d):row=d.iloc[-1].to_dict();row['id']=p.name[:-16];diagnostics.append(row)
  except pd.errors.EmptyDataError:pass
 pd.DataFrame(diagnostics).to_csv(R/'solver_diagnostics.csv',index=False)
 (R/'constants.json').write_text(json.dumps(dict(ceq=CEQ,rho=RHO,rho_cr=CR,epsilon_birth=EPS,ell=ELL,X=X,phi_core=PHI),indent=2))
 print(summary[['id','status','last_s']].to_string(index=False));print(pd.DataFrame(conv).to_string(index=False));print(windows)
if __name__=='__main__':main()
