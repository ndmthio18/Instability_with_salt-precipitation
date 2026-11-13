from analyze import *

def final_case(lab,n=512):
 if lab in ['u00','u03']:return f'fast_production_{lab}_n{n}'
 if lab in ['u06','u09'] and n==256:return f'fast_p2_{lab}_n256_k1'
 return f'fast_grid_{lab}_n{n}'

def main():
 grid=[];summaries=[];identities=[];kinetics=[];bridge=[]
 for lab in ['u00','u03','u06','u09','x008']:
  names=[final_case(lab,n) for n in [256,512]]
  if not all((ROOT/'raw'/f'{c}_status.txt').exists() for c in names):continue
  a,b=map(read,names)
  row=dict(case=lab,coarse=names[0],fine=names[1],**compare(a,b,'long-time grid',['vt','vb','vnet','zt','zb','hc','Fplus','Fminus']))
  row['core_phi_final_abs_difference']=float(abs(a.core_phi.iloc[-1]-b.core_phi.iloc[-1]));row['mean_phi_final_abs_difference']=float(abs(a.phi_mean.iloc[-1]-b.phi_mean.iloc[-1]));grid.append(row);h=b;last=h.iloc[-1];late=h[h.t>=.8*h.t.max()]
  row={'case':lab,'id':names[1],'end_h':last.t/3600}
  for k in ['vt','vb','vnet','vt_corr','vnet_corr','delta_top','delta_bottom','delta_net','zt','zb','hc','hc_corr','height_error_relative','phi_mean','core_phi','core_std','phib','solid','crust_salt','liquid_salt','mass_term','structure_term','supply_term','storage_term','density_term','solid_rate','P_corr','Jsoil']:
   row['final_'+k]=float(last[k])
  row['final_top_error_percent']=100*last.delta_top/last.vt
  row['final_net_error_percent']=100*last.delta_net/last.vnet
  row['final_net_error_reference_percent']=100*last.delta_net/last.vref
  row['full_top_error_reference_percent']=float((abs(h.delta_top)/h.vref*100).max())
  row['full_net_error_reference_percent']=float((abs(h.delta_net)/h.vref*100).max())
  row['late_core_drift_relative']=float((late.core_phi.max()-late.core_phi.min())/late.core_phi.mean())
  row['late_core_spread_relative']=float((late.core_std/late.core_phi).max())
  row['final_H']=float(last.hc/ELL)
  row['final_core_theory_relative_difference']=float(last.core_phi/PHI-1)
  row['final_core_diffusivity_to_birth']=float(((1-last.core_phi)/EPS)**1.5)
  row['final_core_area_to_birth']=float(last.core_phi/(1-EPS))
  row['final_cancellation_fraction']=float(1-abs(last.delta_net)/(abs(last.delta_top)+abs(last.delta_bottom)))
  summaries.append(row)
  for _,r in h.iterrows():
   identities.append(dict(case=lab,t=r.t,mass_term=r.mass_term,structure_term=r.structure_term,vnet=r.vnet,phi_rate_identity=r.phi_rate_identity,phi_rate_fd=r.get('phi_rate_fd',np.nan),fd_residual=r.get('identity_fd_residual',np.nan),supply_term=r.supply_term,storage_term=r.storage_term,density_term=r.density_term))
 for lab in ['u06','u09']:
  base=read(f'fast_p2_{lab}_n256_k1')
  for kap in [.5,2]:
   c=f'fast_p2_{lab}_n256_k{kap:g}'
   if not (ROOT/'raw'/f'{c}_history.csv').exists():continue
   h=read(c);r=compare(h,base,'kinetic sensitivity',['vt','vb','vnet','zt','zb','hc'])
   r.update(case=lab,kappa=kap,final_thickness_delta=float(h.hc.iloc[-1]-base.hc.iloc[-1]),final_net_delta=float(h.vnet.iloc[-1]-base.vnet.iloc[-1]))
   delta=np.abs((h.vnet-h.vnet_corr)-(base.vnet-base.vnet_corr));r['max_net_law_error_sensitivity_reference']=float((delta/base.vref).max())
   kinetics.append(r)
 for lab in ['u06','u09','x008']:
  for rt in ['1e-09','1e-10']:
   a=f'fast_p0_ld_{lab}_r{rt}';b=f'setup_p0_ld_{lab}_r{rt}'
   bridge.append(dict(id=a,reference=b,**compare(read(a),read(b),'power evaluation bridge',['vt','vb','vnet','zt','zb','hc','Fplus','Fminus'])))
 pd.DataFrame(grid).to_csv(R/'long_time_grid_convergence.csv',index=False)
 hierarchy=list(grid)
 for lab in ['u06','u09','x008']:
  ca=f'fast_p0_ld_{lab}_r1e-09';cb=final_case(lab,256)
  if (ROOT/'raw'/f'{cb}_status.txt').exists():hierarchy.append(dict(case=lab,coarse=ca,fine=cb,**compare(read(ca),read(cb),'128 to 256',['vt','vb','vnet','zt','zb','hc','Fplus','Fminus'])))
 pd.DataFrame(hierarchy).to_csv(R/'long_time_grid_hierarchy.csv',index=False)
 pd.DataFrame(summaries).to_csv(R/'growth_validation.csv',index=False)
 pd.DataFrame(identities).to_csv(R/'geometric_identity.csv.gz',index=False,compression={'method':'gzip','mtime':0})
 pd.DataFrame(kinetics).to_csv(R/'kinetic_sensitivity.csv',index=False)
 pd.DataFrame(bridge).to_csv(R/'implementation_bridge.csv',index=False)
 if (ROOT/'raw/fast_production_u06_n512_tight_history.csv').exists():
  a=read('fast_grid_u06_n512');b=read('fast_production_u06_n512_tight')
  pd.DataFrame([compare(a,b,'512-cell temporal',['vt','vb','vnet','zt','zb','hc','Fplus','Fminus'])]).to_csv(R/'fine_grid_temporal.csv',index=False)
 print(pd.DataFrame(summaries).to_string(index=False));print(pd.DataFrame(grid).to_string(index=False))
if __name__=='__main__':main()
