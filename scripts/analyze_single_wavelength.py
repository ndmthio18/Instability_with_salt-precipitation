"""Exact-time comparison gates; no interpolation or relative error for unresolved pairs."""
import argparse,importlib.util,json,re,hashlib
from pathlib import Path
import numpy as np
from analyze_post_onset import readcsv,writecsv,savejson,num
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'analysis/single_wavelength_validation_20261103';CACHE=ROOT/'tmp/single_wavelength_validation_20261103';BASE=ROOT/'analysis/pilot0_integrated_20261102/sal_migration/runs';PRE=json.load(open(OUT/'preregistration.json'));TH=PRE['thresholds']
def match(a,b):
 ta=np.array([num(r['time_s']) for r in a]);tb=np.array([num(r['time_s']) for r in b]);pairs=[]
 for i,t in enumerate(ta):
  if t>tb.max()+1e-7:continue
  j=int(np.argmin(abs(tb-t)));assert abs(tb[j]-t)<1e-7,('time mismatch',t,tb[j]);pairs.append((i,j))
 return pairs
def valdiff(a,b,key):
 ii=match(a,b);av=np.array([num(a[i][key]) for i,j in ii]);bv=np.array([num(b[j][key]) for i,j in ii]);return av,bv

def loop_time(log,n):
 pattern=r'Time step\s+'+str(n)+r'\s+done in [\d.eE+-]+ seconds\. Wall clock time: ([\d.eE+-]+)'
 m=re.search(pattern,log);return float(m.group(1)) if m else None

def main(stage):
 names=[p.name for p in (OUT/'runs').iterdir() if p.name.endswith('-'+stage)];results={};diffs=[];performance=[];pairrows=[];source=ROOT/'analysis/pilot0_integrated_20261102/provenance/extract_onset_2025.py';spec=importlib.util.spec_from_file_location('frozen',source);ext=importlib.util.module_from_spec(spec);spec.loader.exec_module(ext)
 for name in sorted(names):
  d=OUT/'runs'/name;role=name.split('-')[1];p=json.load(open(d/'plan.json'));finished=json.load(open(d/'finished.json'));audit=json.load(open(d/'post/field_audit.json'));model=name.split('-')[0];w=readcsv(d/'post/weighted_diagnostics.csv');ref=readcsv(OUT/'reference'/role/'weighted_diagnostics.csv');checks={};metrics={}
  checks['finished_horizon']=finished['solver']['returncode']==0 and finished['solver']['stop']['reason']=='requested_horizon' and abs(num(w[-1]['time_s'])-p['T_s'])<1e-7
  checks['fields_match_input']=audit['all_rho_pass'] and audit['all_sw_pass'] and audit['initial_field_error']<1e-12
  for key,tol in [('weighted_top_mean',TH['mean_and_storage_relative']),('whole_salt_m',TH['mean_and_storage_relative']),('whole_water_m',TH['mean_and_storage_relative']),('top_max',TH['local_max_relative'])]:
   a,b=valdiff(ref,w,key);err=abs(b-a);metrics[key+'_max_relative']=float(np.max(err/abs(a)));checks[key]=metrics[key+'_max_relative']<=tol
   for (i,j),delta in zip(match(ref,w),b-a):diffs.append(dict(model=model,role=role,stage=stage,time_s=num(ref[i]['time_s']),metric=key,reference=num(ref[i][key]),reduced=num(w[j][key]),difference=delta))
  for key in ['weighted_top_sigma','weighted_fourier','weighted_projection']:
   a,b=valdiff(ref,w,key)
   if role=='M1':metrics[key+'_max_relative']=float(np.max(abs(b-a)/abs(a)));checks[key]=metrics[key+'_max_relative']<=TH['weighted_mode_relative_full']
   else:metrics[key+'_max_absolute']=float(np.max(abs(b)));checks[key]=metrics[key+'_max_absolute']<1e-11
   for i,j in match(ref,w):diffs.append(dict(model=model,role=role,stage=stage,time_s=num(ref[i]['time_s']),metric=key,reference=num(ref[i][key]),reduced=num(w[j][key]),difference=num(w[j][key])-num(ref[i][key])))
  a=np.load(CACHE/'reference'/role/'post/profiles.npz');b=np.load(CACHE/'runs'/name/'post/profiles.npz');assert np.allclose(a['z_m'],b['z_m'],rtol=0,atol=1e-14);z=a['z_m'];dz=np.diff(np.r_[z[0],(z[1:]+z[:-1])/2,z[-1]])
  meanerr=[];sderr=[]
  for i,j in match(ref,w):
   meanerr.append(np.max(abs(b['x_mean'][j]-a['x_mean'][i])/a['x_mean'][i]));den=np.sum(a['x_sigma'][i]**2*dz)
   if role=='M1':sderr.append(np.sqrt(np.sum((b['x_sigma'][j]-a['x_sigma'][i])**2*dz)/den))
  metrics['profile_mean_max_relative']=float(max(meanerr));checks['profile_mean']=max(meanerr)<TH['profile_mean_relative']
  if sderr:metrics['profile_sigma_max_relative_L2']=float(max(sderr));checks['profile_sigma']=max(sderr)<TH['depth_sigma_L2_relative']
  budget=readcsv(d/'transport_budget.csv');refb=readcsv(BASE/role/'transport_budget.csv');metrics['max_salt_residual_m']=max(abs(num(r['salt_residual_m'])) for r in budget);metrics['max_water_residual_m']=max(abs(num(r['water_residual_m'])) for r in budget);checks['salt_closure']=metrics['max_salt_residual_m']<TH['salt_residual_absolute_m'];checks['water_closure']=metrics['max_water_residual_m']<TH['water_residual_absolute_m']
  checks['whole_budget']=max(abs(num(r['whole_salt_residual_m'])) for r in budget)<TH['salt_residual_absolute_m'] and max(abs(num(r['whole_water_residual_m'])) for r in budget)<TH['water_residual_absolute_m']
  for region in ['upper_0.900000','upper_0.750000','upper_0.500000']:
   rr=[r for r in refb if r['region']==region];bb=[r for r in budget if r['region']==region]
   assert abs(num(rr[0]['actual_z_over_H'])-num(bb[0]['actual_z_over_H']))<1e-13
   for key,atol,rtol in [('inventory_salt_m',0,TH['mean_and_storage_relative']),('cut_salt_total_down_m_s',TH['cut_flux_absolute_m_s'],TH['cut_flux_relative']),('cumulative_cut_salt_m',TH['cumulative_cut_absolute_m'],TH['cumulative_cut_relative'])]:
    av,bv=valdiff(rr,bb,key);err=abs(bv-av);metrics[region+'/'+key+'_max_abs']=float(err.max());checks[region+'/'+key]=bool(np.all(err<=atol+rtol*abs(av)))
  metrics['max_endpoint_difference']=max(num(r['max_endpoint_difference']) for r in w);checks['endpoint_symmetry']=metrics['max_endpoint_difference']<1e-6*TH['endpoint_difference_relative_to_seed']
  odd=[num(r['odd_neumann_energy_fraction']) for r in w if np.isfinite(num(r['odd_neumann_energy_fraction']))];metrics['max_odd_neumann_energy_fraction']=max(odd,default=None);checks['odd_neumann_modes']=not odd or max(odd)<TH['odd_neumann_mode_energy_fraction']
  checks['below_solubility']=max(num(r['domain_max']) for r in w)<.1
  onset=ext.extract_onset(str(next(d.glob('*sigmatop.csv')))) if role=='M1' else None
  if role=='M1' and stage=='full':
   refon=ext.extract_onset(str(next((BASE/role).glob('*sigmatop.csv'))));ot=onset.get('onset_refined_s');rt=refon.get('onset_refined_s');metrics['onset_relative']=abs(ot-rt)/rt if ot and rt else None;checks['onset']=onset.get('onset_class')=='resolved_onset' and metrics['onset_relative'] is not None and metrics['onset_relative']<TH['onset_relative_full_window']
  results[name]=dict(passed=all(checks.values()),checks=checks,metrics=metrics,onset=onset,frames=audit['frames'])
  log=(d/'solver.log').read_text();performance.append(dict(run=name,host='frodo',stage=stage,cells_x=p['nx_cells'],vertices=p['DOF_vertices'],scalar_unknowns=p['scalar_unknowns'],accepted_steps=finished['solver']['stop']['accepted_steps'],wall_s=finished['solver']['wall_s'],cpu_s=finished['solver']['cpu_s'],peak_rss_kib=finished['solver']['peak_rss_kib'],newton_iterations=finished['solver']['newton_iterations'],linear_iterations=None,linear_iterations_status='not printed by unchanged ELF',logged_failure_retry_lines=sum(not line.lstrip().startswith('--') for line in finished['solver']['failure_or_retry_log_lines']),output_bytes=finished['solver']['solver_output_bytes'],loop_first10_s=loop_time(log,10),loop_total_s=loop_time(log,finished['solver']['stop']['accepted_steps']),post_and_archive_s=finished['post_wall_s']))
 # Pair initial correction and absolute resolution assessment, never divide by near-zero Delta C.
 models=sorted(set(n.split('-')[0] for n in names));pairchecks={}
 for model in models:
  pair={role:readcsv(OUT/'runs'/f'{model}-{role}-{stage}'/'transport_budget.csv') for role in ['M0','M1']}
  for region in ['upper_0.900000','upper_0.750000','upper_0.500000']:
   a=[r for r in pair['M0'] if r['region']==region];b=[r for r in pair['M1'] if r['region']==region];init=num(b[0]['inventory_salt_m'])-num(a[0]['inventory_salt_m']);pairchecks[model+'/'+region+'/initial_inventory']=abs(init)<TH['initial_inventory_difference_m']
   ref0=[r for r in readcsv(BASE/'M0/transport_budget.csv') if r['region']==region];ref1=[r for r in readcsv(BASE/'M1/transport_budget.csv') if r['region']==region]
   for i,j in match(a,b):
    tt=num(a[i]['time_s']);fi=int(np.argmin([abs(num(r['time_s'])-tt) for r in ref0]));assert abs(num(ref0[fi]['time_s'])-tt)<1e-7
    dc=num(b[j]['cumulative_cut_salt_m'])-num(a[i]['cumulative_cut_salt_m']);di=num(b[j]['inventory_salt_m'])-num(a[i]['inventory_salt_m'])-init;dr=num(b[j]['salt_residual_m'])-num(a[i]['salt_residual_m']);reference_dc=num(ref1[fi]['cumulative_cut_salt_m'])-num(ref0[fi]['cumulative_cut_salt_m']);env=abs(num(b[j]['salt_residual_m']))+abs(num(a[i]['salt_residual_m']));floor=max(env,TH['salt_residual_absolute_m']);pairrows.append(dict(model=model,stage=stage,region=region,time_s=tt,initial_inventory_difference_m=init,corrected_inventory_difference_m=di,delta_C_m=dc,reference_delta_C_m=reference_dc,pair_difference_from_full_m=dc-reference_dc,residual_envelope_m=env,engineering_absolute_floor_m=floor,identity_error_m=di+dc-dr,interpretation='unresolved' if abs(dc)<floor else 'requires_grid_and_time_error_assessment'))
 refinement={}
 if stage=='short':
  for role in ['M0','M1']:
   a=readcsv(OUT/'runs'/f'R24-{role}-short/post/weighted_diagnostics.csv');b=readcsv(OUT/'runs'/f'R48-{role}-short/post/weighted_diagnostics.csv')
   for key in ['weighted_top_mean','whole_salt_m','weighted_top_sigma','weighted_fourier','weighted_projection']:
    if role=='M0' and key in ['weighted_top_sigma','weighted_fourier','weighted_projection']:continue
    av,bv=valdiff(a,b,key);err=float(np.max(abs(bv-av)/abs(bv)));tol=TH['mean_and_storage_relative'] if key in ['weighted_top_mean','whole_salt_m'] else TH['weighted_mode_relative_R24_R48'];refinement[role+'/'+key]=dict(max_relative=err,tolerance=tol,passed=err<tol)
 timing=[]
 for role in ['M0','M1']:
  d=OUT/'runs'/f'F600-{role}-bench'
  if (d/'finished.json').exists():
   full=loop_time((d/'solver.log').read_text(),10)
   for n in [24,48]:
    rr=OUT/'runs'/f'R{n}-{role}-short';small=loop_time((rr/'solver.log').read_text(),10);timing.append(dict(role=role,model='R'+str(n),full_10step_loop_s=full,reduced_10step_loop_s=small,speedup=full/small,host='frodo',scope='same first10 accepted steps; VTU every step; excludes startup; not full-window M1 speedup'))
 numerical=all(r['passed'] for r in results.values()) and all(pairchecks.values()) and all(r['passed'] for r in refinement.values());speed=len(timing)==4 and min(r['speedup'] for r in timing)>=TH['speedup_minimum']
 summary=dict(stage=stage,runs=results,pair_checks=pairchecks,refinement=refinement,numerical_passed=numerical,same_hardware_timing=timing,speed_gate_passed=speed,speed_gate_status=('passed' if speed else 'pending_timing' if len(timing)<4 else 'failed'),advance_to_full=(numerical and speed) if stage=='short' else None,transport_effectiveness='unresolved; U_dt/U_z unknown; no t_eff defined',baseline_truth='fixed-grid numerical reference, not converged continuum truth')
 savejson(OUT/(stage+'_validation.json'),summary);writecsv(OUT/(stage+'_differences.csv'),diffs);writecsv(OUT/(stage+'_paired_transport.csv'),pairrows);writecsv(OUT/(stage+'_performance.csv'),performance)
 print(json.dumps(dict(stage=stage,numerical_passed=numerical,speed_gate_passed=speed,failed={n:[k for k,v in r['checks'].items() if not v] for n,r in results.items()},refinement=refinement),indent=2))
if __name__=='__main__':p=argparse.ArgumentParser();p.add_argument('stage',choices=['short','full']);a=p.parse_args();main(a.stage)
