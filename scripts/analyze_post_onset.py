"""Existing-data amplification metrics; frozen onset classes remain authoritative."""
import csv,hashlib,importlib.util,json,math
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]; OUT=ROOT/'analysis/post_onset_amplification_20261103'
def readcsv(p):return list(csv.DictReader(open(p,encoding='utf-8-sig',newline='')))
def writecsv(p,rows):
 if not rows:return
 keys=list(dict.fromkeys(k for r in rows for k in r))
 with open(p,'w',encoding='utf-8',newline='') as f:
  w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
def num(x,default=np.nan):
 try:return float(x) if x not in (None,'','None') else default
 except (ValueError,TypeError):return default
def clean(v):
 if isinstance(v,dict):return {k:clean(x) for k,x in v.items()}
 if isinstance(v,(list,tuple)):return [clean(x) for x in v]
 if isinstance(v,np.generic):v=v.item()
 if isinstance(v,float) and not math.isfinite(v):return None
 return v
def savejson(p,d):p.write_text(json.dumps(clean(d),ensure_ascii=False,indent=2),encoding='utf-8')
def fit(t,s,lo,hi):
 m=(t>=lo)&(t<=hi)&(s>0)&np.isfinite(s);t=t[m];s=s[m]
 if len(t)<5 or np.ptp(t)<=0:return {'n':len(t),'slope_s_inv':None,'r2':None}
 z=(t-t[0])/np.ptp(t);y=np.log(s);coef=np.polyfit(z,y,1);res=y-np.polyval(coef,z);ss=np.sum((y-y.mean())**2)
 return dict(n=len(t),start_s=t[0],end_s=t[-1],slope_s_inv=coef[0]/np.ptp(t),r2=1-np.sum(res**2)/ss if ss>0 else None,log_change=y[-1]-y[0])
def metrics(a,detect,denom):
 assert a.ndim==2 and a.shape[1]>=3 and np.all(np.isfinite(a[:,:3])) and np.all(np.diff(a[:,0])>0)
 t,mean,s=a[:,:3].T;assert np.all(mean>0) and np.all(s>=0)
 rel=s/mean;result=dict(series_n=len(t),sigma_initial=s[0],sigma_initial_time_s=t[0],sigma_final=s[-1],sigma_final_time_s=t[-1],sigma_literal_min=s.min(),sigma_literal_min_time_s=t[np.argmin(s)],sigma_record_max=s.max(),sigma_record_max_time_s=t[np.argmax(s)],A_rel_end=rel[-1],A_rel_record_max=rel.max(),A_rel_record_max_time_s=t[np.argmax(rel)],top_mean_end=mean[-1],top_mean_record_max=mean.max(),numerical_floor_empirical='unknown',floating_floor_flag=bool(s.min()<=100*np.finfo(float).eps*mean.max()))
 if a.shape[1]>=6:
  result.update(fourier_available=True,fourier_initial=a[0,4],fourier_final=a[-1,4],mode_projection_initial=a[0,5],mode_projection_final=a[-1,5],mode_projection_end_over_initial=a[-1,5]/a[0,5] if a[0,5]!=0 else None,fourier_over_sigma_end=a[-1,4]/s[-1] if s[-1]>0 else None)
  valid=(a[:,4]>0)&(s>0)
  result['fourier_sigma_log_correlation']=float(np.corrcoef(np.log(a[valid,4]),np.log(s[valid]))[0,1]) if valid.sum()>5 and np.std(s[valid])>0 else None
 else:result['fourier_available']=False
 if not np.isfinite(detect) or not np.isfinite(denom) or denom<=0:return result
 post=t>=detect
 if not np.any(post):return result
 last=t[-1];interval=last-detect
 result.update(sigma_detector_min=denom,G_end=s[-1]/denom,G_max=s[post].max()/denom,sigma_post_max=s[post].max(),sigma_post_max_time_s=t[post][np.argmax(s[post])],A_rel_post_max=rel[post].max(),post_window_s=interval,end_over_detect=last/detect if detect>0 else None,post_n=int(post.sum()),post_positive_increment_fraction=float(np.mean(np.diff(s[post])>0)) if post.sum()>1 else None,fall_from_peak_fraction=1-s[-1]/s[post].max())
 for name,lo,hi in [('early_half',detect,detect+.5*interval),('late_half',last-.5*interval,last),('late_quarter',last-.25*interval,last)]:
  f=fit(t,s,lo,hi);result.update({name+'_'+k:v for k,v in f.items()})
 for floor in [1e-14,1e-12,1e-10]:result['G_end_floor_'+str(floor)]=s[-1]/max(denom,floor)
 result['sustained_growth_two_windows']=bool((result.get('late_half_slope_s_inv') or 0)>0 and (result.get('late_quarter_slope_s_inv') or 0)>0 and result.get('late_half_log_change',0)>np.log(1.002) and result.get('late_quarter_log_change',0)>np.log(1.002))
 result['tail_shape']='rise_then_fall' if result['fall_from_peak_fraction']>.002 else ('sustained_growth' if result['sustained_growth_two_windows'] else 'weak_or_flat_tail')
 e=result.get('early_half_slope_s_inv');l=result.get('late_half_slope_s_inv')
 result['late_vs_early_slope']='accelerating' if e is not None and l is not None and l>e else 'decelerating_or_unavailable'
 pre=post&(mean<.1)
 result['mean_threshold_first_s']=t[mean>=.1][0] if np.any(mean>=.1) else None
 result['G_max_mean_below_threshold']=s[pre].max()/denom if np.any(pre) else None
 result['local_threshold_status']='unknown_without_field_check'
 return result

def main():
 spec=importlib.util.spec_from_file_location('frozen_extractor',ROOT/'analysis/pilot0_integrated_20261102/provenance/extract_onset_2025.py');mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod)
 sources=json.loads((OUT/'series_source_manifest.json').read_text());bykey={r['run_key']:r for r in sources}
 inv=json.loads((ROOT/'analysis/phase0_20261102/remote_frodo_inventory.json').read_text())['runs'];remote={r['run_key']:r for r in inv}
 case={r['case_id']:r for r in readcsv(ROOT/'data/external/paper1/158395f3d25a/handoff/cases.csv')}
 rows=[]
 for group,file in [('main','run_audit.csv'),('wp1','wp1_supplement_audit.csv')]:
  for row in readcsv(ROOT/'analysis/onset_semantics_20261102/views'/file):
   r=dict(row,source_group=group);key=r['run_key'];raw=remote.get(key,{});params=raw.get('parameters',{})
   for k,v in case.get(r['case_id'],{}).items():
    if k not in r:r[k]=v
   for out,inp in [('permeability_m2','Permeability'),('evaporation_m_s','EvaporationRate'),('porosity','Porosity'),('initial_amplitude','Amplitude')]:r[out]=params.get('Problem.'+inp,'')
   r['protocol_group']=str(r.get('mode_protocol_hash') or r.get('lsa_observable_protocol_hash') or 'unknown')
   r['nx_over_global_critical']=num(r.get('initialized_mode_nx'))/num(r.get('lsa_global_critical_mode_nx'))
   r['frozen_G_end']=r.get('amplification_since_min','')
   r['series_status']='missing_unverified_host_key' if r.get('host')=='cxhpc' else 'missing'
   r['full_field_review']='available_not_selected' if num(r.get('vtu_count'),0)>=2 else 'unavailable_or_unverified'
   if key in bykey:
    p=ROOT/bykey[key]['local'];b=p.read_bytes();assert hashlib.sha256(b).hexdigest()==bykey[key]['sha256']
    a=np.loadtxt(p,ndmin=2);r.update(series_status='hash_verified',series_sha256=bykey[key]['sha256'])
    ex=mod.extract_onset(str(p));r['recomputed_detection_class']=ex.get('onset_class');r['class_matches_frozen']=r['recomputed_detection_class']==r['dumux_detection_class']
    r['sigma_detector_min_time_s']=ex.get('onset_s');r['recomputed_sigma_min']=ex.get('sigma_min')
    formal=r['dumux_detection_class']=='resolved_onset';detect=num(r.get('dumux_t_detect_s')) if formal else num(r.get('dumux_t_detect_estimate_s'))
    r['growth_reference_status']='formal_resolved' if formal else 'exploratory_candidate_only'
    if not np.isfinite(detect):detect=num(ex.get('onset_s'))
    r['growth_reference_time_s']=detect
    r.update(metrics(a,detect,num(ex.get('sigma_min'))))
    r['input_parameters_sha256']=hashlib.sha256(json.dumps(params,sort_keys=True).encode()).hexdigest()
    r['observed_end_s']=a[-1,0];r['requested_end_s']=num(params.get('TimeLoop.TEnd'))
    r['requested_end_reached']=bool(abs(a[-1,0]-r['requested_end_s'])<max(1e-5,1e-9*a[-1,0]))
    r['returncode']=raw.get('case_manifest.json',{}).get('dumux',{}).get('run',{}).get('returncode')
    r['mean_threshold_warning']=bool(a[:,1].max()>=.1)
   rows.append(clean(r))
 writecsv(OUT/'run_metrics.csv',rows)
 for group in ['main','wp1']:
  rr=[r for r in rows if r['source_group']==group];resolved=[r for r in rr if r['dumux_detection_class']=='resolved_onset' and r['series_status']=='hash_verified']
  print(group,len(rr),'raw',sum(r['series_status']=='hash_verified' for r in rr),'resolved raw',len(resolved))
  for r in sorted(resolved,key=lambda r:r.get('G_end') or 0,reverse=True)[:15]:print(r['case_id'],r['initialized_mode_nx'],round(r.get('G_end',0),4),r.get('A_rel_end'),r.get('top_mean_end'),r['run_key'])
 savejson(OUT/'metric_checks.json',dict(rows=len(rows),class_disagreements=[r['run_key'] for r in rows if r.get('class_matches_frozen') is False],source_manifest_sha256=hashlib.sha256((OUT/'series_source_manifest.json').read_bytes()).hexdigest()))
if __name__=='__main__':main()
