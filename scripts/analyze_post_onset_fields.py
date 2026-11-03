"""Depth and spectral diagnostics from verified existing fields, not new simulations."""
from analyze_post_onset import *
def weights(z):return np.diff(np.r_[z[0],.5*(z[1:]+z[:-1]),z[-1]])
def penetration(z,sd):
 wz=weights(z);energy=sd**2*wz;depth=(z[-1]-z)/(z[-1]-z[0]);total=energy.sum()
 if total<=0:return dict(depth_centroid_H=None,depth90_H=None,lateral_rms=0.)
 order=np.argsort(depth);cs=np.cumsum(energy[order])/total
 return dict(depth_centroid_H=float(np.sum(depth*energy)/total),depth90_H=float(depth[order][min(np.searchsorted(cs,.9),len(z)-1)]),lateral_rms=float(np.sqrt(total/(z[-1]-z[0]))))
def spectrum(x,z,nmode):
 # Integer cosine mode on uniform finite width, omit duplicated endpoint for DFT.
 y=x[:,:-1]-x[:,:-1].mean(axis=1)[:,None];N=y.shape[1];fft=np.fft.rfft(y,axis=1)/N;power=abs(fft)**2;power[:,1:]*=2
 if N%2==0:power[:,-1]*=.5
 e=power.T@weights(z);e[0]=0;total=e.sum();nmode=int(nmode)
 return dict(spectral_fundamental_energy_fraction=float(e[nmode]/total) if nmode<len(e) and total>0 else None,spectral_other_energy_fraction=float((total-e[nmode])/total) if nmode<len(e) and total>0 else None,spectral_second_harmonic_energy_fraction=float(e[2*nmode]/total) if 2*nmode<len(e) and total>0 else None,spectral_dominant_mode=int(np.argmax(e)),spectral_nyquist=N//2,spectral_other_rms=float(np.sqrt(max(0,total-e[nmode])/(z[-1]-z[0]))) if nmode<len(e) else None)
def main():
 rows=readcsv(OUT/'run_metrics.csv');bykey={r['run_key']:r for r in rows};sparse=[];profiles=[];evidence=[]
 selection=json.load(open(OUT/'sparse_field_selection.json'))['runs'];cache=ROOT/'tmp/post_onset_amplification_20261103/fields'
 for i,spec in enumerate(selection):
  info=json.load(open(cache/(str(i)+'.json')));assert info['run_key']==spec['run_key'];p=cache/(str(i)+'.npz');assert hashlib.sha256(p.read_bytes()).hexdigest()==info['npz_sha256'];data=np.load(p);z=data['z_m'];run=bykey[spec['run_key']]
  initial_ok=info['records'][0]['initial_field_error']<1e-12
  for j,rec in enumerate(info['records']):
   tag=str(j>0);x=data['x_'+tag];sd=data['sd_'+tag];mean=data['mean_'+tag];r=dict(run_key=spec['run_key'],case_id=run['case_id'],mode_nx=run['initialized_mode_nx'],**rec,**penetration(z,sd),**spectrum(x,z,num(run['initialized_mode_nx'])))
   r['physical_field_status']='quarantined_constitutive_conflict' if not rec['consistent'] else ('beyond_solubility' if rec['x_max']>=.1 else 'field_input_consistent_binary_unverified')
   r['initial_profile_passed']=initial_ok;sparse.append(r)
   for zi in range(len(z)):profiles.append(dict(run_key=spec['run_key'],frame=rec['frame'],time_s=rec['time_s'],z_m=z[zi],x_mean=mean[zi],sigma_x=sd[zi],sigma_over_mean=sd[zi]/mean[zi]))
  evidence.append({k:v for k,v in info.items() if k not in ['log_tail']})
 writecsv(OUT/'sparse_field_metrics.csv',clean(sparse));writecsv(OUT/'sparse_depth_profiles.csv',clean(profiles));savejson(OUT/'sparse_field_evidence.json',evidence)
 # Reuse audited 1414-frame products; no raw VTU reread.
 base=ROOT/'analysis/dense_history_20261102';manifest=json.load(open(base/'evidence/product_manifest.json'));plans=json.load(open(base/'plan.json'))['runs'];reconciliation=readcsv(base/'time_reconciliation.csv');assess={r['id']:r for r in readcsv(base/'run_assessment.csv')}
 dense=[];spectra=[];summary=[];hashchecks=[]
 for spec in plans:
  rid=spec['id'];d=base/'evidence/runs'/rid;run=bykey[spec['run_key']]
  for name in ['profiles.npz','representative_fields.npz','frame_metrics.csv','inventory_timeseries.csv','audit.json']:
   rel='runs/'+rid+'/'+name;actual=hashlib.sha256((d/name).read_bytes()).hexdigest();assert actual==manifest[rel];hashchecks.append(dict(path=str((d/name).relative_to(ROOT)).replace('\\','/'),sha256=actual))
  pp=np.load(d/'profiles.npz');ff=np.load(d/'representative_fields.npz');z=pp['z_m'];fm=readcsv(d/'frame_metrics.csv');tt=[r for r in reconciliation if r['id']==rid];assert len(tt)==len(pp['time_s'])==101 and all(r['accepted']=='True' for r in tt)
  local=[]
  for j in range(101):
   t=num(tt[j]['matched_diagnostic_time_s']);r=dict(id=rid,run_key=spec['run_key'],case_id=spec['case_id'],mode_nx=run['initialized_mode_nx'],frame=j,time_s=t,**penetration(z,pp['x_sigma'][j]),max_layer_sigma=float(pp['x_sigma'][j].max()),max_layer_relative=float(np.max(pp['x_sigma'][j]/pp['x_mean'][j])),top_sigma=num(fm[j]['top_sigma_node']),top_mean=num(fm[j]['top_x_mean_node']),domain_max=num(fm[j]['domain_x_max']),post_detect=bool(t>=num(run['dumux_t_detect_s'])),G_from_detector=num(fm[j]['top_sigma_node'])/num(run['recomputed_sigma_min']))
   r['threshold_eligible_saved_frame']=bool(r['domain_max']<.1);local.append(r);dense.append(r)
   for q in [.99,.95,.9,.75,.5]:
    idx=np.argmin(abs(z/z[-1]-q));r['sigma_zH_'+str(q)]=pp['x_sigma'][j,idx]
   if j in [0,25,50,75,100]:spectra.append(dict(id=rid,case_id=spec['case_id'],frame=j,time_s=t,threshold_eligible_saved_frame=r['threshold_eligible_saved_frame'],**spectrum(ff['x_'+str(j)],z,num(run['initialized_mode_nx'])),scope='horizontal display subsample; unresolved high modes may alias'))
  eligible=[r for r in local if r['threshold_eligible_saved_frame']];end=eligible[-1];initial=local[0];ss=[r for r in spectra if r['id']==rid and r['threshold_eligible_saved_frame']]
  summary.append(dict(id=rid,run_key=spec['run_key'],case_id=spec['case_id'],mode_nx=run['initialized_mode_nx'],last_pre_observed_threshold_s=end['time_s'],depth_centroid_initial_H=initial['depth_centroid_H'],depth_centroid_last_pre_H=end['depth_centroid_H'],depth90_initial_H=initial['depth90_H'],depth90_last_pre_H=end['depth90_H'],max_relative_pre=max(r['max_layer_relative'] for r in eligible),max_G_pre=max(r['G_from_detector'] for r in eligible if r['post_detect']) if any(r['post_detect'] for r in eligible) else None,max_other_energy_pre_sampled=max(r['spectral_other_energy_fraction'] for r in ss),post_detection_frames_pre=sum(r['post_detect'] for r in eligible),status=assess[rid]['status']))
 # Attach selection and quarantine status to the complete inventory without removing observations.
 sparse_by={}
 for r in sparse:sparse_by.setdefault(r['run_key'],[]).append(r)
 dense_by={r['run_key']:r for r in summary}
 for r in rows:
  rr=sparse_by.get(r['run_key'],[])
  if rr:
   r['full_field_review']='quarantined_constitutive_conflict' if any(q['physical_field_status']=='quarantined_constitutive_conflict' for q in rr) else ('reviewed_initial_only_final_missing' if len(rr)==1 else 'reviewed_initial_final_binary_unverified')
  elif r['run_key'] in dense_by:r['full_field_review']='reviewed_dense_101_frames_binary_unverified'
 writecsv(OUT/'run_metrics.csv',rows)
 writecsv(OUT/'dense_depth_timeseries.csv',clean(dense));writecsv(OUT/'dense_spectral_snapshots.csv',clean(spectra));writecsv(OUT/'dense_summary.csv',clean(summary));savejson(OUT/'dense_reuse_hashes.json',hashchecks)
 print('sparse frames',len(sparse),'dense',len(dense));print('sparse final',[(r['case_id'],r['frame'],r['physical_field_status'],r['spectral_other_energy_fraction'],r['max_layer_relative']) for r in sparse if r['time_s']>0]);print('dense',[(r['id'],r['max_other_energy_pre_sampled'],r['max_relative_pre']) for r in summary])
if __name__=='__main__':main()
