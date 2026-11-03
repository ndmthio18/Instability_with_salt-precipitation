"""Descriptive cohort summaries, one-case representatives and auditable candidate screening."""
from collections import Counter
import re
from analyze_post_onset import *
def stats(rr,key):
 x=np.array([num(r.get(key)) for r in rr]);x=x[np.isfinite(x)]
 return dict(n=len(x),min=float(x.min()),q25=float(np.quantile(x,.25)),median=float(np.median(x)),q75=float(np.quantile(x,.75)),max=float(x.max())) if len(x) else dict(n=0)
def corr(rr,x,y):
 vals=np.array([(num(r.get(x)),num(r.get(y))) for r in rr]);vals=vals[np.all(np.isfinite(vals),axis=1)]
 def rank(a):
  u,inv,count=np.unique(a,return_inverse=True,return_counts=True);r=np.cumsum(count)-.5*(count-1);return r[inv]
 return dict(n=len(vals),spearman=float(np.corrcoef(rank(vals[:,0]),rank(vals[:,1]))[0,1])) if len(vals)>3 and np.std(vals[:,0])>0 and np.std(vals[:,1])>0 else dict(n=len(vals),spearman=None)
def main():
 rows=readcsv(OUT/'run_metrics.csv');main=[r for r in rows if r['source_group']=='main'];supp=[r for r in rows if r['source_group']=='wp1'];summ={};groups=[];cases=[]
 for label,rr in [('main',main),('wp1',supp)]:
  resolved=[r for r in rr if r['dumux_detection_class']=='resolved_onset'];raw=[r for r in resolved if r['series_status']=='hash_verified'];allraw=[r for r in rr if r['series_status']=='hash_verified']
  summ[label]=dict(runs=len(rr),cases=len(set(r['case_id'] for r in rr)),classes=dict(Counter(r['dumux_detection_class'] for r in rr)),incomplete=sum(r.get('record_complete')=='0' for r in rr),raw_count=len(allraw),raw_cases=len(set(r['case_id'] for r in allraw)),resolved_raw=len(raw),G_end_frozen=stats(resolved,'frozen_G_end'),G_end_recomputed=stats(raw,'G_end'),G_max=stats(raw,'G_max'),A_rel_end_all_classes=stats(allraw,'A_rel_end'),A_rel_max_all_classes=stats(allraw,'A_rel_record_max'),sigma_min=stats(raw,'recomputed_sigma_min'),post_window_s=stats(raw,'post_window_s'),end_over_detect=stats(raw,'end_over_detect'),sustained_two_windows=sum(r.get('sustained_growth_two_windows')=='True' for r in raw),rise_then_fall=sum(r.get('tail_shape')=='rise_then_fall' for r in raw),mean_threshold=sum(r.get('mean_threshold_warning')=='True' for r in allraw),fourier_available=sum(r.get('fourier_available')=='True' for r in allraw),fourier_correlation=stats(allraw,'fourier_sigma_log_correlation'),floating_floor=sum(r.get('floating_floor_flag')=='True' for r in allraw),requested_end_reached=sum(r.get('requested_end_reached')=='True' for r in allraw),counts_by_G={str(t):sum(num(r.get('frozen_G_end') if label=='main' else r.get('G_end'))>=t for r in resolved) for t in [1.002,1.1,1.5,2,10,100]})
  for cls in sorted(set(r['dumux_detection_class'] for r in rr)):
   pool=[r for r in rr if r['dumux_detection_class']==cls];groups.append(dict(source_group=label,grouping='detection_class',group=cls,runs=len(pool),cases=len(set(r['case_id'] for r in pool)),raw=sum(r['series_status']=='hash_verified' for r in pool),G_end_median=stats(pool,'G_end').get('median'),frozen_G_end_median=stats(pool,'frozen_G_end').get('median')))
  for field in ['initialized_mode_nx','protocol_group','soil_group','campaign']:
   for v in sorted(set(r.get(field,'unknown') for r in rr)):
    pool=[r for r in rr if r.get(field,'unknown')==v];formal=[r for r in pool if r['dumux_detection_class']=='resolved_onset'];groups.append(dict(source_group=label,grouping=field,group=v,runs=len(pool),cases=len(set(r['case_id'] for r in pool)),resolved=len(formal),G_end_median=stats(formal,'G_end').get('median'),frozen_G_end_median=stats(formal,'frozen_G_end').get('median')))
  for cid in sorted(set(r['case_id'] for r in rr),key=int):
   pool=[r for r in rr if r['case_id']==cid];rep=sorted(pool,key=lambda r:(r.get('is_primary')!='1',r['run_key']))[0];formal=[r for r in pool if r['dumux_detection_class']=='resolved_onset'];cases.append(dict(source_group=label,case_id=cid,runs=len(pool),resolved_runs=len(formal),representative_run_key=rep['run_key'],representative_rule='primary first, then lexicographic run_key; no growth selection',representative_class=rep['dumux_detection_class'],representative_G_end=rep.get('frozen_G_end') or rep.get('G_end'),max_resolved_G_end=max([num(r.get('frozen_G_end') or r.get('G_end')) for r in formal],default=np.nan),max_recorded_A_rel_end=max([num(r.get('A_rel_end')) for r in pool if r.get('A_rel_end')],default=np.nan)))
  complete=[r for r in raw if r.get('record_complete')!='0']
  eligible=[r for r in raw if r.get('full_field_review')!='quarantined_constitutive_conflict']
  summ[label]['resolved_raw_complete']=len(complete)
  summ[label]['G_end_complete_sensitivity']=stats(complete,'G_end')
  summ[label]['known_quarantine_excluded_G_end']=stats(eligible,'G_end')
  summ[label]['late_half_r2']=stats(raw,'late_half_r2')
  summ[label]['late_quarter_r2']=stats(raw,'late_quarter_r2')
  summ[label]['accelerating_tail_count']=sum(r.get('late_vs_early_slope')=='accelerating' for r in raw)
  summ[label]['finite_lsa_sat_raw']=sum(np.isfinite(num(r.get('lsa_t_sat_reference_s'))) for r in allraw)
  representatives=[next(r for r in rr if r['run_key']==c['representative_run_key']) for c in cases if c['source_group']==label and c['representative_class']=='resolved_onset']
  summ[label]['correlations_one_representative_per_case']={x:corr(representatives,x,'G_end') for x in ['end_over_detect','effective_saturation_top','permeability_m2','evaporation_m_s','r_e','nx_over_global_critical']}
  cc=[r for r in cases if r['source_group']==label];summ[label]['case_max_G']=stats(cc,'max_resolved_G_end');summ[label]['case_representative_G']=stats(cc,'representative_G_end');summ[label]['cases_resolved_any']=sum(r['resolved_runs']>0 for r in cc)
  prelsa=[r for r in allraw if np.isfinite(num(r.get('lsa_t_sat_reference_s'))) and num(r['observed_end_s'])<=num(r['lsa_t_sat_reference_s'])];nearlsa=[r for r in prelsa if num(r['observed_end_s'])>=.8*num(r['lsa_t_sat_reference_s'])]
  summ[label]['end_before_lsa_sat']=len(prelsa);summ[label]['end_before_lsa_sat_G_below2']=sum(num(r.get('G_end'))<2 for r in prelsa);summ[label]['end_80to100pct_lsa_sat']=len(nearlsa);summ[label]['end_80to100pct_lsa_sat_G_below2']=sum(num(r.get('G_end'))<2 for r in nearlsa)
  summ[label]['correlations_raw_resolved']={x:corr(raw,x,'G_end') for x in ['post_window_s','end_over_detect','effective_saturation_top','permeability_m2','evaporation_m_s','r_e','nx_over_global_critical']}
 writecsv(OUT/'case_summary.csv',clean(cases));writecsv(OUT/'stratified_summary.csv',clean(groups));savejson(OUT/'population_summary.json',summ)
 # 9 unique physical cases: growth candidates, near-threshold negatives and weak controls.
 choices=[('36877','B','mode comparison / stronger existing amplification'),('25748','B','largest non-36877 growth and lower surface concentration'),('47647','B','different physical tuple with growth and verified final field'),('27710','B','stronger growth but limited solubility headroom'),('31932','C','growth near solubility; no large relative disturbance'),('30858','C','near-threshold negative; final VTU absent'),('45164','D','large relative rank but tiny amplitude; weak-growth comparator'),('34830','D','dense weak-growth control'),('36545','D','dense high-relative-rank control')]
 inv={r['run_key']:r for r in json.load(open(ROOT/'analysis/phase0_20261102/remote_frodo_inventory.json'))['runs']};densekeys={r['run_key']:r['id'] for r in readcsv(OUT/'dense_summary.csv')};sparse=readcsv(OUT/'sparse_field_metrics.csv');cand=[];lsa_rows=[]
 for cid,cat,value in choices:
  rr=[r for r in main if r['case_id']==cid and r['series_status']=='hash_verified'];r=max(([q for q in rr if q['run_key'] in densekeys] if cid in ['34830','36545'] else rr),key=lambda r:num(r.get('G_end'),0));key=r['run_key'];item=inv[key];p=item['parameters'];tail=item.get('log_tail') or '';wall=re.findall(r'Simulation took ([\d.eE+-]+) seconds',tail);wall_s=num(wall[-1]) if wall else np.nan
  found=[x for x in sparse if x['run_key']==key and num(x['time_s'])>0];field_status=found[-1]['physical_field_status'] if found else ('dense_verified' if key in densekeys else 'missing_final_field' if num(r['vtu_count'])==1 else 'available_not_selected')
  rem=num(r.get('lsa_t_sat_reference_s'))-num(r.get('observed_end_s'));rate=num(r.get('late_half_slope_s_inv'));projection=np.log(.01/num(r.get('A_rel_end')))/rate if rate>0 else np.nan
  c={k:r.get(k) for k in ['case_id','run_key','soil_group','climate_group','r_e','beta','gamma','p_hat_b','vg_n','residual_saturation','mualem_l','effective_saturation_top','permeability_m2','evaporation_m_s','porosity','initial_amplitude','initialized_mode_nx','lsa_t_crit_s','lsa_t_seed_theta010_s','lsa_t_initialized_mode_crossing_s','lsa_initialized_mode_source','dumux_detection_class','dumux_t_detect_s','G_end','G_max','A_rel_end','post_window_s','observed_end_s','top_mean_end','late_half_slope_s_inv','late_quarter_slope_s_inv','series_sha256','record_complete','vtu_count']}
  c.update(candidate_category=cat,scientific_value=value,field_status=field_status,lsa_sat_remaining_s=rem,historical_wall_s=wall_s,hypothetical_two_runs_same_window_core_hours=2*wall_s/3600,hypothetical_two_runs_double_window_core_hours=4*wall_s/3600,cost_caveat='historical CPU/wall on one process; host/load/solver differences and long-time stiffness unknown; not authorization',constant_slope_constant_mean_time_to_Arel_1pct_s=projection,nonlinear_evidence='not_established',minimal_future_simulation='matched M0/M1 same mode/profile/amplitude; bounded window selected after common executable and budget gates',lsa_cumulative_growth_status='not_stored; neutral Rc(t) and sharpness are not growth exponents')
  cand.append(c)
  manifest=item.get('case_manifest.json',{});lsa=manifest.get('lsa',{})
  # Snapshot may hold truncated metadata; prefer already cached full manifest.
  if key in densekeys:lsa=json.load(open(ROOT/'analysis/dense_history_20261102/evidence/runs'/densekeys[key]/'case_manifest.json')).get('lsa',{})
  for f in (ROOT/'tmp/post_onset_amplification_20261103/fields').glob('*.json'):
   d=json.load(open(f))
   if d['run_key']==key:lsa=d['lsa']
  lsa_rows.append(dict(case_id=cid,run_key=key,available_fields=';'.join(sorted(lsa)),stored_critical_nx=lsa.get('critical_nx'),stored_crossing_s=lsa.get('t_onset_lsa_s'),stored_stability_margin=lsa.get('stability_margin_log10'),stored_horizon_margin=lsa.get('margin_at_horizon_log10'),scan_points=len(lsa.get('scan_t_hat',[])),finite_growth_eigenvalues_available=False,cumulative_growth_available=False,reason='saved neutral stability scan does not determine temporal growth exponent'))
 writecsv(OUT/'candidate_shortlist.csv',clean(cand));writecsv(OUT/'lsa_growth_availability.csv',clean(lsa_rows));print('population summaries written');print('candidates',[(c['case_id'],c['field_status'],c['historical_wall_s']) for c in cand])
if __name__=='__main__':main()
