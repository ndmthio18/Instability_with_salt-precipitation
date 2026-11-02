"""Reproduce Phase-0 tables from frozen inputs and the read-only remote snapshot.
No solver imports, no simulation, no network. Run prepare_paper1_inputs --verify-only first.
"""
import csv
import hashlib
import json
import math
from collections import Counter
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
PARENT=ROOT.parents[1]
OUT=ROOT/'analysis/phase0_20261102'
DATA=ROOT/'data/external/paper1/158395f3d25a/handoff'
def read(name):
    with (DATA/name).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def number(x):
    try:
        v=float(x)
        return v if math.isfinite(v) else None
    except (ValueError,TypeError):return None
def save(name,rows):
    if not rows:raise ValueError(name+' is empty')
    keys=list(dict.fromkeys(k for row in rows for k in row))
    with (OUT/name).open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader();w.writerows(rows)
def digest(p):return hashlib.sha256(p.read_bytes()).hexdigest()

lock=json.loads((ROOT/'paper1_sources.lock.json').read_text())
for name,expected in lock['release']['handoff_files'].items():
    assert digest(DATA/name)==expected,name
runs=read('dumux_onsets.csv')
cases={r['case_id']:r for r in read('cases.csv')}
obs={r['case_id']:r for r in read('observations.csv')}
seeds={r['case_id']:r for r in read('seeded_modes.csv') if number(r['sharpness_min'])==.1}
snapshot=json.loads((OUT/'remote_frodo_inventory.json').read_text())
remote={r['run_key']:r for r in snapshot['runs']}
field_evidence=json.loads((OUT/'candidate_field_evidence.json').read_text())
fields={r['run_key']:r for r in field_evidence}
assert len(remote)==len(snapshot['runs'])
def availability(key):
    s=remote.get(key,{})
    p=s.get('parameters',{})
    times=s.get('pvd_times_s',[])
    arrays=sorted({v for a in s.get('vtu_header_arrays',{}).values() for v in a})
    files=s.get('files',[])
    onset=s.get('case_manifest.json',{}).get('dumux',{}).get('onset',{})
    return dict(remote_path=s.get('path',''),remote_access=('found' if s.get('exists') else 'path_missing') if s else 'unverified_host_key',
        vtu_count=s.get('vtu_count'),pvd_times_s=json.dumps(times),visible_array_names_in_first_64k=';'.join(arrays),
        array_inspection_scope='prefix_only_not_complete',
        input_available=any(f.endswith('.input') for f in files) if s.get('exists') else None,
        initial_profile_available=any(f.endswith('_eigenprofile.dat') for f in files) if s.get('exists') else None,
        manifest_available='case_manifest.json' in files if s.get('exists') else None,
        vtk_interval=p.get('TimeLoop.OutputInterval'),krw_exponent=p.get('Problem.KrwMualemExponent'),
        perturbation_type=p.get('Problem.PerturbationType'),
        diagnostics_columns=';'.join(s.get('sigmatop_header',[])),
        recorded_top_mean_max=(s.get('sigmatop_column_max') or [None,None])[1],
        first_recorded_top_mean_ge_0p1_s=s.get('first_top_mean_ge_0p1_s'),
        top_local_max_timeseries='not_recorded' if s.get('sigmatop_header') else 'unverified',
        conservative_flux_timeseries='not_found_in_inventory' if s.get('exists') else 'unverified',
        remote_sigmatop_sha256=s.get('sigmatop_sha256'),
        remote_returncode=s.get('case_manifest.json',{}).get('dumux',{}).get('run',{}).get('returncode'),
        remote_final_time_s=(s.get('sigmatop_last') or [None])[0],
        remote_manifest_onset_class=onset.get('onset_class'),
        transport_reuse=('field_diagnostic_candidate_no_conservative_flux' if s.get('vtu_count',0)>2 else 'onset_and_sparse_field_only') if s.get('exists') else 'onset_table_only',
        sampled_field_density_status=('conflicts_with_current_input' if fields[key]['density_check']['max_abs_error_explicit_input_law_kg_m3']>1e-6 else 'sample_consistent') if key in fields else 'not_sampled')

audit=[]
for r in runs:
    cid=r['case_id'];c=cases[cid];o=obs[cid];seed=seeds.get(cid,{})
    a=dict(r);a['scope']='frozen_handoff_20260913'
    a.update({k:c[k] for k in ['case_key','soil_id','forcing_id','soil_group','climate_group','config_key','effective_saturation_top','reference_time_s','application_horizon_s','r_e','beta','gamma','p_hat_b','vg_n','residual_saturation','mualem_l']})
    a.update(lsa_t_crit_s=o['t_onset_s'],lsa_t_seed_010_s=seed.get('t_onset_seeded_s'),lsa_t_sat_reference_s=o['t_precipitation_s'],
             lsa_sat_status=o['precipitation_status'],lsa_quality=o['quality_status'],lsa_regression_eligible=o['onset_regression_eligible'],
             constitutive_branch=o['constitutive_branch'],mode_protocol_hash=seed.get('seeding_protocol_hash'))
    a.update(availability(r['run_key']))
    last=number(r['t_last_s']);det=number(r['onset_s']);sat=number(o['t_precipitation_s']);hor=number(c['application_horizon_s'])
    a['post_detect_record_s']=last-det if last is not None and det is not None else None
    a['record_exceeds_lsa_sat_reference']=last>sat if last is not None and sat is not None else None
    a['record_exceeds_application_horizon']=last>hor*(1+1e-8) if last is not None and hor else None
    a['stop_interpretation']='requested_mode_multiple_end_reached' if r['record_complete']=='1' else 'incomplete_reason_requires_log_review'
    a['physical_sat_stop_verified']=False
    a['remote_record_end_matches_frozen']=math.isclose(a['remote_final_time_s'],last,rel_tol=1e-8,abs_tol=1e-5) if a['remote_final_time_s'] is not None and last is not None else None
    a['post_detect_field_snapshots']=sum(t>det for t in json.loads(a['pvd_times_s'])) if det else None
    audit.append(a)
save('run_audit.csv',audit)

wp=[]
for p in sorted((PARENT/'paper/1/analysis/wp1_paired_v2_raw/runs').glob('*/case_manifest.json')):
    m=json.loads(p.read_text());d=m['dumux'];o=d['onset'];n=m['dumux_numerics']['physical'];l=m['lsa']
    key='frodo:paper1_wp1_paired_v2/runs/'+p.parent.name
    row=dict(run_key=key,case_id=m['case']['case_id'],scope='wp1_paired_v2_separate_supplement',
             local_manifest=str(p),local_manifest_sha256=digest(p),onset_class=o.get('onset_class'),onset_s=o.get('onset_refined_s') or o.get('onset_s'),
             returncode=d['run'].get('returncode'),seeded_nx=l.get('requested_nx'),t_lsa_mode_s=l.get('t_onset_lsa_s'),t_end_s=n.get('t_end_s'),t_last_s=o.get('t_last_s'),
             onset_low_s=o.get('onset_low_s'),onset_high_s=o.get('onset_high_s'),
             onset_time_interpretation='resolved_point_with_interval' if o.get('onset_class')=='resolved_onset' else 'weak_interval_not_validated_point',
             lsa_sat_reference_s=l.get('t_precipitation_s'),lsa_sat_status=l.get('precipitation_status'),
             executable_sha256_16=d['run'].get('executable_sha256_16'))
    row.update(availability(key));wp.append(row)
save('wp1_supplement_audit.csv',wp)

# Candidate rankings are prospective screening evidence, not transport labels.
pool=[r for r in audit if r['is_primary']=='1' and r['realizes_seed_010']=='1' and r['lsa_regression_eligible']=='1']
rank=[]
for r in pool:
    ts=number(r['lsa_t_seed_010_s']);tp=number(r['lsa_t_sat_reference_s']);hor=number(r['application_horizon_s'])
    if not ts:continue
    row=dict(r)
    row['sat_over_seed']=tp/ts if tp else None
    row['available_reference_window_s']=min(tp,hor)-ts if tp else hor-ts
    row['boundary_score']=abs(math.log(tp/ts)) if tp else None
    rank.append(row)
save('pilot_screening.csv',sorted(rank,key=lambda r:(r['boundary_score'] is None,r['boundary_score'] or 0)))
chosen=[]
for label,cid,key in [('A','36877','frodo:c1b/runs_600/shard-001/0017--case-36877--native'),('B','38852','frodo:paper1_wp1_paired_v2/runs/case-38852--both--nx-16')]:
    chosen.append({'label':label,'case_id':int(cid),'status':'candidate_not_transport_validated','case':cases[cid],
        'lsa_observation':obs[cid],'seed_010':seeds[cid],'preferred_existing_run_key':key,
        'actual_input_parameters':remote[key]['parameters'],'existing_source_hashes':remote[key]['small_file_hashes'],
        'field_sample':fields[key]})
(OUT/'pilot_candidates.json').write_text(json.dumps(chosen,ensure_ascii=False,indent=2),encoding='utf-8')
srcdir=OUT/'source_snapshot';srcdir.mkdir(exist_ok=True)
for path,source in snapshot['sources'].items():
    raw=source['text'].encode()
    (srcdir/Path(path).name).write_bytes(raw)
(srcdir/'source_identity.json').write_text(json.dumps({path:{'remote_original_bytes_sha256':s['sha256'],
    'local_newline_normalized_sha256':hashlib.sha256(s['text'].encode()).hexdigest(),
    'normalization':'Python read_text universal newline conversion; source snapshots are review-only, not an executable provenance claim'}
    for path,s in snapshot['sources'].items()},indent=2),encoding='utf-8')

summary=dict(frozen_runs=len(audit),frozen_cases=len({r['case_id'] for r in audit}),
    primary=sum(r['is_primary']=='1' for r in audit),nonprimary=sum(r['is_primary']=='0' for r in audit),
    resolved_primary=sum(r['is_primary']=='1' and r['onset_class']=='resolved_onset' for r in audit),
    frozen_onset_classes=dict(Counter(r['onset_class'] for r in audit)),
    incomplete=sum(r['record_complete']=='0' for r in audit),
    remote_access=dict(Counter(r['remote_access'] for r in audit)),
    frozen_vtu_counts=dict(Counter(str(r['vtu_count']) for r in audit)),
    modes_observation_window_multiples=dict(Counter(round(float(r['window_over_lsa'])) for r in audit)),
    exceeds_lsa_sat_reference=sum(r['record_exceeds_lsa_sat_reference'] is True for r in audit),
    finite_lsa_sat_reference=sum(bool(r['lsa_t_sat_reference_s']) for r in audit),
    exceeds_application_horizon=sum(r['record_exceeds_application_horizon'] is True for r in audit),
    top_mean_threshold_observed=sum(r['first_recorded_top_mean_ge_0p1_s'] is not None for r in audit),
    remote_record_end_mismatch=sum(r['remote_record_end_matches_frozen'] is False for r in audit),
    sampled_field_density_conflicts=[r['run_key'] for r in field_evidence if r['density_check']['max_abs_error_explicit_input_law_kg_m3']>1e-6],
    wp1_runs=len(wp),wp1_cases=len({r['case_id'] for r in wp}),wp1_onset_classes=dict(Counter(r['onset_class'] for r in wp)),
    wp1_vtu_counts=dict(Counter(str(r['vtu_count']) for r in wp)),
    wp1_case_overlap_with_frozen=len({str(r['case_id']) for r in wp}&{r['case_id'] for r in audit}),
    simulations_started=False,remote_coverage_limit='cxhpc 116 rows: host key verification failed; no key acceptance or bypass',
    source_identity_limit='remote source is inspected current checkout, not proof every historical executable matches')
assert summary['frozen_runs']==436 and summary['frozen_cases']==253 and summary['resolved_primary']==293
assert len({r['run_key'] for r in audit})==436
assert len(wp)==84
(OUT/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2),encoding='utf-8')
sources=[DATA/n for n in ['dumux_onsets.csv','cases.csv','observations.csv','seeded_modes.csv']]
sources += [PARENT/'paper/1/analysis_manifest.json',PARENT/'paper/1/protocol_contract.md',OUT/'remote_frodo_inventory.json',OUT/'candidate_field_evidence.json']
(OUT/'audit_manifest.json').write_text(json.dumps({'input_sha256':{str(p):digest(p) for p in sources},
    'script_sha256':{p.name:digest(p) for p in [Path(__file__),ROOT/'scripts/collect_phase0_remote.py',ROOT/'scripts/inspect_phase0_candidates.py']},
    'output_sha256':{p.name:digest(p) for p in list(OUT.glob('*.csv'))+[OUT/'pilot_candidates.json',OUT/'summary.json']},'source_lock_unchanged':True,'new_simulations':False},ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps(summary,indent=2))
