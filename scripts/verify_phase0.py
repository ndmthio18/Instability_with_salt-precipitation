"""Check audit accounting, candidate provenance and saved artifact hashes."""
import csv
import hashlib
import json
from datetime import datetime, timezone
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'analysis/phase0_20261102'
def rows(name):
    with (OUT/name).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
r=rows('run_audit.csv');w=rows('wp1_supplement_audit.csv')
c=json.loads((OUT/'pilot_candidates.json').read_text())
m=json.loads((OUT/'audit_manifest.json').read_text())
checks={}
checks['frozen_accounting']=len(r)==436 and len({x['run_key'] for x in r})==436 and len({x['case_id'] for x in r})==253
checks['resolved_primary_293']=sum(x['is_primary']=='1' and x['onset_class']=='resolved_onset' for x in r)==293
checks['supplement_84_runs_48_cases']=len(w)==84 and len({x['case_id'] for x in w})==48
checks['remote_coverage']=sum(x['remote_access']=='found' for x in r)==320 and all(x['remote_access']=='found' for x in w)
checks['incomplete_retained']=sum(x['record_complete']=='0' for x in r)==12
checks['two_named_candidates']=[x['case_id'] for x in c]==[36877,38852]
checks['candidate_sample_density_and_validity']=all(x['field_sample']['density_check']['max_abs_error_explicit_input_law_kg_m3']<1e-6 and x['field_sample']['final_top_max']<.1 for x in c)
checks['B_weak_class_preserved']=next(x for x in w if 'case-38852--' in x['run_key'])['onset_class']=='weak_interval_onset'
checks['output_hashes']=all(hashlib.sha256((OUT/n).read_bytes()).hexdigest()==h for n,h in m['output_sha256'].items())
checks['source_hashes']=all(hashlib.sha256(Path(n).read_bytes()).hexdigest()==h for n,h in m['input_sha256'].items())
checks['no_simulation_in_audit_manifest']=m['new_simulations'] is False
report={'checked_utc':datetime.now(timezone.utc).isoformat(),'checks':checks,
        'passed':all(checks.values()),'scope':'audit artifacts only; no transport or solver validation claimed'}
(OUT/'verification.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps(report,indent=2))
if not report['passed']:raise SystemExit(1)
