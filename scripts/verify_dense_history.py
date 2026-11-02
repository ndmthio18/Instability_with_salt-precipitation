"""Verify derived-table completeness, frozen inputs, and inventory sign identities."""
import csv,hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'analysis/dense_history_20261102';EV=OUT/'evidence'
def rows(p):
    with p.open(newline='',encoding='utf-8') as f:return list(csv.DictReader(f))
checks={}
plan=json.loads((OUT/'plan.json').read_text())
checks['frozen_plan_matches_remote']=plan==json.loads((EV/'plan.json').read_text())
manifest=json.loads((EV/'product_manifest.json').read_text())
checks['all_collected_product_hashes']=all(hashlib.sha256((EV/name).read_bytes()).hexdigest()==digest for name,digest in manifest.items())
assess=rows(OUT/'run_assessment.csv');inventory=rows(OUT/'inventory_timeseries.csv');times=rows(OUT/'time_reconciliation.csv')
checks['14_runs_10_cases']=len(assess)==14 and len({r['case_id'] for r in assess})==10
checks['1414_time_mappings']=len(times)==1414 and all(r['accepted']=='True' for r in times)
checks['4242_region_frames']=len(inventory)==4242
checks['inputs_and_fields_consistent']=all(r['status'].startswith('field_input_consistent') for r in assess)
checks['sampled_initial_profile_mode_amplitude']=all(r['passed']=='True' for r in rows(OUT/'initial_protocol_check.csv'))
identities=True;monotone_times=True;counts=True
for spec in plan['runs']:
    for region in ['upper_0.9','upper_0.75','upper_0.5']:
        rr=[r for r in inventory if r['id']==spec['id'] and r['region']==region]
        counts &= len(rr)==101
        i0=float(rr[0]['inventory_salt_m'])
        for r in rr:
            i=float(r['inventory_salt_m']);proxy=float(r['inferred_net_downward_m']);relative=float(r['relative_inventory_change'])
            identities &= abs(proxy-(i0-i)) <= 1e-18 and abs(relative-(i/i0-1))<=1e-13
        monotone_times &= all(float(b['time_s'])>float(a['time_s']) for a,b in zip(rr,rr[1:]))
checks.update(region_counts=counts,time_order=monotone_times,inventory_sign_and_normalization=identities,
              analytic_postprocessor_check=json.loads((OUT/'postprocessing_selfcheck.json').read_text())['passed'])
result=dict(passed=all(checks.values()),checks=checks,simulations_started=False,
            scope='postprocessing verification only; not independent discrete-budget or physical validation')
(OUT/'verification.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
raise SystemExit(0 if result['passed'] else 1)
