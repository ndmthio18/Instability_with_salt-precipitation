"""Delivery integrity is independent of the failed full-4T scientific gate."""
import hashlib,io,json,unittest
from wp8_strict_campaign import ROOT,OUT,CACHE,sha,save,T
from analyze_wp8 import rows
from test_wp8_strict import TestStrictDiagnostics
from test_wp8 import TestWP8
def main():
    checks={};prior=json.loads((ROOT/'analysis/wp8_zero_mode_20261103/provenance.json').read_text())
    for n,v in prior['files'].items():checks['historical_unchanged/'+n]=sha(ROOT/n)==v['sha256']
    gt=json.loads((OUT/'T/validation.json').read_text());gf=json.loads((OUT/'4T/validation.json').read_text());checks['T_passed']=gt['passed'];checks['full_4T_failure_preserved']=not gf['passed'];checks['accepted_prefix_audited']=gf['accepted_prefix_passed'];checks['rejected_step_not_saved']=gf['accepted_prefix_last_time_s']==501695.8170227767
    plan=json.loads((OUT/'4T/plan.json').read_text());checks['T_gate_frozen_before_4T']=plan['previous_T_gate_sha256']==sha(OUT/'T/validation.json')
    for stage in ['T','4T']:
        d=OUT/stage;p=json.loads((d/'plan.json').read_text());f=json.loads((d/'finished.json').read_text());rec=json.loads((d/'archive_receipt.json').read_text());manifest=json.loads((d/'post/field_manifest.json').read_text())
        checks[stage+'/authorization']=p['authorization_sha256']==sha(ROOT/'command/13.md');checks[stage+'/source_manifest']=p['source_hashes']=={k:v['sha256'] for k,v in json.loads((ROOT/'model/transport_wp8_background_v1/source_manifest.json').read_text()).items()}
        checks[stage+'/ELF']=p['elf_sha256']==f['elf_sha256'];checks[stage+'/input']=p['input_sha256']==sha(d/'params.input');checks[stage+'/worker']=f['worker_sha256']==sha(ROOT/'scripts/wp8_strict_worker.py');checks[stage+'/archive']=sha(CACHE/stage/'fields.tar.gz')==rec['archive_sha256']==f['archive']['sha256'];checks[stage+'/profiles']=sha(CACHE/stage/'post/profiles.npz')==rec['profile_sha256'];checks[stage+'/frames']=rec['verified_frames']==len(manifest)==len(rows(d/'surface_events.csv'));checks[stage+'/no_restart']=not p['restart']
    r=json.loads((OUT/'results.json').read_text());end=r['valid_background_end_s'];checks['no_4T_background_claim']=not r['background_4T_completed'] and end<4*T;checks['no_labels']=r['t_eff'] is None and r['crystallization_delay'] is None;checks['budget_identity']=r['max_budget_identity_error_m']<1e-18
    for n,h in r['source_hashes'].items():checks['comparison_source/'+n]=sha(ROOT/n)==h
    for key,v in r['summary'].items():
        checks[key+'/common_window']=v['shared_states']==4414 and v['last']['time_s']==end
        checks[key+'/own_inventory_window']=all(x['last_time_s']==end for x in v['own_inventory_history'].values())
        checks[key+'/own_inventory_grows']=all(x['minimum_accepted_increment_m']>0 for x in v['own_inventory_history'].values())
    tr=rows(OUT/'strict_transport_final.csv');checks['all_28_regions']=len(tr)==28 and all(float(x['time_s'])==end for x in tr)
    stream=io.StringIO();suite=unittest.TestSuite([unittest.defaultTestLoader.loadTestsFromTestCase(c) for c in [TestWP8,TestStrictDiagnostics]]);tests=unittest.TextTestRunner(stream=stream,verbosity=2).run(suite);checks['numerical_tests']=tests.wasSuccessful();(OUT/'unit_tests.txt').write_text(stream.getvalue(),encoding='utf-8')
    for name in ['EXECUTION_AND_GATES_ZH.md','FAILURE_DIAGNOSIS_ZH.md','STRICT_RESPONSE_ZH.md','COST_AND_PROVENANCE_ZH.md','SCIENTIFIC_REVIEW_ZH.md']:checks['report/'+name]=(OUT/name).exists()
    for stem in ['Figure_A_strict_mean','Figure_B_surface_concentrations','Figure_C_strict_inventory_transport','Figure_D_depth_response']:
        for suffix in ['.png','.pdf']:checks[stem+suffix]=(OUT/(stem+suffix)).exists()
    scripts=[p for p in (ROOT/'scripts').glob('*wp8*strict*.py')]
    for p in scripts:compile(p.read_bytes(),str(p),'exec')
    audit=dict(delivery_integrity_passed=all(checks.values()),scientific_T_passed=True,scientific_4T_passed=False,accepted_prefix_end_s=end,numerical_tests=tests.testsRun,checks=checks)
    save(OUT/'delivery_verification.json',audit)
    files=[ROOT/'command/13.md']+[p for p in OUT.rglob('*') if p.is_file() and p.name!='provenance.json']+scripts
    save(OUT/'provenance.json',dict(reviewed_commit='14e15c7edb32a6a98faf7d915052590d35039ed6',files={p.relative_to(ROOT).as_posix():dict(sha256=sha(p),bytes=p.stat().st_size) for p in files}))
    print(json.dumps(dict(delivery_integrity_passed=audit['delivery_integrity_passed'],checks=len(checks),numerical_tests=tests.testsRun,scientific_4T_passed=False,failures=[k for k,v in checks.items() if not v]),indent=2));assert audit['delivery_integrity_passed']
if __name__=='__main__':main()
