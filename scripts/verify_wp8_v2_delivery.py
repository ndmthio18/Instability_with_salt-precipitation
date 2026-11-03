"""Final provenance and regression gate, with read-only historical checks."""
import json,subprocess,sys
import wp8_v2_regression as v
c=v.c;O=c.OUT
def main():
    checks={}
    for label in ['wp8_zero_mode_20261103','wp8_strict_background_20261103']:
        j=json.loads((c.ROOT/'analysis'/label/'provenance.json').read_text())
        for n,r in j['files'].items():checks['unchanged/'+n]=c.sha(c.ROOT/n)==r['sha256']
    for name in ['transport_wp8_background_v1','transport_wp8_background_v2','transport_wp8_background_v2_probe','transport_wp8_background_v2_probe_fine']:
        d=c.ROOT/'model'/name;m=json.loads((d/'source_manifest.json').read_text())
        for n,r in m.items():checks[name+'/'+n]=c.sha(d/'source'/n)==r['sha256']
    j=json.loads((O/'regression_summary.json').read_text());checks['all_nine_regressions']=len(j['runs'])==9 and j['all_regression_checks_passed']
    eq=json.loads((O/'equation_diagnostics.json').read_text());checks['FD_scale_study']=eq['FD_scale_study_pass']
    comparison=json.loads((O/'late_protocol_comparison.json').read_text());checks.update({k:comparison[k] for k in ['late_bridge_pass','factor16_vs32_fields_exact','probe_does_not_change_accepted_field']})
    p=json.loads((O/'full4T_protocol.json').read_text());checks['no_full_launch']=not p['launched'] and not p['launch_authorized'] and not p['scientific_4T_passed'];checks['original_IC_only']=p['original_initial_state'] and not p['restart'] and not p['replay'];checks['planned_input_hash']=c.sha(O/'pending_full4T.input')==p['input_sha256'];checks['frozen_gates']=p['gates']==json.loads((c.ROOT/'analysis/wp8_zero_mode_20261103/background_protocol.json').read_text())['gates']
    for n in v.CASES:
        d=O/n;plan=json.loads((d/'plan.json').read_text());r=json.loads((d/'archive_receipt.json').read_text());checks[n+'/archive_local']=c.sha(c.CACHE/n/'fields.tar.gz')==r['archive_sha256'];checks[n+'/profiles_local']=c.sha(c.CACHE/n/'post/profiles.npz')==r['profile_sha256'];checks[n+'/auth']=plan['authorization_sha256']==c.sha(c.ROOT/'command/14.md');checks[n+'/validation']=json.loads((d/'validation.json').read_text())['passed']
    for name in ['NUMERICAL_FAILURE_ANALYSIS_ZH.md','SOURCE_REVIEW_AND_CONTRACT_ZH.md','FULL4T_GO_NOGO_AND_COST_ZH.md','SCIENTIFIC_BOUNDARY_ZH.md','README_ZH.md','Figure_Newton_diagnostics.png','Figure_Newton_diagnostics.pdf']:checks['report/'+name]=(O/name).exists()
    # These existing semantic tests run in a separate interpreter, so the bounded
    # campaign's module-local output bindings cannot redirect historical tests.
    test=subprocess.run([sys.executable,'-m','unittest','discover','-s','scripts','-p','test_wp6_paired_stopping.py'],cwd=c.ROOT,capture_output=True,text=True);(O/'stopping_semantics_tests.txt').write_text(test.stdout+test.stderr,encoding='utf-8');checks['stopping_tests_22']=test.returncode==0 and 'Ran 22 tests' in test.stderr
    paper1=subprocess.run([sys.executable,'scripts/prepare_paper1_inputs.py','--verify-only'],cwd=c.ROOT,capture_output=True,text=True);(O/'paper1_readonly_verification.txt').write_text(paper1.stdout+paper1.stderr,encoding='utf-8');checks['paper1_frozen_sources']=paper1.returncode==0
    scripts=[c.ROOT/'scripts'/n for n in ['prepare_wp8_v2.py','wp8_v2_regression.py','validate_wp8_v2.py','report_wp8_v2.py','plot_wp8_v2.py','verify_wp8_v2_delivery.py','wp8_convergence_v2.hh']]
    for path in scripts:
        if path.suffix=='.py':compile(path.read_bytes(),str(path),'exec')
    checks['python_static_compile']=True
    # No long campaign directory or executable invocation is supplied here.
    report=dict(delivery_passed=all(checks.values()),checks=checks,checks_count=len(checks),scientific_full4T_passed=False,full4T_launched=False,diagnostic_replay_not_restart=True)
    c.save(O/'delivery_verification.json',report)
    files=[c.ROOT/'command/14.md']+scripts+[x for x in O.rglob('*') if x.is_file() and x.name!='provenance.json']
    for name in ['transport_wp8_background_v2','transport_wp8_background_v2_probe','transport_wp8_background_v2_probe_fine']:files.extend(x for x in (c.ROOT/'model'/name).rglob('*') if x.is_file())
    c.save(O/'provenance.json',dict(reviewed_commit='87678761f4fa0d64a7205add3f1f0139a1b716e8',files={x.relative_to(c.ROOT).as_posix():dict(sha256=c.sha(x),bytes=x.stat().st_size) for x in files}))
    print(json.dumps(dict(delivery_passed=report['delivery_passed'],checks=len(checks),failures=[k for k,b in checks.items() if not b],full4T_launched=False),indent=2));assert report['delivery_passed']
if __name__=='__main__':main()
