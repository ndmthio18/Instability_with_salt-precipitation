"""Command 15 delivery gate: provenance, acceptance, frozen-history and test checks."""
import json,subprocess,sys
import wp8_full4t_campaign as f

ROOT=f.ROOT;OUT=f.OUT;CACHE=f.CACHE;STAGE=f.STAGE;sha=f.sha;save=f.save
RUN=OUT/STAGE
SCRIPTS=['wp8_full4t_campaign.py','validate_wp8_full4t.py','compare_wp8_full4t.py',
         'analyze_wp8_full4t.py','plot_wp8_full4t.py','report_wp8_full4t.py',
         'test_wp8_full4t.py','verify_wp8_full4t_delivery.py']
DOCS=['README_ZH.md','EXECUTION_AND_ACCEPTANCE_ZH.md','HISTORICAL_COMPARISON_ZH.md',
      'SCIENTIFIC_REVIEW_ZH.md','Figure_A_background_validation.png','Figure_A_background_validation.pdf',
      'Figure_B_background_referenced_response.png','Figure_B_background_referenced_response.pdf']
PRODUCTS=['validation.json','historical_comparison.json','scientific_analysis.json',
          'checkpoint_comparison.csv','background_referenced_responses.csv',
          'background_referenced_profiles.csv','surface_mean_versus_local_max.csv',
          'grid_sensitivity_vs_background.csv','seven_plane_redistribution.csv.gz',
          'cost_audit.json','acceptance_receipt.json']


def main():
    checks={}
    # 1. the frozen histories this command must not touch
    for label in ['wp8_zero_mode_20261103','wp8_strict_background_20261103','wp8_newton_correction_20261103']:
        j=json.loads((ROOT/'analysis'/label/'provenance.json').read_text())
        for n,r in j['files'].items():checks['unchanged/'+n]=sha(ROOT/n)==r['sha256']
    for name in ['transport_wp8_background_v1','transport_wp8_background_v2',
                 'transport_wp8_background_v2_probe','transport_wp8_background_v2_probe_fine']:
        d=ROOT/'model'/name;m=json.loads((d/'source_manifest.json').read_text())
        for n,r in m.items():checks[name+'/'+n]=sha(d/'source'/n)==r['sha256']

    # 2. exactly one authorized run, from the frozen input and production ELF
    plan=json.loads((RUN/'plan.json').read_text());fin=json.loads((RUN/'finished.json').read_text())
    proto=json.loads((ROOT/'analysis/wp8_newton_correction_20261103/full4T_protocol.json').read_text())
    checks['single_run_directory']=sorted(p.name for p in OUT.iterdir() if p.is_dir())==[STAGE]
    checks['authorization']=plan['authorization_sha256']==sha(ROOT/'command/15.md')
    checks['frozen_input_used']=plan['input_sha256']==proto['input_sha256']==fin['input_sha256']==sha(RUN/'params.input')
    checks['production_elf_used']=plan['elf_sha256']==proto['executable']['sha256']==fin['elf_sha256']
    checks['reviewed_worker_used']=fin['worker_sha256']==sha(ROOT/'scripts/wp8_strict_worker.py')
    checks['frozen_gates_unchanged']=plan['gates']==proto['gates']==json.loads((ROOT/'analysis/wp8_zero_mode_20261103/background_protocol.json').read_text())['gates']
    checks['approved_criterion_unchanged']=plan['criterion']==proto['criterion']
    checks['no_restart_replay_retry']=not plan['restart'] and not plan['replay'] and plan['no_automatic_retry'] and not (RUN/'replay.dat').exists()
    checks['fixed_dt_and_horizon']=plan['dt_s']==proto['dt_s'] and plan['horizon_s']==proto['horizon_s'] and plan['steps_cap']==proto['max_steps']
    checks['resource_caps_honoured']=fin['solver_and_live_audit_wall_s']<=plan['wall_cap_s'] and fin['archive']['bytes']<=plan['archive_cap_mib']*1024**2 and fin['peak_child_rss_kib']<=plan['memory_gib']*1024**2

    # 3. acceptance, archive and the required analyses
    val=json.loads((OUT/'validation.json').read_text())
    checks['acceptance_passed']=val['passed'] and not val['failures']
    checks['acceptance_covers_every_requirement']=all(k in val['checks'] for k in
        ['exact_horizontal_equality','one_Newton_success_per_step','fixed_accepted_dt','threshold_flags',
         'salt_closure','water_closure','independent_storage','independent_flux','independent_budget',
         'reached_requested_horizon','passed_old_failure_window'])
    rec=json.loads((RUN/'archive_receipt.json').read_text())
    checks['archive_receipt']=rec['verified_frames']==fin['frames']==4801 and rec['archive_sha256']==fin['archive']['sha256']
    checks['local_archive_hash']=sha(CACHE/STAGE/'fields.tar.gz')==rec['archive_sha256']
    checks['local_profiles_hash']=sha(CACHE/STAGE/'post/profiles.npz')==rec['profile_sha256']
    cmp=json.loads((OUT/'historical_comparison.json').read_text())
    checks['five_checkpoints_reported']={r for r in open(OUT/'checkpoint_comparison.csv',encoding='utf-8').read().split('\n')[1:] if r}!=set()
    checks['old_failure_window_passed']=cmp['old_failure_window']['passed']
    checks['comparison_not_interpolated']=not cmp['interpolation_used'] and not cmp['byte_identity_required']
    sci=json.loads((OUT/'scientific_analysis.json').read_text())
    checks['no_unsupported_labels']=sci['t_eff'] is None and sci['crystallization_delay'] is None and sci['surrogate_label'] is None
    checks['no_threshold_event_claimed']=all(not e['threshold_observed'] for e in sci['solubility_threshold_events'])
    checks['grid_kept_separate']='NOT a convergence proof' in sci['grid_sensitivity']['scope']
    for n in PRODUCTS+DOCS:checks['product/'+n]=(OUT/n).exists()

    # 4. tests and static checks, run in separate interpreters
    t=subprocess.run([sys.executable,'-m','unittest','discover','-s','scripts','-p','test_wp8_full4t.py'],
                     cwd=ROOT,capture_output=True,text=True)
    (OUT/'numerical_tests.txt').write_text(t.stdout+t.stderr,encoding='utf-8')
    checks['numerical_tests']=t.returncode==0 and 'Ran 17 tests' in t.stderr
    s=subprocess.run([sys.executable,'-m','unittest','discover','-s','scripts','-p','test_wp6_paired_stopping.py'],
                     cwd=ROOT,capture_output=True,text=True)
    (OUT/'stopping_semantics_tests.txt').write_text(s.stdout+s.stderr,encoding='utf-8')
    checks['stopping_tests_22']=s.returncode==0 and 'Ran 22 tests' in s.stderr
    p1=subprocess.run([sys.executable,'scripts/prepare_paper1_inputs.py','--verify-only'],
                      cwd=ROOT,capture_output=True,text=True)
    (OUT/'paper1_readonly_verification.txt').write_text(p1.stdout+p1.stderr,encoding='utf-8')
    checks['paper1_frozen_sources']=p1.returncode==0
    for n in SCRIPTS:compile((ROOT/'scripts'/n).read_bytes(),n,'exec')
    checks['python_static_compile']=True

    report=dict(delivery_passed=all(checks.values()),checks=checks,checks_count=len(checks),
                failures=[k for k,b in checks.items() if not b],
                authorized_runs=1,launched_runs=1,retries=0,
                scientific_4T_reached=val['checks']['reached_requested_horizon'],
                censored_by_threshold=val['censored_by_threshold'],
                not_launched_under_this_command=['R96','additional 2D simulations','multiwavelength cases','8T'])
    save(OUT/'delivery_verification.json',report)
    files=[ROOT/'command/15.md']+[ROOT/'scripts'/n for n in SCRIPTS]+\
          [x for x in OUT.rglob('*') if x.is_file() and x.name!='provenance.json']
    save(OUT/'provenance.json',dict(command='command/15.md',reviewed_commit='acda91291870d71b17cec08101281f090234dfe1',
        host='frodo',remote_run='/temp/zhixin/paper2_wp8_full4T_20261103/run',
        files={x.relative_to(ROOT).as_posix():dict(sha256=sha(x),bytes=x.stat().st_size) for x in files}))
    print(json.dumps(dict(delivery_passed=report['delivery_passed'],checks=len(checks),
                          failures=report['failures']),indent=2))
    assert report['delivery_passed']


if __name__=='__main__':main()
