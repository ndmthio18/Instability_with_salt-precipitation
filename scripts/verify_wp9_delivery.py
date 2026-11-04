"""WP9 delivery gate: provenance, frozen histories, products and tests."""
import json,subprocess,sys
from pathlib import Path
import wp9_common as c

SCRIPTS=['wp9_common.py','wp9_lsa_audit.py','wp9_lsa_vs_dns.py','wp9_corrections.py',
         'wp9_r96_protocol.py','wp9_r96_campaign.py','wp9_pair_frames.py','wp9_compare.py',
         'wp9_validate_r96.py','wp9_post.py','wp9_decision.py','wp9_plot.py',
         'wp9_hpcc_portability.py','wp9_hpcc_job.sbatch','wp9_hpcc_threads.sbatch',
         'wp9_hpcc_worker.py','wp9_hpcc_m0_job.sbatch','wp9_r96_m0_campaign.py',
         'test_wp9.py',
         'verify_wp9_delivery.py']
DOCS=['README_ZH.md','CORRECTION_v2_ZH.md','LSA_WAVENUMBER_AUDIT_ZH.md',
      'R24_R48_R96_COMPARISON_ZH.md','DECISION_ZH.md','HPCC_PORTABILITY_ZH.md',
      'Figure_A_lsa_wavenumber_audit.png','Figure_A_lsa_wavenumber_audit.pdf',
      'Figure_B_resolution_comparison.png','Figure_B_resolution_comparison.pdf',
      'Figure_C_background_and_transport.png','Figure_C_background_and_transport.pdf',
      'Figure_D_corrections.png','Figure_D_corrections.pdf']
PRODUCTS=['corrections_v2.json','decision_criteria.json','amendment_20261104_E01.json',
          'wp9_post_equivalence.json','frame_adapter_selftest.json',
          'execution_policy_20261104.json','background_inventory_attribution.csv','solubility_censoring.csv',
          'lsa_audit.json','lsa_wavenumber_scan.csv','lsa_box_modes.csv','lsa_mode_amplification.csv',
          'lsa_versus_dns.json','lsa_versus_dns.csv',
          'r96_protocol.json','r96_admission.json','user_authorization.json','stage.json','job.json',
          'resolution_series.csv','resolution_spectra.csv','resolution_checkpoint_comparison.csv',
          'resolution_seven_plane.csv','cascade_test.csv','resolution_comparison.json','decision.json']
R96_PRODUCTS=['r96_validation.json','archive_receipt.json','reduction_receipt.json',
              'frame_pairing_receipt.json']
M0_PRODUCTS=['r96_m0/protocol.json','r96_m0/admission.json','r96_m0/user_authorization.json',
             'r96_m0/stage.json','r96_m0/job.json','r96_m0/archive_receipt.json',
             'r96_m0/reduction_receipt.json','r96_m0/frame_pairing_receipt.json',
             'r96_m0/validation.json']
FROZEN=['analysis/wp8_full4T_20261103','analysis/wp8_zero_mode_20261103',
        'analysis/wp8_strict_background_20261103','analysis/wp8_newton_correction_20261103',
        'analysis/wp7_shortwave_validation_20261103','analysis/wp6_4T_exploratory_20261103']


def main():
    checks={}
    # 1 frozen histories, by their own recorded hashes
    for label in FROZEN:
        p=c.ROOT/label/'provenance.json'
        if not p.exists():
            p=c.ROOT/label/'artifact_manifest.json'
            if not p.exists():continue
            man=json.loads(p.read_text())
            for n,r in man.items():
                if (c.ROOT/n).exists():checks['unchanged/'+n]=c.sha(c.ROOT/n)==r['sha256']
            continue
        j=json.loads(p.read_text())
        for n,r in j.get('files',{}).items():
            if (c.ROOT/n).exists():checks['unchanged/'+n]=c.sha(c.ROOT/n)==r['sha256']

    # 2 the superseded command/15 reports are byte-identical
    corr=json.loads((c.OUT/'corrections_v2.json').read_text(encoding='utf-8'))
    for n,r in corr['superseded_files'].items():
        checks['superseded_unchanged/'+Path(n).name]=c.sha(c.ROOT/n)==r['sha256']

    # 3 the original command authorized exactly one M1 leg.  The later M0 leg must
    # carry its own, separate user authorization rather than being silently folded
    # into command/16.
    auth=json.loads((c.OUT/'user_authorization.json').read_text(encoding='utf-8'))
    proto=json.loads((c.OUT/'r96_protocol.json').read_text(encoding='utf-8'))
    checks['authorization_matches_protocol']=auth['protocol_sha256']==c.sha(c.OUT/'r96_protocol.json')
    checks['authorization_is_one_leg']=auth['granted']['runs']==1 and auth['granted']['role']=='M1'
    checks['authority_hash']=proto['authority_sha256']==c.sha(c.ROOT/'command/16.md')
    checks['single_job_record']=len(list(c.OUT.glob('job*.json')))==1
    checks['no_forbidden_run_directory']=not any((c.OUT/n).exists() for n in
        ['run_M0','run_R192','run_8T','run_multiwave'])

    # 4 products
    for n in PRODUCTS:checks['product/'+n]=(c.OUT/n).exists()
    r96=(c.OUT/'r96_validation.json').exists()
    for n in R96_PRODUCTS:checks['r96_product/'+n]=(c.OUT/n).exists() if r96 else True
    for n in DOCS:checks['doc/'+n]=(c.OUT/n).exists()
    m0=(c.OUT/'r96_m0/validation.json').exists()
    for n in M0_PRODUCTS:checks['m0_product/'+n]=(c.OUT/n).exists() if m0 else True
    if m0:
        v0=json.loads((c.OUT/'r96_m0/validation.json').read_text(encoding='utf-8'))
        u0=json.loads((c.OUT/'r96_m0/user_authorization.json').read_text(encoding='utf-8'))
        checks['m0_has_separate_user_authorization']=(
            u0['id']=='WP9-A02' and u0['granted']['runs']==1
            and u0['granted']['role']=='M0')
        checks['m0_acceptance_passed']=v0['passed'] and not v0['failures']
        checks['m0_uses_the_same_executable']=v0['checks']['elf_is_the_R24_R48_executable']
        checks['m0_ran_on_hpcc']=v0['host']=='hpcc'
        a0=json.loads((c.OUT/'r96_m0/admission.json').read_text(encoding='utf-8'))
        checks['m0_differs_from_m1_only_in_the_perturbation']=(
            set(a0['input_key_changes'])==set(a0['r48_pair_key_changes']))
    if r96:
        val=json.loads((c.OUT/'r96_validation.json').read_text(encoding='utf-8'))
        checks['r96_acceptance_passed']=val['passed'] and not val['failures']
        checks['r96_uses_the_R24_R48_executable']=val['checks']['elf_is_the_R24_R48_executable']
        checks['r96_protocol_unchanged']=val['checks']['protocol_unchanged']

    # 5 scientific guards
    comp=json.loads((c.OUT/'resolution_comparison.json').read_text(encoding='utf-8'))
    dec=json.loads((c.OUT/'decision.json').read_text(encoding='utf-8'))
    checks['no_t_eff']=dec['t_eff'] is None
    checks['no_crystallization_label']=dec['crystallization_delay_label'] is None
    checks['no_training_labels']=dec['training_labels_created'] is False
    checks['censoring_recorded']=corr['correction_2_censoring']['status']=='right_censored'
    checks['timestep_confound_declared']='simultaneous timestep change' in [
        k['name'] for k in proto['known_confounds']]
    checks['uncertainty_declared']=all(k in comp['uncertainty'] for k in
        ['x_resolvability_floor','modal_energy_floor','timestep_prior_on_growth'])
    lsa_doc=(c.OUT/'LSA_WAVENUMBER_AUDIT_ZH.md').read_text(encoding='utf-8')
    decision_doc=(c.OUT/'DECISION_ZH.md').read_text(encoding='utf-8')
    readme_doc=(c.OUT/'README_ZH.md').read_text(encoding='utf-8')
    checks['lsa_doc_uses_current_agreement_summary']='中位数 5.44%、最大 9.78%' in lsa_doc
    checks['decision_doc_matches_machine_outcome']=f'**{dec["outcome"]}**' in decision_doc
    checks['readme_records_separate_m0_authorization']='单独授权' in readme_doc

    # 6 tests and static checks
    t=subprocess.run([sys.executable,'-m','unittest','discover','-s','scripts','-p','test_wp9.py'],
                     cwd=c.ROOT,capture_output=True,text=True)
    (c.OUT/'numerical_tests.txt').write_text(t.stdout+t.stderr,encoding='utf-8')
    checks['numerical_tests']=t.returncode==0
    s=subprocess.run([sys.executable,'-m','unittest','discover','-s','scripts','-p','test_wp6_paired_stopping.py'],
                     cwd=c.ROOT,capture_output=True,text=True)
    (c.OUT/'stopping_semantics_tests.txt').write_text(s.stdout+s.stderr,encoding='utf-8')
    checks['stopping_tests']=s.returncode==0
    p1=subprocess.run([sys.executable,'scripts/prepare_paper1_inputs.py','--verify-only'],
                      cwd=c.ROOT,capture_output=True,text=True)
    (c.OUT/'paper1_readonly_verification.txt').write_text(p1.stdout+p1.stderr,encoding='utf-8')
    checks['paper1_frozen_sources']=p1.returncode==0
    for n in SCRIPTS:
        if n.endswith('.py'):compile((c.ROOT/'scripts'/n).read_bytes(),n,'exec')
    checks['python_static_compile']=True
    checks['sbatch_static_syntax']=all(
        subprocess.run(['bash','-n',f'scripts/{n}'],cwd=c.ROOT,
                       capture_output=True,text=True).returncode==0
        for n in ['wp9_hpcc_job.sbatch','wp9_hpcc_threads.sbatch','wp9_hpcc_m0_job.sbatch'])
    ts=c.OUT/'hpcc_portability/thread_scaling.json'
    if ts.exists():
        tj=json.loads(ts.read_text(encoding='utf-8'))
        checks['threads_bitwise_identical']=tj['all_bitwise_identical']
        checks['threads_scope_is_infrastructure']='no scientific result' in tj['scope']
    port=c.OUT/'hpcc_portability/portability.json'
    if port.exists():
        pj=json.loads(port.read_text(encoding='utf-8'))
        checks['hpcc_same_executable']=pj['elf_matches_sal']
        checks['hpcc_bitwise_identical']=pj['numerically_identical']
        checks['hpcc_scope_is_infrastructure']='no scientific result' in pj['scope']

    report=dict(delivery_passed=all(checks.values()),checks_count=len(checks),
                failures=[k for k,v in checks.items() if not v],checks=checks,
                authorized_runs=1+(1 if m0 else 0),
                launched_runs=(1 if (c.OUT/'job.json').exists() else 0)+
                              (1 if (c.OUT/'r96_m0/job.json').exists() else 0),
                launched_under_command=['R96-M1'] if (c.OUT/'job.json').exists() else [],
                launched_under_separate_user_authorization=['R96-M0'] if m0 else [],
                retries=0,
                r96_present=r96,r96_m0_present=m0,
                not_launched=['R192','multiwavelength','full-domain','8T'])
    c.save(c.OUT/'delivery_verification.json',report)
    files=[c.ROOT/'command/16.md']+[c.ROOT/'scripts'/n for n in SCRIPTS]+\
          [x for x in c.OUT.rglob('*') if x.is_file() and x.name!='provenance.json']
    c.save(c.OUT/'provenance.json',dict(command='command/16.md',
        reviewed_baseline='fb4a91b2646d09394e428d408e2ffc005d4c7713',
        hosts=['sal','hpcc'],
        remote_runs=dict(
            R96_M1='/tmp/zhixin_paper2_wp9_r96_20261104/runs/WP9-R96-M1-4T',
            R96_M0='/WORK/hulm_work/zhixin/paper2_wp9_r96_m0_20261104/run'),
        remote_archives=dict(R96_M1='~/paper2_wp9_r96_20261104',
                             R96_M0='/WORK/hulm_work/zhixin/paper2_wp9_r96_m0_20261104'),
        paper1_reference=dict(module='code/lsa_unsaturated.py',
                              sha256=c.sha(c.PAPER1/'lsa_unsaturated.py'),access='read-only'),
        files={x.relative_to(c.ROOT).as_posix():dict(sha256=c.sha(x),bytes=x.stat().st_size)
               for x in files}))
    print(json.dumps(dict(delivery_passed=report['delivery_passed'],checks=len(checks),
                          failures=report['failures']),indent=2))
    assert report['delivery_passed']


if __name__=='__main__':main()
