"""Verify the bounded WP7 delivery and write its artifact manifest."""
import hashlib,json,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'analysis/wp7_shortwave_validation_20261103'
MODEL=ROOT/'model/transport_wp7_localflux_v1'
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def load(path):return json.loads(path.read_text())
checks={}
prereg=sha(OUT/'preregistration.json');e02=sha(OUT/'resource_amendment_E02.json');e03=sha(OUT/'runtime_relocation_E03.json')
for run in (OUT/'runs').iterdir():
    if not run.is_dir():continue
    plan=load(run/'plan.json');checks[run.name+'/input_hash']=sha(run/'params.input')==plan['input_sha256'];checks[run.name+'/profile_hash']=sha(run/'eigenprofile.dat')==plan['profile_sha256'];checks[run.name+'/preregistration_hash']=plan['preregistration_sha256']==prereg
    if plan['stage']=='preflight':checks[run.name+'/resource_amendment']=plan['resource_amendment_sha256']==e02;checks[run.name+'/runtime_relocation']=plan['runtime_relocation_sha256']==e03;checks[run.name+'/completed']=load(run/'finished.json')['solver']['returncode']==0
    else:checks[run.name+'/not_launched']=not (run/'started.json').exists() and not (run/'finished.json').exists()
for run in ['WP7-R48-M1-bench10','WP7-1D-background-T']:
    receipt=load(OUT/'sal_collect_preflight_retry.json')['runs'][run];archive=ROOT/receipt['local_archive'];checks[run+'/archive_hash']=receipt['complete'] and sha(archive)==receipt['archive_sha256'];checks[run+'/independent_checks']=load(OUT/'runs'/run/'post/independent_checks.json')['passed']
failed=list((OUT/'failed_attempts').iterdir());checks['two_zero_step_failures_retained']=len(failed)==2 and all(load(x/'solver_finished.json')['newton_iterations']==0 and load(x/'solver_finished.json')['stop'] is None for x in failed)
manifest=load(MODEL/'source_manifest.json');checks['source_manifest']=all(sha(MODEL/'source'/name)==rec['sha256'] and (MODEL/'source'/name).stat().st_size==rec['bytes'] for name,rec in manifest.items())
build=load(OUT/'build/attempt_01/result.json');checks['isolated_build']=build['configure_exit']==0 and build['build_exit']==0 and build['source_sha256']=={name:rec['sha256'] for name,rec in manifest.items()}
checks['local_flux_static']=load(OUT/'local_flux_static_checks.json')['passed'];checks['local_flux_bridge']=load(OUT/'local_flux_smoke/bridge_result.json')['passed']
checks['background_failure_retained']=not load(OUT/'background_validation.json')['overall_preregistered_gate']
decision=load(OUT/'decision.json');checks['claim_boundary']=decision['R96']=='NO-GO' and decision['multiwave']=='NO-GO' and decision['later_8T']=='NO-GO' and decision['t_eff'] is None and decision['training_labels_created'] is False
old_manifest=ROOT/'analysis/wp6_4T_exploratory_20261103/artifact_manifest.json';head=subprocess.check_output(['git','show','HEAD:analysis/wp6_4T_exploratory_20261103/artifact_manifest.json'],cwd=ROOT);checks['WP6_artifact_manifest_unchanged']=hashlib.sha256(old_manifest.read_bytes()).digest()==hashlib.sha256(head).digest()
required=['R48_PLAN_AND_ADMISSION_ZH.md','R24_R48_COMPARISON_ZH.md','BACKGROUND_1D_VALIDATION_ZH.md','LOCAL_SCVF_DIAGNOSTIC_ZH.md','FINAL_REVIEW_ZH.md','README_ZH.md']
checks['reports_complete']=all((OUT/name).stat().st_size>500 for name in required)
result={'passed':all(checks.values()),'checks':checks};(OUT/'delivery_checks.json').write_text(json.dumps(result,indent=2)+'\n');assert result['passed'],[k for k,v in checks.items() if not v]
include=[p for p in OUT.rglob('*') if p.is_file() and p.name!='artifact_manifest.json']+[ROOT/'command/11.md']
include += [p for p in (ROOT/'model/transport_wp7_localflux_v1').rglob('*') if p.is_file()]
include += [ROOT/'scripts'/name for name in ['prepare_wp7.py','wp7_worker.py','wp7_manager.py','wp7_campaign.py','wp7_benchmark_admission.py','wp7_sal_worker.py','wp7_sal_job.sh','wp7_sal_campaign.py','orted_relocatable.sh','build_wp7_remote.py','collect_wp7_build.py','make_wp7_source_manifest.py','check_wp7_localflux_static.py','run_wp7_localflux_smoke.py','analyze_wp7.py','verify_wp7_delivery.py']]
records={str(p.relative_to(ROOT)).replace('\\','/'):dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(set(include))}
(OUT/'artifact_manifest.json').write_text(json.dumps(records,indent=2)+'\n')
print(json.dumps({'passed':True,'checks':len(checks),'artifacts':len(records)},indent=2))
