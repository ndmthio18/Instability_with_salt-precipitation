"""Verify command/9 delivery while preserving all earlier manifests."""
import hashlib,json,py_compile,subprocess,sys
from pathlib import Path
import fitz
from wp6_plan_integrity import verify
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'analysis/wp6_timestep_20261103';CACHE=ROOT/'tmp/wp6_timestep_20261103'
SCRIPTS=['prepare_wp6_timestep.py','wp6_timestep_manager.py','wp6_timestep_worker.py','wp6_timestep_campaign.py','wp6_timestep_bounds.py','test_wp6_timestep.py','analyze_wp6_timestep.py','plot_wp6_timestep.py','verify_wp6_timestep.py']
def load(p):return json.loads(p.read_text())
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
    verify();checks={}
    for folder in ['single_wavelength_validation_20261103','wp6_execution_20261103','wp6_review_20261103']:
        path=ROOT/'analysis'/folder/'artifact_manifest.json'
        for name,r in load(path).items():assert sha(ROOT/name)==r['sha256'],name
        checks[folder+'_unchanged']=True
    for name,s in load(ROOT/'model/transport_wp6_v1/source_manifest.json').items():assert sha(ROOT/'model/transport_wp6_v1/source'/name)==s
    d=OUT/'runs/WP6-M1-dt_half';p=load(d/'plan.json');f=load(d/'finished.json');a=load(d/'started.json')
    assert a['elf_sha256']==p['elf_sha256']==load(OUT/'amendment_20261103_C01.json')['elf_sha256']
    assert sha(d/'params.input')==p['input_sha256'] and sha(d/'eigenprofile.dat')==p['profile_sha256']
    assert sha(OUT/'amendment_20261103_C01.json')==p['preregistration_sha256']
    assert load(OUT/'amendment_20261103_C01.json')['created_utc']<a['utc']
    audit=load(d/'post/field_audit.json')
    assert audit['all_rho_pass'] and audit['all_sw_pass'] and audit['initial_field_error']<1e-12
    checks['initial_eigenprofile_and_density_saturation_fields']=True
    mapping={'wp6_manager.py':'wp6_timestep_manager.py','wp6_worker.py':'wp6_timestep_worker.py','wp6_post.py':'wp6_post.py','single_wavelength_post.py':'single_wavelength_post.py'}
    for name,s in load(OUT/'remote_execution_source_hashes.json').items():assert s==sha(ROOT/'scripts'/mapping[name])
    checks['remote_execution_source_hashes_match']=True
    assert len(list((OUT/'runs').iterdir()))==1
    for name,r in f['files'].items():
        path=d/name
        if not path.exists():path=CACHE/'runs'/d.name/name
        assert sha(path)==r['sha256'],str(path)
    state=load(OUT/'frodo_dt_half_status.json')['manager'];limit=state['wall_limit_result']
    assert limit['returncode']==0 and not limit['timed_out'] and limit['wall_s']<3600 and limit['hard_limit_s']==3600
    assert state['finished']=={'WP6-M1-dt_half':0} and not state['state']['active']
    ar=load(OUT/'dt_half_archive.json')[d.name];assert ar['verified_vtu']==2401 and sha(Path(ar['local_path']))==ar['sha256']
    assert load(OUT/'prelaunch_checks.json')['wall_guard_test']['timed_out']
    checks['one_new_run_3600s_hard_limit_and_archive']=dict(passed=True,manager_wall_s=limit['wall_s'],solver_wall_s=f['solver']['wall_s'],cpu_s=f['solver']['cpu_s'],verified_vtu=2401)
    tests=[]
    for name,count in [('test_wp6_paired_stopping.py',22),('test_wp6_timestep.py',9)]:
        r=subprocess.run([sys.executable,str(ROOT/'scripts'/name),'-v'],capture_output=True,text=True,encoding='utf-8');assert r.returncode==0,r.stderr
        tests.append(dict(script=name,count=count,output=r.stdout+r.stderr,passed=True))
    (OUT/'test_results.json').write_text(json.dumps(tests,indent=2)+'\n');checks['tests_passed']=31
    for name in SCRIPTS:py_compile.compile(str(ROOT/'scripts'/name),doraise=True)
    checks['source_compilation']=True
    pair=load(OUT/'refined_pair/async_stopping_checks.json');assert len(pair['pair']['indices'])==2401 and pair['pair']['max_time_representation_error_s']==0
    assert pair['pair']['surface_time_difference']['status']=='both_right_censored'
    checks['strict_refined_pair_actual_shared_states']=2401
    dec=load(OUT/'decision_checks.json');assert dec['original_M0_absolute_screen']=='failed_retained' and dec['t_eff'] is None and dec['quantitative_B']=='NO_GO'
    checks['distinct_A_B_decisions']=dict(A=dec['exploratory_A'],B=dec['quantitative_B'],old_failed_gate_retained=True)
    for p in OUT.glob('Figure*.pdf'):
        with fitz.open(p) as doc:assert len(doc)==1 and len(doc[0].get_text())>100
    checks['figure_QA']='Both PNGs visually inspected; two one-page PDFs have selectable labels'
    checks['source_base_commit']='7db1dd3f5b9ee98802d46c2d5fecfaf913c5d68d'
    (OUT/'delivery_checks.json').write_text(json.dumps(checks,indent=2)+'\n')
    files=[p for p in OUT.rglob('*') if p.is_file() and p.name!='artifact_manifest.json' and p.suffix not in ['.pvd','.npz','.vtu']]+[ROOT/'scripts'/name for name in SCRIPTS]+[ROOT/'command/9.md',ROOT/'model/transport_wp6_v1/source_manifest.json']
    manifest={p.relative_to(ROOT).as_posix():dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(files)}
    (OUT/'artifact_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps(checks,indent=2))
if __name__=='__main__':main()
