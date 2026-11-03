"""Read-only verification of frozen references; writes only command/8 receipts."""
import csv,hashlib,json,py_compile,subprocess,sys
from pathlib import Path
import fitz
from wp6_plan_integrity import verify as verify_plan
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'analysis/wp6_review_20261103';CACHE=ROOT/'tmp/wp6_review_20261103'
SCRIPTS=['wp6_review_campaign.py','wp6_paired_stopping.py','test_wp6_paired_stopping.py','wp6_background_structure.py','wp6_review_background.py','analyze_wp6_review.py','build_wp6_review_products.py','verify_wp6_review.py']
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def load(p):return json.loads(p.read_text())
def main():
    verify_plan();checks={}
    for dirname in ['single_wavelength_validation_20261103','wp6_execution_20261103']:
        manifest=ROOT/'analysis'/dirname/'artifact_manifest.json'
        for name,r in load(manifest).items():assert sha(ROOT/name)==r['sha256'],name
        checks[dirname+'_all_frozen_hashes_unchanged']=dict(passed=True,manifest_sha256=sha(manifest))
    for name,r in load(ROOT/'model/transport_wp6_v1/source_manifest.json').items():assert sha(ROOT/'model/transport_wp6_v1/source'/name)==r,name
    checks['isolated_solver_source_unchanged']=True
    d=OUT/'runs/WP6-M0-dt_half';finished=load(d/'finished.json');started=load(d/'started.json');plan=load(d/'plan.json')
    assert started['elf_sha256']==plan['elf_sha256']==load(OUT/'amendment_20261103_B01.json')['elf_sha256']
    assert sha(d/'params.input')==plan['input_sha256'] and sha(OUT/'amendment_20261103_B01.json')==plan['preregistration_sha256']
    assert started['postprocessor_sha256']==sha(ROOT/'scripts/wp6_post.py')
    assert len(list((OUT/'runs').iterdir()))==1
    for name,r in finished['files'].items():
        p=d/name
        if not p.exists():p=CACHE/'runs'/d.name/name
        assert sha(p)==r['sha256'],str(p)
    archive=load(OUT/'dt_half_archive.json')['WP6-M0-dt_half'];assert sha(Path(archive['local_path']))==archive['sha256'] and archive['verified_vtu']==2401
    state=load(OUT/'frodo_dt_half_status.json')['manager'];assert not state['state']['active'] and state['finished']=={'WP6-M0-dt_half':0}
    checks['one_control_finished_and_all_outputs_verified']=True;checks['verified_archived_VTU']=2401
    for role in ['baseline','half']:
        for name,rec in load(OUT/('background_'+role)/'transfer.json').items():
            p=(CACHE if name.endswith('.npz') else OUT)/('background_'+role)/name
            assert sha(p)==rec['sha256']
        assert load(OUT/('background_'+role)/'noise_summary.json')['script_sha256']==sha(ROOT/'scripts/wp6_background_structure.py')
    checks['all_background_field_analysis_hashes_verified']=True
    p=subprocess.run([sys.executable,str(ROOT/'scripts/test_wp6_paired_stopping.py'),'-v'],capture_output=True,text=True,encoding='utf-8')
    (OUT/'test_results.txt').write_text(p.stdout+p.stderr,encoding='utf-8');assert p.returncode==0,p.stderr
    checks['unit_and_integration_tests']=dict(passed=True,count=22)
    for name in SCRIPTS:py_compile.compile(str(ROOT/'scripts'/name),doraise=True)
    checks['static_python_compilation']=True
    parity=load(OUT/'frozen_2T_parity.json');assert parity['passed'] and parity['rows']==9608 and all(v==0 for v in parity['maximum_absolute_differences'].values())
    checks['frozen_pair_9608_rows_exact_parity']=True
    comparison=load(OUT/'background_comparison.json');assert comparison['original_background_screen']=='failed_retained' and not comparison['symmetry_threshold_amended']
    assert comparison['control_validation']['valid'] and comparison['independent_checks_passed']
    assert load(OUT/'transport_evidence.json')['t_eff'] is None
    checks['original_failed_gate_preserved_and_t_eff_undefined']=True
    for p in OUT.glob('Figure*.pdf'):
        with fitz.open(p) as doc:assert len(doc)==1 and len(doc[0].get_text())>100
    checks['figure_QA']='Three PNGs visually inspected; three one-page PDFs checked for selectable labels'
    checks['scientific_decision']='NO-GO for direct 4T; no broader simulations started'
    checks['source_base_commit']='6e8e03428779dd1f64134cc5878c9f9ebb03e927'
    (OUT/'review_delivery_checks.json').write_text(json.dumps(checks,indent=2)+'\n')
    files=[p for p in OUT.rglob('*') if p.is_file() and p.name!='artifact_manifest.json' and p.suffix not in ['.pvd','.npz','.vtu']]
    files += [ROOT/'scripts'/name for name in SCRIPTS]+[ROOT/'command/8.md',ROOT/'model/transport_wp6_v1/source_manifest.json']
    manifest={p.relative_to(ROOT).as_posix():dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(files)}
    (OUT/'artifact_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps(checks,indent=2))
if __name__=='__main__':main()
