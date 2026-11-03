"""Verify immutable references, delivered run outputs and archive hashes."""
import hashlib,json,subprocess
from pathlib import Path
import fitz
from wp6_plan_integrity import verify
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'analysis/wp6_execution_20261103';CACHE=ROOT/'tmp/wp6_execution_20261103'
def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def load(p):return json.loads(Path(p).read_text())
def main():
    verify();checks={}
    wp5=ROOT/'analysis/single_wavelength_validation_20261103'
    for name,r in load(wp5/'artifact_manifest.json').items():assert sha(ROOT/name)==r['sha256'],name
    checks['WP5_all_manifest_files_unchanged']=True
    for name,r in load(ROOT/'model/transport_wp6_v1/source_manifest.json').items():assert sha(ROOT/'model/transport_wp6_v1/source'/name)==r
    checks['isolated_source_hashes']=True
    for d in (OUT/'runs').iterdir():
        f=load(d/'finished.json');assert f['solver']['returncode']==0
        for name,r in f['files'].items():
            p=d/name
            if not p.exists():p=CACHE/'runs'/d.name/name
            assert sha(p)==r['sha256'],p
    checks['all_six_run_output_hashes']=True
    count=0
    for stage in ['smoke','bridge','extend']:
        for n,r in load(OUT/(stage+'_archive.json')).items():assert sha(CACHE/(n+'_fields.tar.gz'))==r['sha256'];count+=r['verified_vtu']
        state=load(OUT/('frodo_'+stage+'_status.json'))['manager'];assert not state['state']['active'] and all(v==0 for v in state['finished'].values())
    assert count==3607;checks['verified_archive_frames']=count;checks['all_managers_finished']=True
    assert load(OUT/'bridge_checks.json')['passed']
    e=load(OUT/'extend_checks.json');assert e['numerical_checks_passed'] and not e['passed'] and not e['checks']['M0/absolute_endpoint_screen']
    checks['background_screen_retained_not_hidden']=True
    assert all(not r['persistent_candidate'] for r in e['paired_summary']);checks['no_validated_transport_label']=e['t_eff'] is None
    for p in OUT.glob('Figure*.pdf'):
        with fitz.open(p) as doc:assert len(doc)==1 and len(doc[0].get_text())>80
    checks['three_figure_groups_inspected']=True
    checks['cold_reference_restoration']=load(OUT/'restoration_checks.json')
    checks['execution_limit']='T bridge + 2T pair only; no 4T or broader simulations'
    (OUT/'delivery_checks.json').write_text(json.dumps(checks,indent=2)+'\n')
    scripts=[p for p in (ROOT/'scripts').glob('*wp6*.py')]
    model=[p for p in (ROOT/'model/transport_wp6_v1').rglob('*') if p.is_file()]
    files=[p for p in OUT.rglob('*') if p.is_file() and p.name!='artifact_manifest.json' and p.suffix not in ['.pvd','.npz','.vtu']]+scripts+model+[ROOT/'command/7.md']
    manifest={str(p.relative_to(ROOT)).replace('\\','/'):dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(files)}
    (OUT/'artifact_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n');print(json.dumps(checks,indent=2))
if __name__=='__main__':main()
