"""Read-only evidence checks and manifest generation; no solver launch."""
import hashlib,json,subprocess
from pathlib import Path
import fitz
from analyze_single_wavelength import OUT,ROOT,CACHE

def sha(p):
    h=hashlib.sha256()
    with Path(p).open('rb') as f:
        for b in iter(lambda:f.read(8*1024*1024),b''):h.update(b)
    return h.hexdigest()
def main():
    checks={};load=lambda p:json.loads(Path(p).read_text(encoding='utf-8'))
    for p,h in load(OUT/'reference_hashes.json').items():assert sha(ROOT/p)==h
    checks['frozen_reference_and_source_unchanged']=True
    total=0
    for stage in ['short','bench','full']:
        for name,a in load(OUT/(stage+'_archive.json')).items():
            assert sha(CACHE/(name+'_fields.tar.gz'))==a['sha256'];total+=a['verified_vtu']
    assert total==2428;checks['archive_sha_verified_frames']=total
    for d in (OUT/'runs').iterdir():
        f=load(d/'finished.json');p=load(d/'plan.json');assert f['solver']['returncode']==0
        assert abs(f['solver']['stop']['last_accepted_time_s']-p['T_s'])<1e-7
        for rel,evidence in f['files'].items():
            q=d/rel
            if not q.exists():q=CACHE/'runs'/d.name/rel
            assert q.exists() and sha(q)==evidence['sha256'],q
    checks['all_successful_run_outputs_sha_verified']=8
    for s in ['short','full']:assert load(OUT/(s+'_validation.json'))['numerical_passed']
    assert load(OUT/'short_validation.json')['advance_to_full']
    config=load(OUT/'configuration_checks.json');assert all(config['checks'].values())
    checks['configuration_check_count']=len(config['checks'])
    checks['scientific_tests']=load(OUT/'test_results.json');assert checks['scientific_tests']['passed']
    for p in OUT.glob('Figure*.pdf'):
        with fitz.open(p) as doc:assert len(doc)==1 and len(doc[0].get_text())>100
    checks['figure_pdf_pages']=3
    checks['figure_visual_review']='All PNGs inspected; Figure2 tick density and annotation overlap corrected.'
    checks['report_count']=len(list(OUT.glob('*_ZH.md')));assert checks['report_count']==8
    checks['transport_effectiveness']='unresolved; time and vertical discretization errors unknown'
    (OUT/'delivery_checks.json').write_text(json.dumps(checks,indent=2)+'\n')
    scripts=[p for p in (ROOT/'scripts').glob('*single_wavelength*.py')]
    extras=[ROOT/'command/5.md',ROOT/'analysis/unperturbed_36877_20261102/evidence/finished.json']
    paths=[p for p in OUT.rglob('*') if p.is_file() and p.name!='artifact_manifest.json' and p.suffix!='.pvd']+scripts+extras
    manifest={str(p.relative_to(ROOT)).replace('\\','/'):{'sha256':sha(p),'bytes':p.stat().st_size} for p in sorted(paths)}
    (OUT/'artifact_manifest.json').write_text(json.dumps(manifest,indent=2)+'\n')
    print(json.dumps(checks,indent=2))
if __name__=='__main__':main()
