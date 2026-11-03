"""Verify the frozen plan without writes, or restore ignored WP5 reference cache."""
import argparse,hashlib,json,subprocess,tarfile,shutil
from pathlib import Path
from single_wavelength_post import post
ROOT=Path(__file__).resolve().parents[1]
PLAN=ROOT/'analysis/wp6_bounded_nonlinear_plan_20261103'
WP5=ROOT/'analysis/single_wavelength_validation_20261103'
CACHE=ROOT/'tmp/single_wavelength_validation_20261103'
def sha(p):
    with Path(p).open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def verify():
    manifest=json.loads((PLAN/'artifact_manifest.json').read_text())
    for name,r in manifest.items():
        if name=='scripts/prepare_wp6_plan.py':
            old=subprocess.check_output(['git','show','314d75d9:'+name],cwd=ROOT)
            assert hashlib.sha256(old).hexdigest()==r['sha256']
        else:assert sha(ROOT/name)==r['sha256'],name
    for name,h in json.loads((PLAN/'source_hashes.json').read_text()).items():
        assert (ROOT/name).exists(),f'Missing cache: {name}; use --restore-reference'
        assert sha(ROOT/name)==h,name
    print('Frozen registration, products and source hashes verified; no writes.')
def restore():
    base=ROOT/'analysis/pilot0_integrated_20261102/sal_migration'
    records=json.loads((base/'field_archive.json').read_text())
    for role,r in records.items():
        archive=Path(r['local_path']);assert sha(archive)==r['sha256']
        d=CACHE/'reference'/role;d.mkdir(parents=True,exist_ok=True)
        with tarfile.open(archive,'r|gz') as tar:
            for m in tar:
                target=(d/m.name).resolve();assert target.is_relative_to(d.resolve()) and m.isfile()
                b=tar.extractfile(m).read();assert hashlib.sha256(b).hexdigest()==r['files'][m.name]
                if target.exists():assert sha(target)==r['files'][m.name]
                else:target.parent.mkdir(parents=True,exist_ok=True);target.write_bytes(b)
        for p in (base/'runs'/role).iterdir():
            if p.suffix in ['.input','.dat','.csv']:
                target=d/p.name
                if target.exists():assert sha(target)==sha(p)
                else:shutil.copyfile(p,target)
        profile=base/'runs/M1/eigenprofile.dat';target=d/'eigenprofile.dat'
        if target.exists():assert sha(target)==sha(profile)
        else:shutil.copyfile(profile,target)
        expected=json.loads((WP5/'reference'/role/'field_audit.json').read_text())['npz_sha256']
        output=d/'post'
        if all((output/n).exists() and sha(output/n)==h for n,h in expected.items()):continue
        # Never overwrite even a corrupt old cache: rebuild in a separate directory.
        rebuild=d/'post_wp6_restore';assert not rebuild.exists(),'Inspect previous restoration attempt first'
        post(d,rebuild)
        for n,h in expected.items():assert sha(rebuild/n)==h,(role,n)
        for n,h in expected.items():
            target=output/n
            if target.exists():assert sha(target)==h,'Conflicting existing cache retained'
            else:output.mkdir(exist_ok=True);shutil.copyfile(rebuild/n,target)
    verify()
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--restore-reference',action='store_true');a=p.parse_args()
    restore() if a.restore_reference else verify()
