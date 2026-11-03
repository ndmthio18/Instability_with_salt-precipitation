"""Recompute weighted reference diagnostics from hash-verified local baseline archives."""
import sys,json,tarfile,shutil,hashlib,time
from pathlib import Path
from single_wavelength_post import post
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'analysis/single_wavelength_validation_20261103';CACHE=ROOT/'tmp/single_wavelength_validation_20261103';BASE=ROOT/'analysis/pilot0_integrated_20261102/sal_migration'
def main():
 archives=json.load(open(BASE/'field_archive.json'));results={}
 for role in ['M0','M1']:
  source=BASE/'runs'/role;d=CACHE/'reference'/role;d.mkdir(parents=True,exist_ok=True);ar=Path(archives[role]['local_path']);assert hashlib.sha256(ar.read_bytes()).hexdigest()==archives[role]['sha256']
  with tarfile.open(ar) as tar:
   for member in tar.getmembers():
    target=(d/member.name).resolve();assert target.is_relative_to(d.resolve()) and member.isfile()
    with tar.extractfile(member) as stream:b=stream.read()
    assert hashlib.sha256(b).hexdigest()==archives[role]['files'][member.name]
    if not target.exists():target.write_bytes(b)
    else:assert hashlib.sha256(target.read_bytes()).hexdigest()==archives[role]['files'][member.name]
  for p in source.iterdir():
   if p.suffix in ['.input','.dat','.csv']:shutil.copyfile(p,d/p.name)
  shutil.copyfile(BASE/'runs/M1/eigenprofile.dat',d/'eigenprofile.dat')
  output=CACHE/'reference'/role/'post';results[role]=post(d,output)
  for p in output.iterdir():
   if p.suffix in ['.csv','.json']:
    dest=OUT/'reference'/role/p.name;dest.parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(p,dest)
  print(role,'weighted reference finished',flush=True)
 (OUT/'reference/recompute.json').write_text(json.dumps(results,indent=2))
if __name__=='__main__':main()
