"""Collect small products from the completed historical analysis, without simulation."""
import hashlib
import json
from pathlib import Path
import subprocess
import tarfile

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'analysis/dense_history_20261102'
BRIDGE='/mnt/d/3OutputFile/1Research/Article/22Surrogate_model/paper/2/scripts/frodo_bridge.py'
launch=OUT/'launch_r02.json' if (OUT/'launch_r02.json').exists() else OUT/'launch.json'
remote=json.loads(launch.read_text())['remote_root']
code='''
import hashlib,json,tarfile
from pathlib import Path
root=Path(REMOTE)
finished=json.loads((root/'finished.json').read_text())
files={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in root.rglob('*') if p.is_file()}
(root/'product_manifest.json').write_text(json.dumps(files,indent=2))
with tarfile.open(str(root)+'.tar.gz','w:gz') as t:
 for p in root.iterdir():t.add(p,arcname=p.name)
print(json.dumps(dict(finished=finished,sha256=hashlib.sha256(Path(str(root)+'.tar.gz').read_bytes()).hexdigest())))
'''.replace('REMOTE',repr(remote))
p=subprocess.run(['wsl','python3',BRIDGE,'run','python3 -'],input=code.encode(),capture_output=True,timeout=60)
if p.returncode:raise RuntimeError(p.stderr.decode(errors='replace'))
result=json.loads(p.stdout)
dest=ROOT/'data/generated/dense_history_20261102.tar.gz';dest.parent.mkdir(parents=True,exist_ok=True)
wsl='/mnt/'+dest.drive[0].lower()+dest.as_posix()[2:]
subprocess.run(['wsl','python3',BRIDGE,'pull',remote+'.tar.gz',wsl],check=True,timeout=60)
assert hashlib.sha256(dest.read_bytes()).hexdigest()==result['sha256']
evidence=OUT/'evidence';evidence.mkdir(exist_ok=True)
with tarfile.open(dest) as t:
    for m in t.getmembers():
        if m.isfile() and (evidence/m.name).exists():assert (evidence/m.name).read_bytes()==t.extractfile(m).read(),m.name
    t.extractall(evidence,filter='data')
manifest=json.loads((evidence/'product_manifest.json').read_text())
for name,digest in manifest.items():assert hashlib.sha256((evidence/name).read_bytes()).hexdigest()==digest,name
assert json.loads((evidence/'plan.json').read_text())==json.loads((OUT/'plan.json').read_text())
revision=evidence/'worker_revision.json'
expected_worker=json.loads(revision.read_text())['actual_worker_sha256'] if revision.exists() else json.loads((OUT/'plan.json').read_text())['worker_sha256']
assert hashlib.sha256((evidence/'worker.py').read_bytes()).hexdigest()==expected_worker
result['verified_products']=len(manifest)
(OUT/'transfer_check.json').write_text(json.dumps(result,indent=2))
print(json.dumps(result))
