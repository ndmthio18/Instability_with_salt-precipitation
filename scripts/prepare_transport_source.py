"""Fetch a read-only source snapshot; refuse to overwrite an existing isolated tree."""
import base64, hashlib, json, subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
DEST=ROOT/'model/transport_v1/upstream'
if DEST.exists():raise SystemExit('Existing upstream snapshot: refusing overwrite')
bridge='/mnt/d/3OutputFile/1Research/Article/22Surrogate_model/paper/2/scripts/frodo_bridge.py'
code=r'''
import base64,hashlib,json,subprocess,datetime
from pathlib import Path
root=Path('/temp/zhixin/surrogate/code/dumux_kiemle2023b_unsat/kiemle2023b-main')
files={}
for p in root.rglob('*'):
 rel=p.relative_to(root)
 if any(x in ('build-cmake','.git') for x in rel.parts):continue
 if p.is_file() and p.suffix in ('.hh','.cc','.txt','.cmake','.module','.in','.md','.input'):
  b=p.read_bytes();files[str(rel)]={'sha256':hashlib.sha256(b).hexdigest(),'base64':base64.b64encode(b).decode()}
head=subprocess.check_output(['git','-C','/temp/zhixin/surrogate','rev-parse','HEAD']).decode().strip()
print(json.dumps({'upstream_repository_head':head,'origin':'DaRUS-4610 unsaturated fork; current remote working source',
 'remote_path':str(root),'captured_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'files':files}))
'''
r=subprocess.run(['wsl','python3',bridge,'run','python3 -'],input=code.encode(),capture_output=True,timeout=100)
if r.returncode:raise RuntimeError(r.stderr.decode(errors='replace'))
j=json.loads(r.stdout)
for name,rec in j['files'].items():
 p=DEST/name
 if not p.resolve().is_relative_to(DEST.resolve()):raise ValueError(name)
 b=base64.b64decode(rec.pop('base64'));assert hashlib.sha256(b).hexdigest()==rec['sha256']
 p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(b)
(DEST.parent/'upstream_manifest.json').write_text(json.dumps(j,indent=2),encoding='utf-8')
print('Captured',len(j['files']),'source/config files; no build artifacts, no simulation')
