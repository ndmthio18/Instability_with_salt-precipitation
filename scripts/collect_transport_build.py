"""Retrieve compile evidence only; never execute the compiled application."""
import base64,hashlib,json,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'analysis/transport_build_20261102'
launch=json.loads(sorted(OUT.glob('attempt_*_launch.json'))[-1].read_text())
code=r'''
import json,base64,hashlib
from pathlib import Path
attempt=Path(ATTEMPT);root=attempt.parent
if not (attempt/'result.json').exists():raise RuntimeError('Build still running or failed before report; inspect worker.log')
files={}
for name in ['result.json','build.log','configure.log','worker.log','staged_source_sha256.json','build_worker.py']:
 p=attempt/name
 if p.exists():
  b=p.read_bytes();files[name]={'sha256':hashlib.sha256(b).hexdigest(),'base64':base64.b64encode(b).decode()}
original=Path('/temp/zhixin/surrogate/code/dumux_kiemle2023b_unsat/kiemle2023b-main')
original_hashes={str(p.relative_to(original)):hashlib.sha256(p.read_bytes()).hexdigest() for p in original.rglob('*')
 if p.is_file() and 'build-cmake' not in p.parts and '.git' not in p.parts and p.suffix in ('.hh','.cc','.txt','.cmake','.module','.in','.md','.input')}
print(json.dumps({'files':files,'original_source_sha256_after_build':original_hashes}))
'''.replace('ATTEMPT',repr(launch['attempt']))
bridge='/mnt/d/3OutputFile/1Research/Article/22Surrogate_model/paper/2/scripts/frodo_bridge.py'
r=subprocess.run(['wsl','python3',bridge,'run','python3 -'],input=code.encode(),capture_output=True,timeout=60)
if r.returncode:raise RuntimeError(r.stderr.decode(errors='replace'))
j=json.loads(r.stdout);dest=OUT/Path(launch['attempt']).name;dest.mkdir(exist_ok=True)
for name,rec in j['files'].items():
 b=base64.b64decode(rec['base64']);assert hashlib.sha256(b).hexdigest()==rec['sha256'];(dest/name).write_bytes(b)
(OUT/'shared_source_after_build.json').write_text(json.dumps(j['original_source_sha256_after_build'],indent=2),encoding='utf-8')
report=json.loads((dest/'result.json').read_text())
print(json.dumps({'configure_exit':report['configure_exit'],'build_exit':report.get('build_exit'),'executable':report.get('executable'),'new_simulations':False},indent=2))
