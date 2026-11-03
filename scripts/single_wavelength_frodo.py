"""Frodo migration after proven cancellation of the unstarted sal array."""
import base64,hashlib,json,subprocess,tarfile,argparse
from pathlib import Path
from collect_post_onset import remote
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'analysis/single_wavelength_validation_20261103';CACHE=ROOT/'tmp/single_wavelength_validation_20261103';REMOTE='/temp/zhixin/paper2_wp5_20261103_r02';BRIDGE='/mnt/d/3OutputFile/1Research/Article/22Surrogate_model/paper/2/scripts/frodo_bridge.py'
def save(p,v):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(v,indent=2))
def stage(stage):
 assert json.load(open(OUT/'sal_cancelled_before_start.json'))['started_files']==0
 remote('import json\nfrom pathlib import Path\nPath('+repr(REMOTE)+').mkdir(exist_ok=True)\nprint(json.dumps(dict(ok=True)))')
 flag=OUT/'frodo_bundle_staged_r02.json'
 if not flag.exists():
  subprocess.run(['wsl','python3',BRIDGE,'push','/mnt/d/3OutputFile/1Research/Article/22Surrogate_model/paper/2/tmp/sal_migration/bundle.tar.gz','/mnt/d/3OutputFile/1Research/Article/22Surrogate_model/paper/2/tmp/sal_migration/runtime_support.tar.gz',REMOTE+'/'],check=True,timeout=180)
 files={}
 for p in [OUT/'preregistration.json',ROOT/'scripts/single_wavelength_frodo_worker.py',ROOT/'scripts/single_wavelength_manager.py',ROOT/'scripts/single_wavelength_post.py']:
  files[p.name]=base64.b64encode(p.read_bytes()).decode()
 for d in (OUT/'runs').iterdir():
  if d.name.endswith('-'+stage):
   for name in ['params.input','plan.json','eigenprofile.dat']:files['runs/'+d.name+'/'+name]=base64.b64encode((d/name).read_bytes()).decode()
 launcher=(ROOT/'analysis/pilot0_integrated_20261102/sal_migration/execution/launch_solver.py').read_text()
 code='''import json,base64,tarfile,hashlib,os
from pathlib import Path
r=Path(REMOTE)
if not (r/'bundle_ready.json').exists():
 for name,expected in ARCHIVES.items():assert hashlib.sha256((r/name).read_bytes()).hexdigest()==expected
 with tarfile.open(r/'bundle.tar.gz') as tar:tar.extractall(r,filter='data')
 with tarfile.open(r/'runtime_support.tar.gz') as tar:tar.extractall(r/'bundle',filter='data')
 (r/'bundle/launch_solver.py').write_text(LAUNCHER)
 (r/'bundle/bin/orted').write_text(ORTED_WRAPPER);(r/'bundle/bin/orted').chmod(0o755)
 (r/'bundle_ready.json').write_text(json.dumps(ARCHIVES))
for name,b in FILES.items():
 p=r/name;p.parent.mkdir(parents=True,exist_ok=True);data=base64.b64decode(b)
 assert not p.exists() or p.read_bytes()==data,'refuse overwrite '+name
 if not p.exists():p.write_bytes(data)
print(json.dumps(dict(root=str(r),load=os.getloadavg(),bundle_sha256=hashlib.sha256((r/'bundle/bin/paper2_transport').read_bytes()).hexdigest())))
'''.replace('REMOTE',repr(REMOTE)).replace('FILES',repr(files)).replace('LAUNCHER',repr(launcher)).replace('ORTED_WRAPPER',repr('#!/bin/sh\nexec '+REMOTE+'/bundle/lib/ld-linux-x86-64.so.2 --library-path '+REMOTE+'/bundle/lib '+REMOTE+'/bundle/bin/orted.real \"$@\"\n')).replace('ARCHIVES',repr({n:hashlib.sha256((ROOT/'tmp/sal_migration'/n).read_bytes()).hexdigest() for n in ['bundle.tar.gz','runtime_support.tar.gz']}))
 res=remote(code);save(flag,res);print(res)
def launch(stage):
 assert json.load(open(OUT/'configuration_checks.json'))['passed']
 if stage=='full':assert json.load(open(OUT/'short_validation.json'))['advance_to_full']
 names=sorted(p.name for p in (OUT/'runs').iterdir() if p.name.endswith('-'+stage))
 code='''import json,subprocess,sys,os
from pathlib import Path
r=Path(REMOTE)
assert (r/'bundle_ready.json').exists() and (r/'single_wavelength_manager.py').exists()
assert all((r/'runs'/n/'plan.json').exists() and not (r/'runs'/n/'started.json').exists() for n in NAMES)
with (r/(STAGE+'_launch_claim.json')).open('x') as f:f.write('{}')
log=(r/(STAGE+'_manager.log')).open('w');p=subprocess.Popen([sys.executable,str(r/'single_wavelength_manager.py'),str(r),STAGE,json.dumps(NAMES)],stdout=log,stderr=subprocess.STDOUT,start_new_session=True,stdin=subprocess.DEVNULL)
res=dict(pid=p.pid,stage=STAGE,host='frodo',names=NAMES);(r/(STAGE+'_launch_claim.json')).write_text(json.dumps(res));print(json.dumps(res))
'''.replace('REMOTE',repr(REMOTE)).replace('STAGE',repr(stage)).replace('NAMES',repr(names))
 res=remote(code);save(OUT/('frodo_launch_'+stage+'.json'),res);print(res)
def collect(stage):
 names=sorted(p.name for p in (OUT/'runs').iterdir() if p.name.endswith('-'+stage))
 code='''import json,base64,hashlib
from pathlib import Path
r=Path(REMOTE);result=dict(runs={},manager={})
for n in ['state','resources','finished']:
 p=r/(STAGE+'_manager_'+n+'.json')
 if p.exists():result['manager'][n]=json.loads(p.read_text())
for name in NAMES:
 d=r/'runs'/name;info={}
 for n in ['started.json','stop_record.json','solver_finished.json','finished.json','failure.json']:
  if (d/n).exists():info[n]=json.loads((d/n).read_text())
 if (d/'solver_finished.json').exists():info['files']={str(p.relative_to(d)):dict(data=base64.b64encode(p.read_bytes()).decode(),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in d.rglob('*') if p.is_file() and p.suffix in ['.csv','.json','.log','.pvd','.npz']}
 result['runs'][name]=info
print(json.dumps(result))
'''.replace('REMOTE',repr(REMOTE)).replace('STAGE',repr(stage)).replace('NAMES',repr(names))
 res=remote(code)
 for name,info in res['runs'].items():
  for rel,v in info.pop('files',{}).items():
   b=base64.b64decode(v['data']);assert hashlib.sha256(b).hexdigest()==v['sha256'];p=(CACHE/'runs'/name/rel) if rel.endswith('.npz') else (OUT/'runs'/name/rel);p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(b)
 save(OUT/('frodo_'+stage+'_status.json'),res);print(res['manager']);print({k:v.get('stop_record.json',{}).get('accepted_steps') for k,v in res['runs'].items()})
def archive(stage):
 result={}
 for d in sorted((OUT/'runs').iterdir()):
  if not d.name.endswith('-'+stage) or not (d/'finished.json').exists():continue
  rec=json.load(open(d/'finished.json'))['archive'];dest=CACHE/(d.name+'_fields.tar.gz')
  if not dest.exists() or hashlib.sha256(dest.read_bytes()).hexdigest()!=rec['sha256']:
   destwsl='/mnt/d/3OutputFile/1Research/Article/22Surrogate_model/paper/2/tmp/single_wavelength_validation_20261103/'+dest.name
   subprocess.run(['wsl','python3',BRIDGE,'pull',rec['path'],destwsl],check=True,timeout=300)
  assert hashlib.sha256(dest.read_bytes()).hexdigest()==rec['sha256'];expected={Path(q['path']).name:q['sha256'] for q in json.load(open(d/'post/field_hashes.json'))}
  with tarfile.open(dest,'r|gz') as tar:
   seen=set()
   for member in tar:
    if member.name in expected:
     with tar.extractfile(member) as f:assert hashlib.file_digest(f,'sha256').hexdigest()==expected[member.name]
     seen.add(member.name)
   assert seen==set(expected)
  result[d.name]=dict(**rec,local_path=str(dest),verified_vtu=len(expected))
 save(OUT/(stage+'_archive.json'),result);print('archived',len(result))
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('action',choices=['stage','launch','collect','archive']);p.add_argument('stage',choices=['short','full','bench']);a=p.parse_args();globals()[a.action](a.stage)
