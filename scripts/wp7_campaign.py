"""Stage, launch once, collect and archive the bounded WP7 runs on frodo."""
import argparse,base64,hashlib,json,subprocess,tarfile
from pathlib import Path
from collect_post_onset import remote
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'analysis/wp7_shortwave_validation_20261103';CACHE=ROOT/'tmp/wp7_shortwave_validation_20261103'
REMOTE='/temp/zhixin/paper2_wp7_shortwave_20261103_r01';BRIDGE='/mnt/d/3OutputFile/1Research/Article/22Surrogate_model/paper/2/scripts/frodo_bridge.py'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,v):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(v,indent=2)+'\n')
def names(stage):return sorted(p.name for p in (OUT/'runs').iterdir() if json.loads((p/'plan.json').read_text())['stage']==stage)
def stage(which):
    files={}
    mapping={'wp7_worker.py':'wp7_worker.py','wp7_manager.py':'wp7_manager.py','wp6_post.py':'wp6_post.py','single_wavelength_post.py':'single_wavelength_post.py','wp6_timestep_manager.py':'wp6_timestep_manager.py'}
    for dst,src in mapping.items():files[dst]=base64.b64encode((ROOT/'scripts'/src).read_bytes()).decode()
    files['preregistration.json']=base64.b64encode((OUT/'preregistration.json').read_bytes()).decode()
    for n in names(which):
        for f in ['params.input','plan.json','eigenprofile.dat']:files[f'runs/{n}/{f}']=base64.b64encode((OUT/'runs'/n/f).read_bytes()).decode()
    code=r'''import base64,hashlib,json,os,shutil
from pathlib import Path
r=Path(REMOTE);r.mkdir(exist_ok=True)
if not (r/'bundle_ready.json').exists():
 shutil.copytree('/temp/zhixin/paper2_wp6_4T_exploratory_20261103_r01/bundle',r/'bundle')
 (r/'bundle_ready.json').write_text(json.dumps(dict(elf=hashlib.sha256((r/'bundle/bin/paper2_transport').read_bytes()).hexdigest())))
for name,s in FILES.items():
 p=r/name;b=base64.b64decode(s);p.parent.mkdir(parents=True,exist_ok=True);assert not p.exists() or p.read_bytes()==b,name
 if not p.exists():p.write_bytes(b)
avail=int(next(x.split()[1] for x in Path('/proc/meminfo').read_text().splitlines() if x.startswith('MemAvailable:')))*1024;free=os.statvfs(r).f_bavail*os.statvfs(r).f_frsize
assert avail>16*1024**3 and free>200*1024**3 and os.getloadavg()[0]<8
print(json.dumps(dict(root=str(r),available_memory_bytes=avail,free_disk_bytes=free,load=os.getloadavg(),elf=hashlib.sha256((r/'bundle/bin/paper2_transport').read_bytes()).hexdigest())))
'''.replace('REMOTE',repr(REMOTE)).replace('FILES',repr(files))
    result=remote(code);save(OUT/f'stage_{which}.json',result);print(result)
def launch(which):
    ns=names(which);assert json.loads((OUT/'static_checks.json').read_text())['passed']
    if which=='r48pair':assert json.loads((OUT/'benchmark_admission.json').read_text())['admitted']
    code=r'''import json,subprocess,sys
from pathlib import Path
r=Path(REMOTE);ns=NAMES
assert all(not (r/'runs'/n/'started.json').exists() for n in ns)
claim=r/(STAGE+'_launch_claim.json');claim.open('x').write('{}')
log=(r/(STAGE+'_manager.log')).open('w');p=subprocess.Popen([sys.executable,str(r/'wp7_manager.py'),str(r),STAGE,json.dumps(ns)],stdout=log,stderr=subprocess.STDOUT,start_new_session=True,stdin=subprocess.DEVNULL)
result=dict(pid=p.pid,stage=STAGE,names=ns);claim.write_text(json.dumps(result));print(json.dumps(result))
'''.replace('REMOTE',repr(REMOTE)).replace('NAMES',repr(ns)).replace('STAGE',repr(which))
    result=remote(code);save(OUT/f'launch_{which}.json',result);print(result)
def collect(which):
    ns=names(which);code=r'''import base64,hashlib,json
from pathlib import Path
r=Path(REMOTE);result=dict(runs={},manager={})
for suffix in ['state','resources','finished','wall_limit_result']:
 p=r/(STAGE+('_wall_limit_result.json' if suffix=='wall_limit_result' else '_manager_'+suffix+'.json'))
 if p.exists():result['manager'][suffix]=json.loads(p.read_text())
for n in NAMES:
 d=r/'runs'/n;info={}
 for x in ['started.json','solver_finished.json','finished.json','failure.json','stop_record.json']:
  if (d/x).exists():info[x]=json.loads((d/x).read_text())
 if (d/'solver_finished.json').exists():info['files']={str(p.relative_to(d)):dict(data=base64.b64encode(p.read_bytes()).decode(),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in d.rglob('*') if p.is_file() and p.suffix in ['.csv','.json','.log','.pvd','.npz']}
 result['runs'][n]=info
print(json.dumps(result))
'''.replace('REMOTE',repr(REMOTE)).replace('NAMES',repr(ns)).replace('STAGE',repr(which))
    result=remote(code,timeout=300)
    for n,info in result['runs'].items():
        for rel,v in info.pop('files',{}).items():
            b=base64.b64decode(v['data']);assert hashlib.sha256(b).hexdigest()==v['sha256'];p=(CACHE/'runs'/n/rel) if rel.endswith('.npz') else (OUT/'runs'/n/rel);p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(b)
    save(OUT/f'status_{which}.json',result);print(result['manager'])
def archive(which):
    receipt={};CACHE.mkdir(parents=True,exist_ok=True)
    for n in names(which):
        d=OUT/'runs'/n
        if not (d/'finished.json').exists():continue
        rec=json.loads((d/'finished.json').read_text())['archive'];dest=CACHE/f'{n}_fields.tar.gz'
        if not dest.exists() or sha(dest)!=rec['sha256']:
            subprocess.run(['wsl','python3',BRIDGE,'pull',rec['path'],f'/mnt/d/3OutputFile/1Research/Article/22Surrogate_model/paper/2/tmp/wp7_shortwave_validation_20261103/{dest.name}'],check=True,timeout=900)
        assert sha(dest)==rec['sha256'];expected={Path(q['path']).name:q['sha256'] for q in json.loads((d/'post/field_hashes.json').read_text())}
        with tarfile.open(dest,'r|gz') as tar:
            seen=set()
            for m in tar:
                if m.name in expected:
                    with tar.extractfile(m) as f:assert hashlib.file_digest(f,'sha256').hexdigest()==expected[m.name]
                    seen.add(m.name)
        assert seen==set(expected);receipt[n]=dict(**rec,local_path=str(dest),verified_vtu=len(seen))
    save(OUT/f'archive_{which}.json',receipt);print('archived',len(receipt))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['stage','launch','collect','archive']);p.add_argument('stage',choices=['preflight','r48pair']);a=p.parse_args();globals()[a.action](a.stage)
