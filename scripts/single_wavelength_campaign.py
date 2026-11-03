"""WP5 staging, explicit submission and collection; no implicit retry or extension."""
import argparse,base64,hashlib,json,subprocess,tarfile
from pathlib import Path
from pilot0_sal_campaign import remote
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'analysis/single_wavelength_validation_20261103';REMOTE='/tmp/zhixin_paper2_wp5_20261103_r01';CACHE=ROOT/'tmp/single_wavelength_validation_20261103';CACHE.mkdir(exist_ok=True,parents=True)
def save(p,v):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(v,indent=2),encoding='utf-8')
def prepare(stage):
 runs=sorted(p for p in (OUT/'runs').iterdir() if p.name.endswith('-'+stage));files={}
 for p in [OUT/'preregistration.json',*[ROOT/'scripts'/n for n in ['single_wavelength_worker.py','single_wavelength_post.py']]]:files[p.name]=base64.b64encode(p.read_bytes()).decode()
 for d in runs:
  for p in d.iterdir():
   if p.name in ['params.input','plan.json','eigenprofile.dat']:files['runs/'+d.name+'/'+p.name]=base64.b64encode(p.read_bytes()).decode()
 names=[p.name for p in runs]
 code='''import base64,hashlib,json,os,subprocess
from pathlib import Path
r=Path(REMOTE);r.mkdir(exist_ok=True)
for n,b in FILES.items():
 p=r/n;p.parent.mkdir(parents=True,exist_ok=True);data=base64.b64decode(b)
 if p.exists():assert p.read_bytes()==data,'Refusing overwrite: '+n
 else:p.write_bytes(data)
resource=dict(load=os.getloadavg(),node=subprocess.check_output(['scontrol','show','node'],text=True),memory=Path('/proc/meminfo').read_text(),queue=subprocess.check_output(['squeue','-o','%.12i %.15j %.10T %.6C %.20R'],text=True),filesystem_free_bytes=os.statvfs(r).f_bavail*os.statvfs(r).f_frsize)
print(json.dumps(resource))
'''.replace('REMOTE',repr(REMOTE)).replace('FILES',repr(files))
 result=remote(code);save(OUT/(stage+'_resource_preflight.json'),result)
 batch='#!/bin/bash\n#SBATCH --partition=debug\n#SBATCH --nodes=1\n#SBATCH --ntasks=1\n#SBATCH --cpus-per-task=1\n#SBATCH --hint=nomultithread\n#SBATCH --time=12:00:00\n#SBATCH --array=0-'+str(len(names)-1)+'%2\n#SBATCH --job-name=P2-WP5-'+stage+'\n#SBATCH --output='+REMOTE+'/'+stage+'-%A_%a.log\nset -eu\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1\nruns=('+ ' '.join(names)+')\nsrun --mpi=none -n 1 -c 1 --cpu-bind=cores python3 '+REMOTE+'/single_wavelength_worker.py '+REMOTE+' "${runs[$SLURM_ARRAY_TASK_ID]}"\n'
 (OUT/(stage+'.sbatch')).write_text(batch,encoding='utf-8');remote('from pathlib import Path\nimport json\np=Path('+repr(REMOTE+'/'+stage+'.sbatch')+')\ns='+repr(batch)+'\nassert not p.exists() or p.read_text()==s\np.write_text(s)\nprint(json.dumps(dict(staged=True)))')
 print('prepared',stage,names)
def submit(stage):
 if stage=='full':assert json.load(open(OUT/'short_validation.json'))['advance_to_full']
 code='''import json,subprocess
from pathlib import Path
r=Path(REMOTE);p=r/('submission_'+STAGE+'.json')
with p.open('x') as f:f.write('{}')
job=subprocess.check_output(['sbatch','--parsable',str(r/(STAGE+'.sbatch'))],text=True).strip().split(';')[0]
rec=dict(job_id=job,stage=STAGE,host='sal',root=str(r));p.write_text(json.dumps(rec));print(json.dumps(rec))
'''.replace('REMOTE',repr(REMOTE)).replace('STAGE',repr(stage))
 result=remote(code);save(OUT/('submission_'+stage+'.json'),result);print(result)
def collect(stage):
 job=json.load(open(OUT/('submission_'+stage+'.json')))['job_id'];names=[p.name for p in (OUT/'runs').iterdir() if p.name.endswith('-'+stage)]
 code='''import json,base64,hashlib,subprocess
from pathlib import Path
r=Path(REMOTE);res=dict(queue=subprocess.check_output(['squeue','-j',JOB,'-o','%.15i %.10T %.10M %.6C %.20R'],text=True),runs={})
for name in NAMES:
 d=r/'runs'/name;info={}
 for n in ['started.json','solver_finished.json','finished.json','failure.json','stop_record.json']:
  if (d/n).exists():info[n]=json.loads((d/n).read_text())
 if (d/'solver_finished.json').exists():
  info['files']={str(p.relative_to(d)):dict(data=base64.b64encode(p.read_bytes()).decode(),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in d.rglob('*') if p.is_file() and p.suffix in ['.csv','.json','.log','.pvd','.npz']}
 res['runs'][name]=info
print(json.dumps(res))
'''.replace('REMOTE',repr(REMOTE)).replace('JOB',repr(job)).replace('NAMES',repr(names))
 result=remote(code)
 for name,info in result['runs'].items():
  for rel,v in info.pop('files',{}).items():
   data=base64.b64decode(v['data']);assert hashlib.sha256(data).hexdigest()==v['sha256'];p=(CACHE/'runs'/name/rel) if rel.endswith('.npz') else (OUT/'runs'/name/rel);p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(data)
 save(OUT/(stage+'_latest_status.json'),result);print(result['queue']);print({k:('finished' if 'finished.json' in v else 'solver_finished' if 'solver_finished.json' in v else 'running' if 'started.json' in v else 'pending') for k,v in result['runs'].items()})
def archive(stage):
 result={}
 for d in sorted((OUT/'runs').iterdir()):
  if not d.name.endswith('-'+stage) or not (d/'finished.json').exists():continue
  record=json.load(open(d/'finished.json'))['archive'];dest=CACHE/(d.name+'_fields.tar.gz')
  if not dest.exists() or hashlib.sha256(dest.read_bytes()).hexdigest()!=record['sha256']:
   subprocess.run(['scp','-o','BatchMode=yes','sal:'+record['path'],str(dest)],check=True,timeout=600)
  assert hashlib.sha256(dest.read_bytes()).hexdigest()==record['sha256']
  expected={Path(r['path']).name:r['sha256'] for r in json.load(open(d/'post/field_hashes.json'))}
  with tarfile.open(dest) as tar:
   for n,h in expected.items():
    with tar.extractfile(n) as stream:assert hashlib.file_digest(stream,'sha256').hexdigest()==h
  result[d.name]=dict(**record,local_path=str(dest),raw_vtu_hashes_verified=len(expected))
 save(OUT/(stage+'_archive.json'),result);print('archived',len(result),'runs')
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('action',choices=['prepare','submit','collect','archive']);ap.add_argument('stage',choices=['short','full']);a=ap.parse_args();globals()[a.action](a.stage)
