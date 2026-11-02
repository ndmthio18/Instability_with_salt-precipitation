"""Slurm staging/monitoring for the authorized sal migration. Never resubmit implicitly."""
import argparse
import base64
import hashlib
import json
from pathlib import Path
import subprocess
import tarfile

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'analysis/pilot0_integrated_20261102/sal_migration'
REMOTE='/tmp/zhixin_paper2_pilot0_20261102_r01'


def remote(code):
    p=subprocess.run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=8','sal','python3 -'],input=code.encode(),capture_output=True,timeout=120)
    if p.returncode:raise RuntimeError(p.stderr.decode(errors='replace')+p.stdout.decode(errors='replace'))
    return json.loads(p.stdout)


def save(p,r):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(r,indent=2)+'\n',encoding='utf-8')


def prepare():
    support=hashlib.sha256((ROOT/'tmp/sal_migration/runtime_support.tar.gz').read_bytes()).hexdigest()
    code='''import hashlib,json,os,tarfile
from pathlib import Path
r=Path(REMOTE);b=r/'bundle'
assert hashlib.sha256((r/'runtime_support.tar.gz').read_bytes()).hexdigest()==SUPPORT
with tarfile.open(r/'runtime_support.tar.gz') as tar:tar.extractall(b,filter='data')
launcher=(b/'launch_solver.py').read_text().replace("OMPI_MCA_plm='isolated',","OMPI_MCA_plm='isolated',OMPI_MCA_btl='self',PMIX_MCA_mca_base_component_path=str(r/'pmix'),")
(b/'launch_solver.py').write_text(launcher)
plans={};inputs={}
for role in ('M0','M1'):
 d=r/role;d.mkdir(exist_ok=False)
 p=json.loads((b/role/'plan.json').read_text())
 p.update(status='registered_authorized_sal',role=role,execution_host='sal',wall_timeout_s=43200,
   original_frodo_plan_sha256=hashlib.sha256((b/role/'plan.json').read_bytes()).hexdigest(),
   executable=dict(p['executable'],path=str(b/'bin/paper2_transport')),
   remote_source=str(b/'source'),M0_remote=str(r/'M0'),
   launch_gate='Slurm 1 rank/1 CPU allocation; M1 after successful verified M0',
   runtime_launcher=str(b/'launch_solver.py'),runtime_manifest=str(b/'runtime_manifest.json'),
   comparison_family='new_same_host_sal_M0_M1',solver_executed=False)
 for name in ('params.input','eigenprofile.dat'):
  src=b/role/name
  if src.exists():(d/name).write_bytes(src.read_bytes())
 assert hashlib.sha256((d/'params.input').read_bytes()).hexdigest()==p['input_sha256']
 (d/'plan.json').write_text(json.dumps(p,indent=2)+'\\n');plans[role]=p
 inputs[role]={q.name:base64.b64encode(q.read_bytes()).decode() for q in d.iterdir() if q.is_file()}
 batch='#!/bin/bash\\n#SBATCH --partition=debug\\n#SBATCH --nodes=1\\n#SBATCH --ntasks=1\\n#SBATCH --cpus-per-task=1\\n#SBATCH --time=12:00:00\\n#SBATCH --job-name=p2-'+role+'\\n#SBATCH --output='+str(d/'slurm-%j.log')+'\\nset -eu\\nexport OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1\\nsrun --mpi=none -n 1 -c 1 python3 '+str(r/'pilot0_sal_worker.py')+' '+str(r)+' '+role+'\\n'
 (r/(role+'.sbatch')).write_text(batch)
result=dict(root=str(r),status='prepared_sal_not_submitted',same_elf=plans['M0']['executable']['sha256']==plans['M1']['executable']['sha256'],
 inputs=inputs,scripts={n:(r/n).read_text() for n in ('M0.sbatch','M1.sbatch','pilot0_sal_worker.py')},launcher=launcher,
 support_sha256=SUPPORT,memory_policy='Slurm RealMemory=1M is a site configuration limitation; no memory request; 8GiB address-space limit per solver',
 storage_policy='Isolated /tmp tmpfs; staged inputs backed up locally, outputs must be collected and archived after completion')
print(json.dumps(result))
'''.replace('REMOTE',repr(REMOTE)).replace('SUPPORT',repr(support)).replace('import hashlib,json,os,tarfile','import base64,hashlib,json,os,tarfile')
    result=remote(code)
    for role,files in result.pop('inputs').items():
        for name,b in files.items():
            p=OUT/'runs'/role/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(base64.b64decode(b))
    for n,s in result.pop('scripts').items():
        p=OUT/'execution'/n;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(s.encode())
    (OUT/'execution/launch_solver.py').write_bytes(result.pop('launcher').encode())
    save(OUT/'registration.json',result);print(json.dumps(result,indent=2))


def submit():
    result=remote('''import json,subprocess
from pathlib import Path
r=Path(REMOTE)
with (r/'slurm_submission_claim.json').open('x') as f:f.write('{}')
m0=subprocess.check_output(['sbatch','--parsable',str(r/'M0.sbatch')],text=True).strip().split(';')[0]
record=dict(host='sal',remote=str(r),M0_job_id=m0)
(r/'slurm_jobs.json').write_text(json.dumps(record,indent=2))
m1=subprocess.check_output(['sbatch','--parsable','--dependency=afterok:'+m0,'--kill-on-invalid-dep=yes',str(r/'M1.sbatch')],text=True).strip().split(';')[0]
record['M1_job_id']=m1;record['policy']='M0 then M1, one large solver at a time on sal, each <=12h'
(r/'slurm_jobs.json').write_text(json.dumps(record,indent=2));print(json.dumps(record))
'''.replace('REMOTE',repr(REMOTE)))
    save(OUT/'slurm_jobs.json',result);print(json.dumps(result,indent=2))


def collect():
    result=remote('''import base64,hashlib,json,subprocess
from pathlib import Path
r=Path(REMOTE);jobs=json.loads((r/'slurm_jobs.json').read_text());ids=','.join(str(jobs[k]) for k in ('M0_job_id','M1_job_id'))
result=dict(jobs=jobs,squeue=subprocess.check_output(['squeue','-j',ids,'-o','%.18i %.10j %.9T %.10M %.20R'],text=True),runs={})
accounting=subprocess.run(['sacct','-j',ids,'--format=JobID,State,ExitCode,Elapsed','-n','-P'],capture_output=True,text=True)
result['accounting']=dict(returncode=accounting.returncode,stdout=accounting.stdout,stderr=accounting.stderr)
for role in ('M0','M1'):
 d=r/role;rec={}
 for name in ('started.json','runtime_state.json','finished.json','wrapper_failure.json','wrapper_exception.json','initial_field_check.json','bridge_validation.json'):
  if (d/name).exists():rec[name]=json.loads((d/name).read_text())
 if (d/'stop_record.json').exists():rec['stop']=json.loads((d/'stop_record.json').read_text())
 if (d/'finished.json').exists():
  rec['files']={p.name:dict(base64=base64.b64encode(p.read_bytes()).decode(),sha256=hashlib.sha256(p.read_bytes()).hexdigest()) for p in d.iterdir() if p.is_file() and p.suffix in ('.json','.input','.dat','.csv','.log','.pvd')}
 result['runs'][role]=rec
print(json.dumps(result))
'''.replace('REMOTE',repr(REMOTE)))
    for role,rec in result['runs'].items():
        for name,v in rec.pop('files',{}).items():
            data=base64.b64decode(v['base64']);assert hashlib.sha256(data).hexdigest()==v['sha256']
            (OUT/'runs'/role/name).write_bytes(data)
    save(OUT/'latest_status.json',result)
    print(result['squeue'])
    for role,rec in result['runs'].items():print(role,json.dumps(dict(stop=rec.get('stop'),failure=rec.get('wrapper_failure.json'),finished='finished.json' in rec)))


def archive():
    """Back up completed raw fields from volatile scratch; never delete remote files."""
    result=remote("""import hashlib,json,tarfile
from pathlib import Path
r=Path(REMOTE);result={}
for role in ('M0','M1'):
 d=r/role
 if not (d/'finished.json').exists():continue
 names=sorted(p.name for p in d.iterdir() if p.is_file() and p.suffix in ('.vtu','.pvd'))
 if not names:continue
 target=r/(role+'_fields.tar.gz')
 if not target.exists():
  temp=target.with_suffix('.partial')
  with tarfile.open(temp,'w:gz') as tar:
   for name in names:tar.add(d/name,arcname=name)
  temp.replace(target)
 result[role]=dict(remote_path=str(target),bytes=target.stat().st_size,sha256=hashlib.sha256(target.read_bytes()).hexdigest(),files={n:hashlib.sha256((d/n).read_bytes()).hexdigest() for n in names})
print(json.dumps(result))
""".replace('REMOTE',repr(REMOTE)))
    destination=ROOT/'tmp/sal_migration';destination.mkdir(parents=True,exist_ok=True)
    for role,rec in result.items():
        path=destination/(role+'_fields.tar.gz')
        if not path.exists() or hashlib.sha256(path.read_bytes()).hexdigest()!=rec['sha256']:
            partial=path.with_suffix('.partial')
            subprocess.run(['scp','-o','BatchMode=yes','sal:'+rec['remote_path'],str(partial)],check=True,timeout=600)
            if hashlib.sha256(partial.read_bytes()).hexdigest()!=rec['sha256']:
                raise ValueError('Transferred archive hash mismatch')
            partial.replace(path)
        with tarfile.open(path,'r:gz') as tar:
            assert set(tar.getnames())==set(rec['files'])
            for name,digest in rec['files'].items():
                with tar.extractfile(name) as stream:
                    if hashlib.file_digest(stream,'sha256').hexdigest()!=digest:
                        raise ValueError('Archived member hash mismatch: '+name)
        rec['local_path']=str(path);rec['verified']=True
    if result:save(OUT/'field_archive.json',result)
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    ap=argparse.ArgumentParser();ap.add_argument('action',choices=['prepare','submit','collect','archive']);args=ap.parse_args()
    OUT.mkdir(exist_ok=True);{'prepare':prepare,'submit':submit,'collect':collect,'archive':archive}[args.action]()
