"""Stage, submit, collect, and verify the bounded WP7 preflight on sal."""
import argparse,hashlib,json,shutil,subprocess,tarfile,tempfile
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'analysis/wp7_shortwave_validation_20261103'
CACHE=ROOT/'tmp/wp7_shortwave_validation_20261103'
REMOTE='/tmp/zhixin_paper2_wp7_shortwave_20261103_r01'
RUNS=['WP7-R48-M1-bench10','WP7-1D-background-T']
R48_RUNS=['WP7-R48-M0-4T','WP7-R48-M1-4T']

def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()
def run(args,**kwargs): return subprocess.run(args,check=True,text=True,**kwargs)
def ssh(command,capture=False):
    q=run(['ssh','sal',command],capture_output=capture)
    return q.stdout if capture else ''
def scp(source,destination): run(['scp',str(source),destination])
def save(path,value): path.write_text(json.dumps(value,indent=2)+'\n')

def stage():
    bundle=CACHE/'bundle.tar.gz'
    assert sha(bundle)=='84240d783b68b9e051c7a31af8da6c20c53978c856c7919fead686794b0ae58a'
    amendment=OUT/'resource_amendment_E02.json'; amend_sha=sha(amendment)
    for name in RUNS:
        plan_path=OUT/'runs'/name/'plan.json'; plan=json.loads(plan_path.read_text())
        plan['minimum_free_disk_bytes']=20*1024**3
        plan['resource_amendment_sha256']=amend_sha
        save(plan_path,plan)
    ssh(f"test ! -e {REMOTE} && mkdir -p {REMOTE}/runs")
    scp(bundle,f'sal:{REMOTE}/bundle.tar.gz')
    ssh(f"cd {REMOTE} && tar -xzf bundle.tar.gz && test -x bundle/bin/paper2_transport")
    for filename in ['wp7_worker.py','wp7_sal_worker.py','wp7_sal_job.sh','wp6_post.py','single_wavelength_post.py']:
        scp(ROOT/'scripts'/filename,f'sal:{REMOTE}/{filename}')
    for filename in ['preregistration.json','resource_amendment_E02.json']:
        scp(OUT/filename,f'sal:{REMOTE}/{filename}')
    for name in RUNS:
        ssh(f'mkdir -p {REMOTE}/runs/{name}')
        for filename in ['params.input','plan.json','eigenprofile.dat']:
            scp(OUT/'runs'/name/filename,f'sal:{REMOTE}/runs/{name}/{filename}')
    probe=json.loads(ssh(f"python3 - <<'PY'\nimport hashlib,json,os\nfrom pathlib import Path\nr=Path('{REMOTE}'); p=r/'bundle/bin/paper2_transport'\nprint(json.dumps(dict(root=str(r),elf_sha256=hashlib.sha256(p.read_bytes()).hexdigest(),free_disk_bytes=os.statvfs(r).f_bavail*os.statvfs(r).f_frsize)))\nPY",True))
    assert probe['elf_sha256']=='112076759fdff0942000307f729c2bd3c7e96963229b43c1b67c1f12bbf9df81'
    probe.update(bundle_archive_sha256=sha(bundle),resource_amendment_sha256=amend_sha,runs=RUNS)
    save(OUT/'sal_stage_preflight.json',probe);print(json.dumps(probe,indent=2))

def submit():
    assert (OUT/'sal_stage_preflight.json').exists()
    jobs={}
    for name in RUNS:
        raw=ssh(f"cd {REMOTE} && sbatch --parsable wp7_sal_job.sh {REMOTE} {name}",True).strip()
        jobs[name]=raw.split(';')[0]
    save(OUT/'sal_jobs_preflight.json',jobs);print(json.dumps(jobs,indent=2))

def repair():
    first=json.loads((OUT/'sal_jobs_preflight.json').read_text())
    relocation=OUT/'runtime_relocation_E03.json'; relocation_sha=sha(relocation)
    wrapper=ROOT/'scripts/orted_relocatable.sh'
    for name,job in first.items():
        source=OUT/'runs'/name; failed=OUT/'failed_attempts'/f'{name}_job{job}'
        failed.parent.mkdir(parents=True,exist_ok=True)
        assert (source/'solver_finished.json').exists() and not failed.exists()
        solver=json.loads((source/'solver_finished.json').read_text())
        assert solver['returncode']==1 and solver['newton_iterations']==0 and solver['stop'] is None
        shutil.move(source,failed);source.mkdir(parents=True)
        for filename in ['params.input','plan.json','eigenprofile.dat']:
            shutil.copy2(failed/filename,source/filename)
        plan=json.loads((source/'plan.json').read_text());plan['runtime_relocation_sha256']=relocation_sha;save(source/'plan.json',plan)
        scp(f'sal:{REMOTE}/slurm-{job}.out',failed/f'slurm-{job}.out')
    scp(wrapper,f'sal:{REMOTE}/orted.relocatable')
    scp(relocation,f'sal:{REMOTE}/runtime_relocation_E03.json')
    remote_names=json.dumps(first)
    command=f"""python3 - <<'PY'
import hashlib,json,os,shutil
from pathlib import Path
r=Path('{REMOTE}'); jobs=json.loads({remote_names!r}); wrapper=r/'orted.relocatable'; target=r/'bundle/bin/orted'
old=hashlib.sha256(target.read_bytes()).hexdigest(); real=hashlib.sha256((r/'bundle/bin/orted.real').read_bytes()).hexdigest(); elf=hashlib.sha256((r/'bundle/bin/paper2_transport').read_bytes()).hexdigest()
for name,job in jobs.items():
 src=r/'runs'/name; dst=r/'runs'/(name+'_attempt1_failed_job'+job)
 assert not dst.exists(); assert json.loads((src/'solver_finished.json').read_text())['newton_iterations']==0
 src.rename(dst);src.mkdir()
shutil.copy2(wrapper,target);target.chmod(0o755)
manifest=json.loads((r/'bundle/runtime_manifest.json').read_text());manifest['bin/orted']['sha256']=hashlib.sha256(target.read_bytes()).hexdigest();(r/'bundle/runtime_manifest.json').write_text(json.dumps(manifest,indent=2)+'\\n')
receipt=dict(old_orted_sha256=old,new_orted_sha256=manifest['bin/orted']['sha256'],orted_real_sha256=real,elf_sha256=elf,failed_jobs=jobs)
(r/'runtime_relocation_receipt.json').write_text(json.dumps(receipt,indent=2)+'\\n');print(json.dumps(receipt))
PY"""
    receipt=json.loads(ssh(command,True));assert receipt['elf_sha256']=='112076759fdff0942000307f729c2bd3c7e96963229b43c1b67c1f12bbf9df81'
    for name in RUNS:
        for filename in ['params.input','plan.json','eigenprofile.dat']:
            scp(OUT/'runs'/name/filename,f'sal:{REMOTE}/runs/{name}/{filename}')
    save(OUT/'runtime_relocation_receipt.json',receipt);print(json.dumps(receipt,indent=2))

def submit_retry():
    assert (OUT/'runtime_relocation_receipt.json').exists();jobs={}
    for name in RUNS:
        raw=ssh(f"cd {REMOTE} && sbatch --parsable wp7_sal_job.sh {REMOTE} {name}",True).strip();jobs[name]=raw.split(';')[0]
    save(OUT/'sal_jobs_preflight_retry.json',jobs);print(json.dumps(jobs,indent=2))

def status_retry():
    jobs=json.loads((OUT/'sal_jobs_preflight_retry.json').read_text());ids=','.join(jobs.values())
    queue=ssh(f"squeue -h -j {ids} -o '%i|%T|%M|%R' || true",True)
    details={job:ssh(f'scontrol show job {job} 2>/dev/null || true',True) for job in jobs.values()}
    result=dict(jobs=jobs,queue=queue.splitlines(),details=details);save(OUT/'sal_status_preflight_retry.json',result);print(json.dumps(result,indent=2))

def collect_retry():
    retry=OUT/'sal_jobs_preflight_retry.json'; original=OUT/'sal_jobs_preflight.json'; backup=OUT/'sal_jobs_preflight_initial.json'
    shutil.copy2(original,backup);shutil.copy2(retry,original)
    try:collect()
    finally:shutil.copy2(backup,original)
    shutil.move(OUT/'sal_collect_preflight.json',OUT/'sal_collect_preflight_retry.json')

def stage_r48():
    authorization=OUT/'user_authorization_E05.json';auth_sha=sha(authorization)
    relocation_sha=sha(OUT/'runtime_relocation_E03.json')
    probe=json.loads(ssh(f"python3 - <<'PY'\nimport json,os\nfrom pathlib import Path\nr=Path('{REMOTE}');print(json.dumps(dict(free_disk_bytes=os.statvfs(r).f_bavail*os.statvfs(r).f_frsize,available_memory_bytes=int(next(x.split()[1] for x in Path('/proc/meminfo').read_text().splitlines() if x.startswith('MemAvailable:')))*1024)))\nPY",True))
    assert probe['free_disk_bytes']>60*1024**3 and probe['available_memory_bytes']>16*1024**3
    scp(authorization,f'sal:{REMOTE}/user_authorization_E05.json')
    for name in R48_RUNS:
        local=OUT/'runs'/name;plan=json.loads((local/'plan.json').read_text())
        plan['minimum_free_disk_bytes']=60*1024**3;plan['runtime_relocation_sha256']=relocation_sha;plan['user_authorization_sha256']=auth_sha
        save(local/'plan.json',plan)
        exists=ssh(f"test -e {REMOTE}/runs/{name} && echo yes || true",True).strip();assert not exists,name
        ssh(f'mkdir -p {REMOTE}/runs/{name}')
        for filename in ['params.input','plan.json','eigenprofile.dat']:
            scp(local/filename,f'sal:{REMOTE}/runs/{name}/{filename}')
    probe.update(runs=R48_RUNS,user_authorization_sha256=auth_sha,runtime_relocation_sha256=relocation_sha)
    save(OUT/'sal_stage_r48.json',probe);print(json.dumps(probe,indent=2))

def submit_r48():
    assert (OUT/'sal_stage_r48.json').exists();jobs={}
    for name in R48_RUNS:
        raw=ssh(f"cd {REMOTE} && sbatch --parsable --time=02:00:00 --job-name={name} wp7_sal_job.sh {REMOTE} {name}",True).strip();jobs[name]=raw.split(';')[0]
    save(OUT/'sal_jobs_r48.json',jobs);print(json.dumps(jobs,indent=2))

def status_r48():
    jobs=json.loads((OUT/'sal_jobs_r48.json').read_text());ids=','.join(jobs.values())
    queue=ssh(f"squeue -h -j {ids} -o '%i|%T|%M|%R' || true",True)
    details={job:ssh(f'scontrol show job {job} 2>/dev/null || true',True) for job in jobs.values()}
    progress={}
    for name in R48_RUNS:
        raw=ssh(f"python3 - <<'PY'\nimport json\nfrom pathlib import Path\nd=Path('{REMOTE}/runs/{name}');out={{}}\nfor n in ['started.json','solver_finished.json','finished.json','failure.json','stop_record.json']:\n p=d/n\n if p.exists():\n  try:out[n]=json.loads(p.read_text())\n  except Exception:out[n]={{'invalid':True}}\nprint(json.dumps(out))\nPY",True)
        progress[name]=json.loads(raw)
    result=dict(jobs=jobs,queue=queue.splitlines(),details=details,progress=progress);save(OUT/'sal_status_r48.json',result);print(json.dumps(result,indent=2))

def collect_r48():
    jobs=json.loads((OUT/'sal_jobs_r48.json').read_text());receipt_path=OUT/'sal_collect_r48.json'
    receipts=load_json(receipt_path) if receipt_path.exists() else {'runs':{}}
    for name in R48_RUNS:
        if receipts['runs'].get(name,{}).get('complete'):continue
        remote_run=f'{REMOTE}/runs/{name}'
        state=ssh(f"test -f {remote_run}/finished.json && echo finished || (test -f {remote_run}/solver_finished.json && echo solver_finished || echo active)",True).strip()
        if state=='active':receipts['runs'][name]={'complete':False,'state':'active'};continue
        with tempfile.TemporaryDirectory() as td:
            meta=Path(td)/f'{name}_metadata.tar.gz';ssh(f"cd {remote_run} && tar --exclude='*.vtu' -czf {REMOTE}/{name}_metadata.tar.gz .");scp(f'sal:{REMOTE}/{name}_metadata.tar.gz',meta)
            local=OUT/'runs'/name
            with tarfile.open(meta) as tar:tar.extractall(local,filter='data')
        if not (local/'finished.json').exists():receipts['runs'][name]={'complete':False,'state':'failed','solver':json.loads((local/'solver_finished.json').read_text())};continue
        archive_info=json.loads((local/'finished.json').read_text())['archive'];target=CACHE/f'{name}_fields.tar.gz';scp(f"sal:{archive_info['path']}",target);assert sha(target)==archive_info['sha256']
        expected={Path(x['path']).name:x['sha256'] for x in json.loads((local/'post/field_hashes.json').read_text())};seen=set()
        with tarfile.open(target,'r|gz') as tar:
            for item in tar:
                if item.name in expected:
                    f=tar.extractfile(item);assert f and hashlib.file_digest(f,'sha256').hexdigest()==expected[item.name];seen.add(item.name)
        assert seen==set(expected);receipts['runs'][name]=dict(complete=True,state='archived',archive_sha256=sha(target),verified_vtu=len(seen),local_archive=str(target))
        save(receipt_path,receipts)
    save(receipt_path,receipts);print(json.dumps(receipts,indent=2))

def load_json(path):return json.loads(path.read_text())

def status():
    jobs=json.loads((OUT/'sal_jobs_preflight.json').read_text())
    ids=','.join(jobs.values())
    queue=ssh(f"squeue -h -j {ids} -o '%i|%T|%M|%R' || true",True)
    accounting=ssh(f"sacct -n -P -j {ids} --format=JobID,JobName,State,Elapsed,TotalCPU,MaxRSS,ExitCode || true",True)
    result=dict(jobs=jobs,queue=queue.splitlines(),accounting=accounting.splitlines())
    save(OUT/'sal_status_preflight.json',result);print(json.dumps(result,indent=2))

def collect():
    jobs=json.loads((OUT/'sal_jobs_preflight.json').read_text()); ids=','.join(jobs.values())
    queue=ssh(f"squeue -h -j {ids} -o '%i|%T' || true",True).strip()
    if queue: raise RuntimeError('preflight jobs are still active: '+queue)
    accounting=ssh(f"sacct -n -P -j {ids} --format=JobID,JobName,State,Elapsed,TotalCPU,MaxRSS,ExitCode 2>/dev/null || true",True)
    receipts={}
    for name in RUNS:
        remote_run=f'{REMOTE}/runs/{name}'
        with tempfile.TemporaryDirectory() as td:
            meta=Path(td)/f'{name}_metadata.tar.gz'
            ssh(f"cd {remote_run} && tar --exclude='*.vtu' -czf {REMOTE}/{name}_metadata.tar.gz .")
            scp(f'sal:{REMOTE}/{name}_metadata.tar.gz',meta)
            local=OUT/'runs'/name;local.mkdir(parents=True,exist_ok=True)
            with tarfile.open(meta) as tar:tar.extractall(local,filter='data')
        finished=local/'finished.json'
        if not finished.exists():
            receipts[name]=dict(complete=False,failure=(local/'failure.json').read_text() if (local/'failure.json').exists() else None)
            continue
        archive_info=json.loads(finished.read_text())['archive']; target=CACHE/f'{name}_fields.tar.gz'
        scp(f"sal:{archive_info['path']}",target)
        assert sha(target)==archive_info['sha256']
        expected={Path(x['path']).name:x['sha256'] for x in json.loads((local/'post/field_hashes.json').read_text())}
        seen=set()
        with tarfile.open(target,'r|gz') as tar:
            for item in tar:
                if item.name in expected:
                    f=tar.extractfile(item);assert f and hashlib.file_digest(f,'sha256').hexdigest()==expected[item.name];seen.add(item.name)
        assert seen==set(expected)
        receipts[name]=dict(complete=True,archive_sha256=sha(target),verified_vtu=len(seen),local_archive=str(target))
    save(OUT/'sal_collect_preflight.json',dict(accounting=accounting.splitlines(),runs=receipts));print(json.dumps(receipts,indent=2))

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('action',choices=['stage','submit','status','collect','repair','submit_retry','status_retry','collect_retry','stage_r48','submit_r48','status_r48','collect_r48']);args=parser.parse_args()
    globals()[args.action]()
