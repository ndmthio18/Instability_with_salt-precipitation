"""Submit and collect post-only WP7 spatial analysis on sal."""
import argparse,hashlib,json,subprocess,tempfile,tarfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'analysis/wp7_shortwave_validation_20261103';CACHE=ROOT/'tmp/wp7_shortwave_validation_20261103';REMOTE='/tmp/zhixin_paper2_wp7_shortwave_20261103_r01'
RUNS={'M0':'WP7-R48-M0-4T','M1':'WP7-R48-M1-4T'}
def run(args,**kw):return subprocess.run(args,check=True,text=True,**kw)
def ssh(command,capture=False):
    result=run(['ssh','sal',command],capture_output=capture)
    return result.stdout if capture else ''
def scp(src,dst):run(['scp',str(src),dst])
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
def submit():
    for name in ['wp7_spatial.py','wp6_4t_metrics.py','single_wavelength_post.py','wp7_sal_spatial_job.sh']:scp(ROOT/'scripts'/name,f'sal:{REMOTE}/{name}')
    jobs={}
    for role,name in RUNS.items():
        assert ssh(f'test -f {REMOTE}/runs/{name}/finished.json && echo yes || true',True).strip()=='yes'
        out=f'spatial_{role}';assert not ssh(f'test -e {REMOTE}/{out} && echo yes || true',True).strip()
        jobs[role]=ssh(f"cd {REMOTE} && sbatch --parsable --job-name=WP7-spatial-{role} wp7_sal_spatial_job.sh {REMOTE} {name} {out}",True).strip().split(';')[0]
    save(OUT/'sal_jobs_spatial.json',jobs);print(json.dumps(jobs,indent=2))
def status():
    jobs=json.loads((OUT/'sal_jobs_spatial.json').read_text());ids=','.join(jobs.values());queue=ssh(f"squeue -h -j {ids} -o '%i|%T|%M|%R' || true",True);details={j:ssh(f'scontrol show job {j} 2>/dev/null || true',True) for j in jobs.values()};result=dict(jobs=jobs,queue=queue.splitlines(),details=details);save(OUT/'sal_status_spatial.json',result);print(json.dumps(result,indent=2))
def collect():
    jobs=json.loads((OUT/'sal_jobs_spatial.json').read_text());ids=','.join(jobs.values());assert not ssh(f"squeue -h -j {ids} || true",True).strip()
    receipt={}
    for role in RUNS:
        remote=f'{REMOTE}/spatial_{role}'
        with tempfile.TemporaryDirectory() as td:
            arc=Path(td)/f'spatial_{role}.tar.gz';ssh(f'cd {remote} && tar -czf {REMOTE}/spatial_{role}.tar.gz .');scp(f'sal:{REMOTE}/spatial_{role}.tar.gz',arc)
            dest=OUT/f'spatial_{role}';dest.mkdir(exist_ok=True)
            with tarfile.open(arc) as tar:tar.extractall(dest,filter='data')
        receipt[role]={p.name:dict(sha256=sha(p),bytes=p.stat().st_size) for p in dest.iterdir() if p.is_file()}
        npz=dest/'selected_fields.npz';target=CACHE/f'spatial_{role}_selected_fields.npz';target.write_bytes(npz.read_bytes());npz.unlink();receipt[role]['selected_fields.npz']['local_cache']=str(target)
    save(OUT/'sal_collect_spatial.json',receipt);print(json.dumps(receipt,indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['submit','status','collect']);a=p.parse_args();globals()[a.action]()
