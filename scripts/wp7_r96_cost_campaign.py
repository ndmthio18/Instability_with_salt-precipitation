"""Conditionally execute only the preregistered R96 ten-step cost benchmark."""
import argparse,hashlib,json,subprocess,tarfile,tempfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'analysis/wp7_shortwave_validation_20261103';CACHE=ROOT/'tmp/wp7_shortwave_validation_20261103';REMOTE='/tmp/zhixin_paper2_wp7_shortwave_20261103_r01';NAME='WP7-R96-M1-bench10'
def run(args,**kw):return subprocess.run(args,check=True,text=True,**kw)
def ssh(cmd,capture=False):
    result=run(['ssh','sal',cmd],capture_output=capture)
    return result.stdout if capture else ''
def scp(src,dst):run(['scp',str(src),dst])
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
def triggered():
    d=json.loads((OUT/'r48_comparison.json').read_text())['decision'];return d['near_resolution_limit'] or not d['spectral_convergence']
def stage():
    assert triggered();local=OUT/'runs'/NAME;plan=json.loads((local/'plan.json').read_text());plan['runtime_relocation_sha256']=sha(OUT/'runtime_relocation_E03.json');save(local/'plan.json',plan)
    assert not ssh(f'test -e {REMOTE}/runs/{NAME} && echo yes || true',True).strip();ssh(f'mkdir -p {REMOTE}/runs/{NAME}')
    for filename in ['params.input','plan.json','eigenprofile.dat']:scp(local/filename,f'sal:{REMOTE}/runs/{NAME}/{filename}')
    rec=dict(triggered=True,trigger_decision=json.loads((OUT/'r48_comparison.json').read_text())['decision'],input_sha256=plan['input_sha256'],trigger_sha256=plan['trigger_sha256']);save(OUT/'sal_stage_r96_cost.json',rec);print(json.dumps(rec,indent=2))
def submit():
    job=ssh(f"cd {REMOTE} && sbatch --parsable --time=00:10:00 --job-name=WP7-R96-cost wp7_sal_job.sh {REMOTE} {NAME}",True).strip().split(';')[0];save(OUT/'sal_job_r96_cost.json',{'job':job});print(job)
def status():
    job=json.loads((OUT/'sal_job_r96_cost.json').read_text())['job'];result=dict(job=job,queue=ssh(f"squeue -h -j {job} -o '%i|%T|%M|%R' || true",True).splitlines(),details=ssh(f'scontrol show job {job} 2>/dev/null || true',True));save(OUT/'sal_status_r96_cost.json',result);print(json.dumps(result,indent=2))
def collect():
    job=json.loads((OUT/'sal_job_r96_cost.json').read_text())['job'];assert not ssh(f'squeue -h -j {job} || true',True).strip();remote=f'{REMOTE}/runs/{NAME}'
    with tempfile.TemporaryDirectory() as td:
        arc=Path(td)/'meta.tar.gz';ssh(f"cd {remote} && tar --exclude='*.vtu' -czf {REMOTE}/{NAME}_metadata.tar.gz .");scp(f'sal:{REMOTE}/{NAME}_metadata.tar.gz',arc);local=OUT/'runs'/NAME
        with tarfile.open(arc) as tar:tar.extractall(local,filter='data')
    finished=json.loads((local/'finished.json').read_text());info=finished['archive'];target=CACHE/f'{NAME}_fields.tar.gz';scp(f"sal:{info['path']}",target);assert sha(target)==info['sha256'];expected={Path(x['path']).name:x['sha256'] for x in json.loads((local/'post/field_hashes.json').read_text())};seen=set()
    with tarfile.open(target,'r|gz') as tar:
        for item in tar:
            if item.name in expected:
                f=tar.extractfile(item);assert f and hashlib.file_digest(f,'sha256').hexdigest()==expected[item.name];seen.add(item.name)
    assert seen==set(expected);receipt=dict(complete=True,archive_sha256=sha(target),verified_vtu=len(seen),local_archive=str(target));save(OUT/'sal_collect_r96_cost.json',receipt);print(json.dumps(receipt,indent=2))
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['stage','submit','status','collect']);a=p.parse_args();globals()[a.action]()
