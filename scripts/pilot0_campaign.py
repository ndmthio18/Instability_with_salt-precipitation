"""Prepare, queue, inspect and collect the one authorized Pilot-0 baseline M1.

Run from Windows. No retries or additional simulations are implicit.
"""
import argparse
import base64
import hashlib
import json
from pathlib import Path
import subprocess
from datetime import datetime, timezone
from prepare_transport_run import parse

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'analysis/pilot0_integrated_20261102'
REMOTE='/temp/zhixin/paper2_transport_20261102/pilot0_integrated_20261102_r01'
M0='/temp/zhixin/paper2_transport_20261102/unperturbed_36877_20261102_r01'
BRIDGE='/mnt/d/3OutputFile/1Research/Article/22Surrogate_model/paper/2/scripts/frodo_bridge.py'


def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(path,obj):path.write_text(json.dumps(obj,indent=2,ensure_ascii=False)+'\n',encoding='utf-8')
def remote(code):
    p=subprocess.run(['wsl','python3',BRIDGE,'run','python3 -'],input=code.encode(),capture_output=True,timeout=150)
    if p.returncode:raise RuntimeError(p.stderr.decode(errors='replace')+'\n'+p.stdout.decode(errors='replace'))
    return json.loads(p.stdout)


def prepare():
    planned=ROOT/'analysis/matched_pair_36877_20261102/planned_M1'
    plan=json.loads((planned/'plan.json').read_text())
    assert sha(planned/'params.input')==plan['input_sha256']
    assert sha(planned/'eigenprofile.dat')==plan['profile_sha256']
    assert all(sha(ROOT/'model/transport_v1/source'/n)==h for n,h in plan['source_sha256'].items())
    actual=ROOT/'analysis/matched_pair_36877_20261102/M0_actual_params.input'
    v0=parse(actual.read_text());v1=parse((planned/'params.input').read_text())
    allowed={'Problem.Name','Problem.PerturbationType','Problem.Amplitude','Problem.StandardDeviationSaltMoleFraction','Problem.EigenprofilePath'}
    differing=[k for k in set(v0)|set(v1) if v0.get(k)!=v1.get(k)]
    assert set(differing)<=allowed
    assert v1==plan['parameters']
    m0_plan=json.loads((ROOT/'analysis/unperturbed_36877_20261102/plan.json').read_text())
    plan.update(frozen_preparation_M0_input_sha256=plan['M0_input_sha256'],
        frozen_preparation_M0_plan_sha256=plan['M0_plan_sha256'],
        M0_input_sha256=sha(actual),
        M0_plan_semantic_sha256=hashlib.sha256(json.dumps(m0_plan,sort_keys=True,separators=(',',':')).encode()).hexdigest())
    plan.update(status='registered_authorized_by_command3',solver_executed=False,
        authorization_file='command/3.md',authorization_sha256=sha(ROOT/'command/3.md'),
        source_paper2_commit=subprocess.check_output(['git','-C',str(ROOT),'rev-parse','HEAD'],text=True).strip(),
        frozen_preparation_plan_sha256=sha(planned/'plan.json'),
        wall_timeout_s=43200,queue_timeout_s=43200,
        resource_gate=dict(max_load_1m=20,max_load_5m=22,min_available_memory_GiB=4,min_free_disk_GiB=10,no_other_paper2_solver=True),
        M0_remote=M0,remote_source='/temp/zhixin/paper2_transport_20261102/source',
        launch_gate='command/3.md authorized; worker verifies resource and source gates',
        additional_runs_authorized_by_this_worker=False)
    staged={'M1/params.input':(planned/'params.input').read_bytes(),'M1/eigenprofile.dat':(planned/'eigenprofile.dat').read_bytes(),
        'M1/plan.json':(json.dumps(plan,indent=2)+'\n').encode(),
        'pilot0_remote_worker.py':(ROOT/'scripts/pilot0_remote_worker.py').read_bytes(),
        'pilot0_field_review.py':(ROOT/'scripts/pilot0_field_review.py').read_bytes()}
    payload={name:base64.b64encode(data).decode() for name,data in staged.items()}
    code="""import base64,hashlib,json,shutil
from pathlib import Path
r=Path(REMOTE);m0=Path(M0);payload=PAYLOAD
plan=json.loads(base64.b64decode(payload['M1/plan.json']))
old=json.loads((m0/'plan.json').read_text());finished=json.loads((m0/'finished.json').read_text())
assert finished['returncode']==0 and finished['stop']['reason']=='requested_horizon'
assert hashlib.sha256((m0/'params.input').read_bytes()).hexdigest()==old['input_sha256']==plan['M0_input_sha256']
assert hashlib.sha256(json.dumps(old,sort_keys=True,separators=(',',':')).encode()).hexdigest()==plan['M0_plan_semantic_sha256']
assert old['executable']['sha256']==plan['executable']['sha256']
assert old['source_sha256']==plan['source_sha256']
assert hashlib.sha256(Path(plan['executable']['path']).read_bytes()).hexdigest()==plan['executable']['sha256']
assert all(hashlib.sha256((Path(plan['remote_source'])/n).read_bytes()).hexdigest()==h for n,h in plan['source_sha256'].items())
r.mkdir(exist_ok=False)
for name,b in payload.items():
 p=r/name;p.parent.mkdir(exist_ok=True);p.write_bytes(base64.b64decode(b))
small=['params.input','plan.json','finished.json','stop_record.json','transport_budget.csv','surface_events.csv','transport_metadata.json']
small += [p.name for p in m0.glob('*_sigmatop.csv')]
original={}
(r/'M0').mkdir()
for name in small:
 b=(m0/name).read_bytes()
 if name in finished['files']:assert hashlib.sha256(b).hexdigest()==finished['files'][name]['sha256']
 (r/'M0'/name).write_bytes(b);original[name]=base64.b64encode(b).decode()
result=dict(status='prepared_not_started',remote=str(r),input_differences=DIFFERING,
 M0_completed_and_verified=True,ELF_and_21_source_hashes_verified=True,
 M0_actual_plan_sha256=hashlib.sha256((m0/'plan.json').read_bytes()).hexdigest(),
 staged_sha256={n:hashlib.sha256((r/n).read_bytes()).hexdigest() for n in payload},M0_files=original)
print(json.dumps(result))
""".replace('REMOTE',repr(REMOTE)).replace('M0)',repr(M0)+')').replace('PAYLOAD',repr(payload)).replace('DIFFERING',repr(sorted(differing)))
    # Replace only the explicit source placeholder; never shell-interpolate file content.
    result=remote(code)
    for name,data in staged.items():
        p=OUT/'runs'/name if name.startswith('M1/') else OUT/'execution'/name
        p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(data)
    for name,encoded in result.pop('M0_files').items():
        p=OUT/'runs/M0'/name;p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(base64.b64decode(encoded))
    result['checked_utc']=datetime.now(timezone.utc).isoformat()
    write(OUT/'preflight.json',result);print(json.dumps(result,indent=2))


def submit():
    result=remote("""import json,subprocess,sys
from pathlib import Path
r=Path(REMOTE)
assert (r/'M1/plan.json').exists()
with (r/'launch_claim.json').open('x') as claim:claim.write('{}')
with (r/'worker.log').open('w') as log:
 p=subprocess.Popen([sys.executable,str(r/'pilot0_remote_worker.py'),str(r)],stdin=subprocess.DEVNULL,stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
record=dict(worker_pid=p.pid,remote=str(r),status='resource_gated_worker_submitted_not_solver_launch')
(r/'launch.json').write_text(json.dumps(record,indent=2));print(json.dumps(record))
""".replace('REMOTE',repr(REMOTE)))
    write(OUT/'launch.json',result);print(json.dumps(result,indent=2))


def collect():
    result=remote("""import base64,hashlib,json
from pathlib import Path
r=Path(REMOTE);result={}
for name in ['worker_state.json','worker_failure.json','resource_at_launch.json','M1/started.json','M1/initial_field_check.json','M1/finished.json']:
 p=r/name
 if p.exists():result[name]=json.loads(p.read_text())
if (r/'M1/finished.json').exists():
 files={}
 for p in (r/'M1').iterdir():
  if p.is_file() and p.suffix in ('.json','.csv','.input','.dat','.pvd','.log'):
   b=p.read_bytes();files[p.name]=dict(base64=base64.b64encode(b).decode(),sha256=hashlib.sha256(b).hexdigest())
 result['files']=files
print(json.dumps(result))
""".replace('REMOTE',repr(REMOTE)))
    for name,rec in result.pop('files',{}).items():
        data=base64.b64decode(rec['base64']);assert hashlib.sha256(data).hexdigest()==rec['sha256']
        (OUT/'runs/M1'/name).write_bytes(data)
    write(OUT/'latest_status.json',result)
    for name in ('M1/initial_field_check.json','M1/finished.json'):
        if name in result:write(OUT/Path(name).name,result[name])
    state=result.get('worker_state.json',{})
    print(json.dumps({k:v for k,v in state.items() if k!='finished'},indent=2))
    if result.get('worker_failure.json'):print(json.dumps(result['worker_failure.json'],indent=2))


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('action',choices=['prepare','submit','collect']);args=ap.parse_args()
    OUT.mkdir(exist_ok=True)
    {'prepare':prepare,'submit':submit,'collect':collect}[args.action]()
