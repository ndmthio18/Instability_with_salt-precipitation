"""One Slurm-owned, 12-hour-limited leg on sal; parallel M1 requires explicit authorization.

The solver runs with the original frodo loader/libraries, one rank and one core.
No retries or changes of physical parameters are implemented here.
"""
import csv
import hashlib
import json
import os
from pathlib import Path
import resource
import subprocess
import sys
import time
import traceback
from datetime import datetime,timezone


def now():return datetime.now(timezone.utc).isoformat()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(p,v):
    tmp=p.with_suffix(p.suffix+'.tmp');tmp.write_text(json.dumps(v,indent=2)+'\n');tmp.replace(p)


def check_background(reference,current):
    metrics={};checks={}
    for name,keys in [('transport_budget.csv',['inventory_salt_m','cumulative_cut_salt_m','salt_residual_m']),
                      ('surface_events.csv',['x_top_area_mean','x_top_max','x_domain_max','sigma_top_node'])]:
        with (reference/name).open() as f:a=list(csv.DictReader(f))
        with (current/name).open() as f:b=list(csv.DictReader(f))
        checks[name+'/row_count']=len(a)==len(b)
        if len(a)!=len(b):continue
        checks[name+'/time_and_region']=all(abs(float(x['time_s'])-float(y['time_s']))<1e-8 and x.get('region')==y.get('region') for x,y in zip(a,b))
        for k in keys:
            diff=max(abs(float(x[k])-float(y[k])) for x,y in zip(a,b));metrics[name+'/'+k]=diff
            tolerance=1e-14 if name=='transport_budget.csv' else (1e-11 if k!='sigma_top_node' else 1e-11)
            checks[name+'/'+k]=diff<=tolerance
    return dict(passed=all(checks.values()),checks=checks,max_abs_differences=metrics,
        scope='Host/runtime bridge only; final physical pair uses two new sal legs. These thresholds are not transport-effectiveness thresholds.')


def check_m1_launch(root, run, plan):
    authorization = root/'parallel_authorization.json'
    if authorization.exists():
        auth = json.loads(authorization.read_text())
        assert auth['allow_parallel_pair'] is True
        assert str(auth['M1_job_id']) == os.environ['SLURM_JOB_ID']
        assert os.environ.get('SLURM_CPUS_PER_TASK') == '1'
        m0 = root/'M0'
        start = json.loads((m0/'started.json').read_text())
        assert str(start['SLURM_JOB_ID']) == str(auth['M0_job_id'])
        assert start['elf_sha256'] == plan['executable']['sha256']
        assert set(os.sched_getaffinity(0)).isdisjoint(start['affinity']), 'M0/M1 CPU affinity overlaps'
        stop = json.loads((m0/'stop_record.json').read_text())
        assert stop['initial_diagnostics_written'] and stop['accepted_steps'] >= 1
        assert list(m0.glob('*-00000.vtu'))
        for name in ('wrapper_failure.json', 'wrapper_exception.json'):
            assert not (m0/name).exists(), 'M0 failure: '+name
        if (m0/'finished.json').exists():
            done = json.loads((m0/'finished.json').read_text())
            assert done['returncode'] == 0 and not done.get('wrapper_stop_reason')
            assert json.loads((m0/'bridge_validation.json').read_text())['passed']
    else:
        assert json.loads((root/'M0/bridge_validation.json').read_text())['passed']
        assert json.loads((root/'M0/finished.json').read_text())['returncode'] == 0
    assert sha(run/'eigenprofile.dat') == plan['profile_sha256']


def main(root,role):
    root=root.resolve();run=root/role;bundle=root/'bundle';p=json.loads((run/'plan.json').read_text())
    started=run/'started.json';state=run/'runtime_state.json'
    if started.exists():raise RuntimeError('Already started; duplicate forbidden')
    if not os.environ.get('SLURM_JOB_ID'):raise RuntimeError('Must run inside Slurm allocation')
    manifest=json.loads((bundle/'runtime_manifest.json').read_text())
    for name,rec in manifest.items():
        assert sha(bundle/name)==rec['sha256'],name
    assert sha(run/'params.input')==p['input_sha256']
    assert sha(bundle/'bin/paper2_transport')==p['executable']['sha256']
    if role=='M1':
        check_m1_launch(root, run, p)
    if os.statvfs(root).f_bavail*os.statvfs(root).f_frsize<4*1024**3:raise RuntimeError('Insufficient scratch')
    command=[sys.executable,str(bundle/'launch_solver.py'),'params.input']
    # This process is already under srun CPU placement. Strip scheduler-discovery
    # variables only from the one-rank MPI child, retaining cgroup/CPU ownership.
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
    def limits():resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,8*1024**3))
    before=time.monotonic();proc=None;reason=None;initial_checked=False
    try:
        with (run/'solver.log').open('w') as log:
            proc=subprocess.Popen(command,cwd=run,stdout=log,stderr=subprocess.STDOUT,env=env,preexec_fn=limits)
            write(started,dict(utc=now(),pid=proc.pid,command=command,SLURM_JOB_ID=os.environ['SLURM_JOB_ID'],
                execution_host='sal',elf_sha256=sha(bundle/'bin/paper2_transport'),single_rank=True,
                solver_address_space_limit_bytes=8*1024**3,affinity=sorted(os.sched_getaffinity(0))))
            while proc.poll() is None:
                stop=json.loads((run/'stop_record.json').read_text()) if (run/'stop_record.json').exists() else None
                if role=='M1' and stop and stop['accepted_steps']>=1 and not initial_checked:
                    check=subprocess.run([sys.executable,str(root/'pilot0_field_review.py'),'--m0',str(root/'M0'),
                        '--m1',str(run),'--out',str(run/'initial_field_check.json')],env=env,timeout=120)
                    initial_checked=True
                    if check.returncode:reason='failed_initial_field_gate';break
                write(state,dict(status='in_progress',role=role,utc=now(),elapsed_s=time.monotonic()-before,stop=stop))
                if time.monotonic()-before>=43200:reason='wall_timeout';break
                time.sleep(20)
    except Exception:
        reason='wrapper_exception';write(run/'wrapper_exception.json',dict(traceback=traceback.format_exc()))
    finally:
        if proc and proc.poll() is None:
            proc.terminate()
            try:proc.wait(timeout=30)
            except subprocess.TimeoutExpired:proc.kill();proc.wait()
    rc=proc.returncode if proc else None
    stop=json.loads((run/'stop_record.json').read_text()) if (run/'stop_record.json').exists() else None
    success=rc==0 and not reason and stop and stop['reason'] in ('requested_horizon','surface_solubility_threshold','interior_solubility_threshold')
    if role=='M0' and success:
        bridge=check_background(bundle/'M0',run);write(run/'bridge_validation.json',bridge)
        success=bridge['passed']
        if not success:reason='failed_host_bridge'
    if role=='M1' and success and not initial_checked:
        reason='initial_check_missing';success=False
    result=dict(returncode=rc,wrapper_stop_reason=reason,wall_s=time.monotonic()-before,stop=stop,
        finished_utc=now(),execution_host='sal',files={})
    for f in run.iterdir():
        if f.is_file() and f.name not in ('runtime_state.json','finished.json'):
            result['files'][f.name]=dict(bytes=f.stat().st_size,sha256=sha(f))
    write(run/'finished.json',result)
    write(state,dict(status='completed' if success else 'failed',role=role,utc=now(),finished=result))
    return 0 if success else 1


if __name__=='__main__':
    try:raise SystemExit(main(Path(sys.argv[1]),sys.argv[2]))
    except Exception:
        write(Path(sys.argv[1])/sys.argv[2]/'wrapper_failure.json',dict(status='failed',traceback=traceback.format_exc(),utc=now()))
        raise
