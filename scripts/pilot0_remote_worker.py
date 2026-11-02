"""Bounded single-run worker on frodo: resource queue, immutable inputs, monitoring.

Only registered M1 is authorized here. No automatic retries, extensions or batches.
"""
import fcntl
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time
import traceback
from datetime import datetime, timezone


def now():return datetime.now(timezone.utc).isoformat()
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()


def resources(root):
    busy=[]
    for p in Path('/proc').iterdir():
        if not p.name.isdigit():continue
        try:
            exe=os.readlink(p/'exe')
            if Path(exe).name=='paper2_transport':
                busy.append(dict(pid=int(p.name),exe=exe,cwd=os.readlink(p/'cwd')))
        except (OSError,PermissionError):pass
    memory={s.split(':')[0]:int(s.split()[1])*1024 for s in Path('/proc/meminfo').read_text().splitlines()}
    return dict(utc=now(),load_1m=os.getloadavg()[0],load_5m=os.getloadavg()[1],cpu=os.cpu_count(),
        available_memory_bytes=memory['MemAvailable'],free_disk_bytes=shutil.disk_usage(root).free,active_paper2_solvers=busy)


def eligible(r):
    return not r['active_paper2_solvers'] and r['load_1m']<=20 and r['load_5m']<=22 and r['available_memory_bytes']>=4*1024**3 and r['free_disk_bytes']>=10*1024**3


def atomic(path,value):
    temp=path.with_suffix(path.suffix+'.tmp');temp.write_text(json.dumps(value,indent=2)+'\n');temp.replace(path)


def main(root):
    root=root.resolve(); run=root/'M1'; p=json.loads((run/'plan.json').read_text())
    state=root/'worker_state.json'; start=time.monotonic()
    lock=(root.parent/'.paper2_pilot0.lock').open('a')
    try:fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    except BlockingIOError:raise RuntimeError('Another Pilot-0 worker owns the launch lock')
    if (run/'started.json').exists():raise RuntimeError('Registered run already started; no duplicate execution')
    while True:
        r=resources(root)
        if eligible(r):break
        elapsed=time.monotonic()-start
        atomic(state,dict(status='queued_resource_gate',elapsed_queue_s=elapsed,resources=r,
            reason='Wait for existing Paper2 solver to finish and load/memory/disk gates to pass'))
        if elapsed>=p['queue_timeout_s']:
            atomic(state,dict(status='pending_resource_queue_expired',resources=r,utc=now()));return
        time.sleep(60)
    atomic(root/'resource_at_launch.json',r)
    assert sha(run/'params.input')==p['input_sha256']
    assert sha(run/'eigenprofile.dat')==p['profile_sha256']
    exe=Path(p['executable']['path']);assert sha(exe)==p['executable']['sha256']
    source=Path(p['remote_source'])
    assert all(sha(source/n)==h for n,h in p['source_sha256'].items())
    # Preserve the exact ELF used; subsequent build changes cannot replace it during this run.
    artifacts=run/'artifacts';artifacts.mkdir(exist_ok=False)
    executable=artifacts/'paper2_transport';shutil.copy2(exe,executable)
    assert sha(executable)==p['executable']['sha256']
    env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
    command=['nice','-n','10',str(executable),'params.input']
    with (run/'solver.log').open('w') as log:
        proc=subprocess.Popen(command,cwd=run,stdout=log,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,env=env,start_new_session=True)
        launched=time.monotonic()
        atomic(run/'started.json',dict(utc=now(),pid=proc.pid,command=command,
            executed_elf_sha256=sha(executable),threads={k:env[k] for k in ('OPENBLAS_NUM_THREADS','OMP_NUM_THREADS','MKL_NUM_THREADS')}))
        initial_checked=False;reason=None
        while proc.poll() is None:
            elapsed=time.monotonic()-launched
            # The first accepted step guarantees t=0 VTK and diagnostics have finished writing.
            stopfile=run/'stop_record.json'
            stop=json.loads(stopfile.read_text()) if stopfile.exists() else None
            if stop and stop['accepted_steps']>=1 and not initial_checked:
                check=subprocess.run([sys.executable,str(root/'pilot0_field_review.py'),'--m0',p['M0_remote'],
                    '--m1',str(run),'--out',str(run/'initial_field_check.json')],env=env,timeout=120)
                initial_checked=True
                if check.returncode!=0:reason='failed_initial_field_gate';break
            atomic(state,dict(status='in_progress',utc=now(),elapsed_solver_s=elapsed,solver_pid=proc.pid,
                stop=stop,initial_field_checked=initial_checked))
            if elapsed>=p['wall_timeout_s']:reason='wall_timeout';break
            time.sleep(30)
        if reason:
            proc.terminate()
            try:proc.wait(timeout=30)
            except subprocess.TimeoutExpired:proc.kill();proc.wait()
        rc=proc.wait()
    stop=json.loads((run/'stop_record.json').read_text()) if (run/'stop_record.json').exists() else None
    result=dict(returncode=rc,wall_s=time.monotonic()-launched,finished_utc=now(),stop=stop,
        wrapper_stop_reason=reason,files={})
    for f in run.iterdir():
        if f.is_file():result['files'][f.name]=dict(bytes=f.stat().st_size,sha256=sha(f))
    atomic(run/'finished.json',result)
    complete=rc==0 and not reason and initial_checked and stop and stop['reason'] in ('requested_horizon','surface_solubility_threshold','interior_solubility_threshold')
    atomic(state,dict(status='completed_pending_pair_analysis' if complete else 'failed',utc=now(),
        finished=result,initial_field_checked=initial_checked))


if __name__=='__main__':
    root=Path(sys.argv[1])
    try:main(root)
    except Exception:
        atomic(root/'worker_failure.json',dict(status='failed',utc=now(),traceback=traceback.format_exc()))
        raise
