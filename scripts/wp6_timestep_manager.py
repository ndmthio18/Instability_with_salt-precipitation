"""Exactly one compute-host worker, with a process-group hard wall limit."""
import json,os,signal,subprocess,sys,time
from pathlib import Path

def bounded_run(command,log,env,limit_s):
    start=time.monotonic()
    p=subprocess.Popen(command,stdout=log,stderr=subprocess.STDOUT,env=env,start_new_session=True)
    timed_out=False
    try:rc=p.wait(timeout=limit_s)
    except subprocess.TimeoutExpired:
        timed_out=True;os.killpg(p.pid,signal.SIGKILL);rc=p.wait()
    finally:
        # Also remove any remaining child after an early launcher/worker failure.
        try:os.killpg(p.pid,signal.SIGKILL)
        except ProcessLookupError:pass
    return dict(returncode=rc,timed_out=timed_out,wall_s=time.monotonic()-start,pid=p.pid,hard_limit_s=limit_s)

def main(root,stage,names):
    root=Path(root);names=json.loads(names);assert names==['WP6-M1-dt_half']
    cpuinfo=subprocess.check_output(['lscpu','-p=CPU,CORE,SOCKET'],text=True);cores={}
    for line in cpuinfo.splitlines():
        if line.startswith('#'):continue
        cpu,core,socket=map(int,line.split(','));cores.setdefault((socket,core),cpu)
    cpu=list(cores.values())[0];assert os.getloadavg()[0]<8
    (root/(stage+'_manager_resources.json')).write_text(json.dumps(dict(cpus=[cpu],load=os.getloadavg(),lscpu=cpuinfo,hard_wall_limit_s=3600),indent=2))
    name=names[0];env=dict(os.environ,WP6_MANAGED='1',WP6_CPU=str(cpu),OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
    (root/(stage+'_manager_state.json')).write_text(json.dumps(dict(active={name:'bounded worker'},waiting=[],finished={})))
    with (root/'runs'/name/'worker.log').open('w') as log:
        result=bounded_run([sys.executable,str(root/'wp6_worker.py'),str(root),name],log,env,3600)
    (root/(stage+'_wall_limit_result.json')).write_text(json.dumps(result,indent=2))
    done={name:result['returncode']}
    if result['timed_out']:(root/'runs'/name/'failure.json').write_text(json.dumps(dict(reason='hard_wall_timeout',**result)))
    (root/(stage+'_manager_state.json')).write_text(json.dumps(dict(active={},waiting=[],finished=done)))
    (root/(stage+'_manager_finished.json')).write_text(json.dumps(done))
if __name__=='__main__':main(*sys.argv[1:])
