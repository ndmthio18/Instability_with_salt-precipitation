"""Two independent single-core legs with existing process-group wall guards."""
import concurrent.futures,json,os,subprocess,sys,time
from pathlib import Path
from wp6_timestep_manager import bounded_run

def select_cores():
    cpuinfo=subprocess.check_output(['lscpu','-p=CPU,CORE,SOCKET'],text=True);groups={}
    for line in cpuinfo.splitlines():
        if line.startswith('#'):continue
        cpu,core,socket=map(int,line.split(','));groups.setdefault((socket,core),[]).append(cpu)
    def sample():return {s.split()[0]:list(map(int,s.split()[1:])) for s in Path('/proc/stat').read_text().splitlines() if s.startswith('cpu')}
    a=sample();time.sleep(1);b=sample();busy={}
    for k in b:
        d=[v-u for u,v in zip(a[k],b[k])];busy[k]=1-(d[3]+d[4])/sum(d) if sum(d) else 0
    choices=sorted((max(busy['cpu'+str(i)] for i in v),v[0],v) for v in groups.values())
    free=[x for x in choices if x[0]<.25];assert len(free)>=2 and os.getloadavg()[0]<8, 'no two idle physical cores'
    return [x[1] for x in free[:2]],dict(physical_core_choices=choices,cpu_busy_fraction=busy,lscpu=cpuinfo,load=os.getloadavg())

def main(root,stage,names):
    root=Path(root);names=json.loads(names);assert names==['WP6-M0-explore4T','WP6-M1-explore4T']
    cpus,resource=select_cores();resource.update(cpus=cpus,hard_wall_limit_s_per_leg=3600,total_allocated_core_hours_cap=2)
    (root/(stage+'_manager_resources.json')).write_text(json.dumps(resource,indent=2))
    def run(name,cpu):
        env=dict(os.environ,WP6_MANAGED='1',WP6_CPU=str(cpu),OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
        with (root/'runs'/name/'worker.log').open('w') as log:result=bounded_run([sys.executable,str(root/'wp6_worker.py'),str(root),name],log,env,3600)
        if result['timed_out']:(root/'runs'/name/'failure.json').write_text(json.dumps(dict(reason='hard_wall_timeout',**result)))
        return result
    results={};done={}
    with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
        pending={pool.submit(run,name,cpu):name for name,cpu in zip(names,cpus)}
        while pending:
            (root/(stage+'_manager_state.json')).write_text(json.dumps(dict(active={n:'bounded worker' for n in pending.values()},waiting=[],finished=done)))
            ready,_=concurrent.futures.wait(pending,timeout=2,return_when=concurrent.futures.FIRST_COMPLETED)
            for f in ready:
                n=pending.pop(f);results[n]=f.result();done[n]=results[n]['returncode']
                (root/(stage+'_wall_limit_result.json')).write_text(json.dumps(results,indent=2))
    (root/(stage+'_manager_state.json')).write_text(json.dumps(dict(active={},waiting=[],finished=done)))
    (root/(stage+'_manager_finished.json')).write_text(json.dumps(done))
if __name__=='__main__':main(*sys.argv[1:])
