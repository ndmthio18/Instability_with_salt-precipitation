"""WP7 bounded one-core workers with measured physical-core placement."""
import concurrent.futures,json,os,subprocess,sys,time
from pathlib import Path
from wp6_timestep_manager import bounded_run
def cores(n):
    groups={}
    for line in subprocess.check_output(['lscpu','-p=CPU,CORE,SOCKET'],text=True).splitlines():
        if not line.startswith('#'):
            cpu,core,socket=map(int,line.split(','));groups.setdefault((socket,core),[]).append(cpu)
    def snap():return {s.split()[0]:list(map(int,s.split()[1:])) for s in Path('/proc/stat').read_text().splitlines() if s.startswith('cpu')}
    a=snap();time.sleep(1);b=snap();busy={}
    for k in b:
        d=[v-u for u,v in zip(a[k],b[k])];busy[k]=1-(d[3]+d[4])/sum(d) if sum(d) else 0
    q=sorted((max(busy['cpu'+str(i)] for i in v),v[0],v) for v in groups.values());free=[x for x in q if x[0]<.25]
    assert len(free)>=n and os.getloadavg()[0]<8
    return [x[1] for x in free[:n]],dict(choices=q,busy=busy,load=os.getloadavg(),lscpu=subprocess.check_output(['lscpu','-p=CPU,CORE,SOCKET'],text=True))
def main(root,stage,names):
    root=Path(root);names=json.loads(names);assert stage in ['preflight','r48pair'];assert len(names) in [1,2]
    cpus,res=cores(len(names));res.update(cpus=cpus,stage=stage);(root/f'{stage}_manager_resources.json').write_text(json.dumps(res,indent=2))
    def one(name,cpu):
        plan=json.loads((root/'runs'/name/'plan.json').read_text());limit=plan['hard_wall_limit_s']
        env=dict(os.environ,WP6_MANAGED='1',WP6_CPU=str(cpu),OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
        with (root/'runs'/name/'worker.log').open('w') as log:return bounded_run([sys.executable,str(root/'wp7_worker.py'),str(root),name],log,env,limit)
    result={};done={}
    with concurrent.futures.ThreadPoolExecutor(max_workers=len(names)) as pool:
        pending={pool.submit(one,n,c):n for n,c in zip(names,cpus)}
        while pending:
            (root/f'{stage}_manager_state.json').write_text(json.dumps(dict(active=list(pending.values()),finished=done)))
            ready,_=concurrent.futures.wait(pending,timeout=2,return_when=concurrent.futures.FIRST_COMPLETED)
            for f in ready:
                n=pending.pop(f);result[n]=f.result();done[n]=result[n]['returncode'];(root/f'{stage}_wall_limit_result.json').write_text(json.dumps(result,indent=2))
    (root/f'{stage}_manager_state.json').write_text(json.dumps(dict(active=[],finished=done)));(root/f'{stage}_manager_finished.json').write_text(json.dumps(done))
if __name__=='__main__':main(*sys.argv[1:])
