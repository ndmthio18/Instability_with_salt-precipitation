"""Managed two-physical-core execution on the authorized frodo compute workstation."""
import os,json,subprocess,sys,time
from pathlib import Path
root=Path(sys.argv[1]);stage=sys.argv[2];names=json.loads(sys.argv[3]);cpuinfo=subprocess.check_output(['lscpu','-p=CPU,CORE,SOCKET'],text=True);cores={}
for line in cpuinfo.splitlines():
 if line.startswith('#'):continue
 cpu,core,socket=map(int,line.split(','));cores.setdefault((socket,core),cpu)
cpus=list(cores.values())[:2];assert len(cpus)==2 and os.getloadavg()[0]<8
(root/(stage+'_manager_resources.json')).write_text(json.dumps(dict(cpus=cpus,load=os.getloadavg(),lscpu=cpuinfo),indent=2))
active={};done={};todo=list(names)
while todo or active:
 for cpu in cpus:
  if cpu not in active and todo:
   name=todo.pop(0);log=(root/'runs'/name/'worker.log').open('w');env=dict(os.environ,WP6_MANAGED='1',WP6_CPU=str(cpu),OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1');p=subprocess.Popen([sys.executable,str(root/'wp6_worker.py'),str(root),name],stdout=log,stderr=subprocess.STDOUT,env=env);active[cpu]=(p,name,log)
 for cpu,(p,name,log) in list(active.items()):
  if p.poll() is not None:
   done[name]=p.returncode;log.close();del active[cpu]
   if p.returncode:
    for pending in todo:done[pending]='not_launched_after_failure'
    todo=[]
 (root/(stage+'_manager_state.json')).write_text(json.dumps(dict(active={name:p.pid for p,name,log in active.values()},waiting=todo,finished=done)))
 time.sleep(2)
(root/(stage+'_manager_finished.json')).write_text(json.dumps(done))
