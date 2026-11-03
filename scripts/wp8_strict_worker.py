"""Bounded immutable-executable background worker; never retries or extends."""
import csv,hashlib,json,os,re,resource,signal,subprocess,sys,tarfile,time,traceback
from pathlib import Path
import xml.etree.ElementTree as ET
import numpy as np

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,v):p.write_text(json.dumps(v,indent=2)+'\n')
def widths(x):return np.diff(np.r_[x[0],(x[:-1]+x[1:])/2,x[-1]])
def field(p):
    raw=p.read_bytes();aa={v.attrib['Name']:np.fromstring(v.text or '',sep=' ') for v in ET.fromstring(raw).iter('DataArray') if v.attrib.get('Name') in ['Coordinates','moleFraction','sw','pw']}
    xyz=aa['Coordinates'].reshape(-1,3);xs=np.unique(xyz[:,0]);zs=np.unique(xyz[:,1]);assert len(xs)==2 and len(zs)==201
    order=np.argsort(np.searchsorted(zs,xyz[:,1])*2+np.searchsorted(xs,xyz[:,0]));x,sw,pw=[aa[k][order].reshape(len(zs),2) for k in ['moleFraction','sw','pw']]
    xr=float(np.max(abs(x[:,1]-x[:,0])));pr=float(np.max(abs(pw[:,1]-pw[:,0])))
    return dict(x_mean=x.mean(axis=1),sw_mean=sw.mean(axis=1),sx_mean=(sw*x).mean(axis=1),pw_mean=pw.mean(axis=1),z_m=zs,x_range=xr,p_range=pr,sha256=hashlib.sha256(raw).hexdigest(),bytes=len(raw),max_x=float(x.max()))
def main():
    run=Path(sys.argv[1]);plan=json.loads((run/'plan.json').read_text());exe=Path(plan['executable']);assert sha(exe)==plan['elf_sha256'];assert sha(run/'params.input')==plan['input_sha256']
    src=Path(plan['remote_source'])
    for n,h in plan['source_hashes'].items():assert sha(src/n)==h,n
    start=time.time();cpu=plan['cpu'];os.sched_setaffinity(0,{cpu});data=[];seen=set();trigger=None;proc=None
    def bind():resource.setrlimit(resource.RLIMIT_AS,(2*1024**3,2*1024**3));os.sched_setaffinity(0,{cpu})
    def scan(final=False):
        for f in sorted(run.glob('*.vtu')):
            if f.name in seen:continue
            try:d=field(f)
            except (ET.ParseError,ValueError):
                if final:raise
                break
            d['name']=f.name;data.append(d);seen.add(f.name)
            if d['x_range']!=0 or d['p_range']!=0:raise RuntimeError('lost exact horizontal equality')
    save(run/'started.json',dict(pid=os.getpid(),cpu=cpu,started_unix=start,plan_sha256=sha(run/'plan.json')))
    try:
        with (run/'solver.log').open('w') as log:
            proc=subprocess.Popen([str(exe),'params.input'],cwd=run,stdout=log,stderr=subprocess.STDOUT,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1',MKL_NUM_THREADS='1'),preexec_fn=bind,start_new_session=True)
            while proc.poll() is None:
                scan()
                if time.time()-start>plan['wall_cap_s']:trigger='wall_limit';break
                save(run/'monitor.json',dict(frames_checked=len(data),elapsed_s=time.time()-start,max_x=max((d['max_x'] for d in data),default=None),exact_uniform=True))
                time.sleep(.5)
            if trigger:os.killpg(proc.pid,signal.SIGTERM);proc.wait(timeout=10)
        scan(final=True)
    except Exception:
        trigger=traceback.format_exc()
        if proc and proc.poll() is None:os.killpg(proc.pid,signal.SIGTERM);proc.wait(timeout=10)
    solver_wall=time.time()-start
    post=run/'post';post.mkdir(exist_ok=True)
    if data:
        np.savez_compressed(post/'profiles.npz',z_m=data[0]['z_m'],**{k:np.array([d[k] for d in data]) for k in ['x_mean','sw_mean','sx_mean','pw_mean']})
        save(post/'field_manifest.json',[{k:d[k] for k in ['name','sha256','bytes','x_range','p_range','max_x']} for d in data])
    archive=run/'fields.tar.gz'
    with tarfile.open(archive,'w:gz') as tar:
        for f in sorted(run.glob('*.vtu')):tar.add(f,arcname=f.name)
    result=dict(returncode=proc.returncode if proc else None,trigger=trigger,solver_and_live_audit_wall_s=solver_wall,total_wall_s=time.time()-start,cpu=cpu,peak_child_rss_kib=resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss,frames=len(data),raw_field_bytes=sum(d['bytes'] for d in data),archive=dict(path=str(archive),sha256=sha(archive),bytes=archive.stat().st_size),archive_budget_pass=archive.stat().st_size<=250*1024**2,elf_sha256=sha(exe),input_sha256=sha(run/'params.input'),worker_sha256=sha(Path(__file__)))
    save(run/'finished.json',result)
if __name__=='__main__':main()
