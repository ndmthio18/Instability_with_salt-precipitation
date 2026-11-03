"""One immutable-input WP7 run; postprocess and archive accepted fields."""
import hashlib,json,os,re,resource,subprocess,sys,tarfile,time,traceback
from datetime import datetime,timezone
from pathlib import Path
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,v):p.write_text(json.dumps(v,indent=2))
def main(root,rid):
    root=Path(root);run=root/'runs'/rid;p=json.loads((run/'plan.json').read_text());bundle=root/'bundle'
    assert os.environ.get('WP6_MANAGED')=='1' and not (run/'started.json').exists();os.sched_setaffinity(0,{int(os.environ['WP6_CPU'])})
    for n,v in json.loads((bundle/'runtime_manifest.json').read_text()).items():assert sha(bundle/n)==v['sha256']
    assert sha(bundle/'bin/paper2_transport')==p['elf_sha256'] and sha(run/'params.input')==p['input_sha256'] and sha(root/'preregistration.json')==p['preregistration_sha256']
    avail=int(next(s.split()[1] for s in Path('/proc/meminfo').read_text().splitlines() if s.startswith('MemAvailable:')))*1024;free=os.statvfs(root).f_bavail*os.statvfs(root).f_frsize
    min_free=p.get('minimum_free_disk_bytes',200*1024**3)
    assert avail>16*1024**3 and free>min_free
    save(run/'started.json',dict(utc=datetime.now(timezone.utc).isoformat(),affinity=sorted(os.sched_getaffinity(0)),available_memory_bytes=avail,free_disk_bytes=free,minimum_free_disk_bytes=min_free,load=os.getloadavg(),elf_sha256=p['elf_sha256']))
    def limit():resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,8*1024**3))
    start=time.monotonic();reason=None
    with (run/'solver.log').open('w') as log:
        try:q=subprocess.run([sys.executable,str(bundle/'launch_solver.py'),'params.input'],cwd=run,stdout=log,stderr=subprocess.STDOUT,timeout=p['hard_wall_limit_s']-300,preexec_fn=limit);rc=q.returncode
        except subprocess.TimeoutExpired:rc=-9;reason='solver_wall_timeout'
    usage=resource.getrusage(resource.RUSAGE_CHILDREN);txt=(run/'solver.log').read_text();stop=json.loads((run/'stop_record.json').read_text()) if (run/'stop_record.json').exists() else None
    solver=dict(returncode=rc,reason=reason,wall_s=time.monotonic()-start,cpu_s=usage.ru_utime+usage.ru_stime,peak_rss_kib=usage.ru_maxrss,stop=stop,newton_iterations=len(re.findall(r'Newton iteration\s+\d+ done',txt)),failure_lines=[s for s in txt.splitlines() if re.search(r'fail|retry|reject|reduc.*time.step',s,re.I)])
    save(run/'solver_finished.json',solver)
    if rc or not stop or stop['reason'] not in ['requested_horizon','surface_solubility_threshold','interior_solubility_threshold']:return 1
    poststart=time.monotonic();q=subprocess.run([sys.executable,str(root/'wp6_post.py'),str(run),str(run/'post')],capture_output=True,text=True,timeout=1800);(run/'post.log').write_text(q.stdout+q.stderr)
    if q.returncode:return 2
    archive=root/f'{rid}_fields.tar.gz'
    with tarfile.open(archive,'w:gz',compresslevel=1) as tar:
        for f in run.iterdir():
            if f.suffix in ['.vtu','.pvd']:tar.add(f,arcname=f.name)
    save(run/'finished.json',dict(solver=solver,post_archive_wall_s=time.monotonic()-poststart,archive=dict(path=str(archive),bytes=archive.stat().st_size,sha256=sha(archive)),files={str(f.relative_to(run)):dict(sha256=sha(f),bytes=f.stat().st_size) for f in run.rglob('*') if f.is_file() and f.suffix!='.vtu'}))
    return 0
if __name__=='__main__':
    try:raise SystemExit(main(sys.argv[1],sys.argv[2]))
    except Exception:
        save(Path(sys.argv[1])/'runs'/sys.argv[2]/'failure.json',dict(traceback=traceback.format_exc()));raise
