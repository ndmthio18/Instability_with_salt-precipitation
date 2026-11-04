"""One immutable-input WP9 run on hpcc, with threaded assembly.

A close fork of `wp7_worker.py`. Only three things differ, each because the leg is
threaded rather than pinned to a single core:

  * the process keeps the whole CPU set Slurm allocated instead of pinning one CPU;
  * DUMUX_NUM_THREADS is set from the allocation, which is the only control the frozen
    executable exposes for its TBB backend;
  * the stride-aware `wp9_post.py` replaces `wp6_post.py`, because this leg archives a
    frame every second accepted step.

Everything else, including every hash check, every resource gate and the archive step,
is unchanged. Threading was measured to be bit-identical to serial on this executable
(analysis/wp9_r96_convergence_20261104/hpcc_portability/thread_scaling.json).
"""
import hashlib,json,os,re,resource,subprocess,sys,tarfile,time,traceback
from datetime import datetime,timezone
from pathlib import Path


def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def save(p,v):p.write_text(json.dumps(v,indent=2))


def main(root,rid):
    root=Path(root);run=root/'runs'/rid
    p=json.loads((run/'plan.json').read_text());bundle=root/'bundle'
    assert os.environ.get('WP9_MANAGED')=='1' and not (run/'started.json').exists()
    assert os.environ.get('SLURM_JOB_ID'),'the worker requires a Slurm allocation'
    cpus=sorted(os.sched_getaffinity(0))
    threads=int(os.environ.get('WP9_THREADS',len(cpus)))
    assert threads<=len(cpus),'more threads than allocated cpus'
    for n,v in json.loads((bundle/'runtime_manifest.json').read_text()).items():
        assert sha(bundle/n)==v['sha256']
    assert sha(bundle/'bin/paper2_transport')==p['elf_sha256']
    assert sha(run/'params.input')==p['input_sha256']
    assert sha(root/'preregistration.json')==p['preregistration_sha256']
    assert sha(root/'user_authorization.json')==p['user_authorization_sha256']
    avail=int(next(s.split()[1] for s in Path('/proc/meminfo').read_text().splitlines()
                   if s.startswith('MemAvailable:')))*1024
    st=os.statvfs(root);free=st.f_bavail*st.f_frsize
    assert avail>p['minimum_available_memory_bytes'] and free>p['minimum_free_disk_bytes']
    save(run/'started.json',dict(utc=datetime.now(timezone.utc).isoformat(),
        hostname=os.uname().nodename,affinity=cpus,threads=threads,
        available_memory_bytes=avail,free_disk_bytes=free,load=os.getloadavg(),
        elf_sha256=p['elf_sha256'],slurm_job_id=os.environ['SLURM_JOB_ID']))

    def limit():resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,8*1024**3))
    env=dict(os.environ,DUMUX_NUM_THREADS=str(threads))
    start=time.monotonic();reason=None
    with (run/'solver.log').open('w') as log:
        try:
            q=subprocess.run([sys.executable,str(bundle/'launch_solver.py'),'params.input'],
                             cwd=run,stdout=log,stderr=subprocess.STDOUT,
                             timeout=p['hard_wall_limit_s']-600,preexec_fn=limit,env=env)
            rc=q.returncode
        except subprocess.TimeoutExpired:rc=-9;reason='solver_wall_timeout'
    usage=resource.getrusage(resource.RUSAGE_CHILDREN)
    txt=(run/'solver.log').read_text()
    stop=json.loads((run/'stop_record.json').read_text()) if (run/'stop_record.json').exists() else None
    solver=dict(returncode=rc,reason=reason,wall_s=time.monotonic()-start,
                cpu_s=usage.ru_utime+usage.ru_stime,peak_rss_kib=usage.ru_maxrss,stop=stop,
                threads=threads,
                newton_iterations=len(re.findall(r'Newton iteration\s+\d+ done',txt)),
                failure_lines=[s for s in txt.splitlines()
                               if re.search(r'fail|retry|reject|reduc.*time.step',s,re.I)])
    save(run/'solver_finished.json',solver)
    if rc or not stop or stop['reason'] not in ['requested_horizon','surface_solubility_threshold',
                                                'interior_solubility_threshold']:
        return 1
    poststart=time.monotonic()
    q=subprocess.run([sys.executable,str(root/'wp9_post.py'),str(run),str(run/'post')],
                     capture_output=True,text=True,timeout=5400)
    (run/'post.log').write_text(q.stdout+q.stderr)
    if q.returncode:return 2
    archive=root/f'{rid}_fields.tar.gz'
    with tarfile.open(archive,'w:gz',compresslevel=1) as tar:
        for f in run.iterdir():
            if f.suffix in ['.vtu','.pvd']:tar.add(f,arcname=f.name)
    save(run/'finished.json',dict(solver=solver,post_archive_wall_s=time.monotonic()-poststart,
        archive=dict(path=str(archive),bytes=archive.stat().st_size,sha256=sha(archive)),
        files={str(f.relative_to(run)).replace(os.sep,'/'):dict(sha256=sha(f),bytes=f.stat().st_size)
               for f in run.rglob('*') if f.is_file() and f.suffix!='.vtu'}))
    return 0


if __name__=='__main__':
    try:raise SystemExit(main(sys.argv[1],sys.argv[2]))
    except Exception:
        save(Path(sys.argv[1])/'runs'/sys.argv[2]/'failure.json',
             dict(traceback=traceback.format_exc()));raise
