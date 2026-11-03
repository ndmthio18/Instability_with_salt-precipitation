"""One bounded Slurm leg. Independent new directory; original bundle is read-only."""
import hashlib,json,os,re,resource,subprocess,sys,time,tarfile,traceback
from pathlib import Path
from datetime import datetime,timezone

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,v):p.write_text(json.dumps(v,indent=2))
def main(root,rid):
 root=Path(root);run=root/'runs'/rid;p=json.loads((run/'plan.json').read_text());bundle=root/'bundle';assert os.environ.get('WP6_MANAGED')=='1';os.sched_setaffinity(0,{int(os.environ['WP6_CPU'])});assert not (run/'started.json').exists()
 for n,v in json.loads((bundle/'runtime_manifest.json').read_text()).items():assert sha(bundle/n)==v['sha256'],n
 assert sha(bundle/'bin/paper2_transport')==p['elf_sha256'];assert sha(run/'params.input')==p['input_sha256'];assert sha(run/'eigenprofile.dat')==p['profile_sha256'];assert sha(root/'preregistration.json')==p['preregistration_sha256']
 avail=int(next(s.split()[1] for s in Path('/proc/meminfo').read_text().splitlines() if s.startswith('MemAvailable:')))*1024
 assert avail>16*1024**3 and os.statvfs(root).f_bavail*os.statvfs(root).f_frsize>200*1024**3
 save(run/'started.json',dict(utc=datetime.now(timezone.utc).isoformat(),job='frodo-managed-'+rid,affinity=sorted(os.sched_getaffinity(0)),load=os.getloadavg(),available_memory_bytes=avail,elf_sha256=p['elf_sha256'],postprocessor_sha256=sha(root/'wp6_post.py'),launcher_sha256=sha(bundle/'launch_solver.py')))
 start=time.monotonic();reason=None
 def limits():resource.setrlimit(resource.RLIMIT_AS,(8*1024**3,8*1024**3))
 with (run/'solver.log').open('w') as log:
  try:proc=subprocess.run([sys.executable,str(bundle/'launch_solver.py'),'params.input'],cwd=run,stdout=log,stderr=subprocess.STDOUT,timeout=3500,preexec_fn=limits);rc=proc.returncode
  except subprocess.TimeoutExpired:rc=-9;reason='wall_timeout'
 wall=time.monotonic()-start;usage=resource.getrusage(resource.RUSAGE_CHILDREN);txt=(run/'solver.log').read_text();stop=json.loads((run/'stop_record.json').read_text()) if (run/'stop_record.json').exists() else None
 solver=dict(returncode=rc,reason=reason,wall_s=wall,cpu_s=usage.ru_utime+usage.ru_stime,peak_rss_kib=usage.ru_maxrss,stop=stop,newton_iterations=len(re.findall(r'Newton iteration\s+\d+ done',txt)),linear_iterations=None,linear_iterations_status='not emitted by unchanged executable',failure_or_retry_log_lines=[s for s in txt.splitlines() if re.search(r'fail|reduc.*time.step|retry|reject',s,re.I)],solver_output_bytes=sum(f.stat().st_size for f in run.iterdir() if f.is_file()))
 save(run/'solver_finished.json',solver)
 if rc or not stop or stop['reason'] not in ['requested_horizon','surface_solubility_threshold','interior_solubility_threshold']:return 1
 assert stop['last_accepted_time_s']<=p['T_s']+1e-7
 if stop['reason']=='requested_horizon':assert abs(stop['last_accepted_time_s']-p['T_s'])<1e-7
 poststart=time.monotonic();post=subprocess.run([sys.executable,str(root/'wp6_post.py'),str(run),str(run/'post')],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,timeout=1800);(run/'post.log').write_text(post.stdout)
 if post.returncode:raise RuntimeError('postprocessing failed; original solver outputs preserved')
 target=root/(rid+'_fields.tar.gz')
 with tarfile.open(target,'w:gz',compresslevel=1) as tar:
  for f in run.iterdir():
   if f.suffix in ('.vtu','.pvd'):tar.add(f,arcname=f.name)
 save(run/'finished.json',dict(solver=solver,post_wall_s=time.monotonic()-poststart,archive=dict(path=str(target),bytes=target.stat().st_size,sha256=sha(target)),files={str(f.relative_to(run)):dict(sha256=sha(f),bytes=f.stat().st_size) for f in run.rglob('*') if f.is_file() and f.suffix!='.vtu'}))
 return 0
if __name__=='__main__':
 try:raise SystemExit(main(sys.argv[1],sys.argv[2]))
 except Exception:
  save(Path(sys.argv[1])/'runs'/sys.argv[2]/'failure.json',dict(traceback=traceback.format_exc()));raise
