import json,subprocess,os,time,hashlib,datetime
from pathlib import Path
r=Path(__file__).resolve().parent
p=json.loads((r/'plan.json').read_text());exe=Path(p['executable']['path'])
assert hashlib.sha256(exe.read_bytes()).hexdigest()==p['executable']['sha256']
assert hashlib.sha256((r/'params.input').read_bytes()).hexdigest()==p['input_sha256']
start=time.monotonic()
with (r/'solver.log').open('w') as f:
 try:
  proc=subprocess.Popen(['nice','-n','10',str(exe),'params.input'],cwd=r,stdout=f,stderr=subprocess.STDOUT,env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1'))
  (r/'solver_pid.json').write_text(json.dumps({'pid':proc.pid}))
  rc=proc.wait(timeout=p['wall_timeout_s'])
 except subprocess.TimeoutExpired:
  proc.terminate();proc.wait();rc='wall_timeout'
result={'returncode':rc,'wall_s':time.monotonic()-start,'finished_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'stop':json.loads((r/'stop_record.json').read_text()) if (r/'stop_record.json').exists() else None}
result['files']={q.name:{'bytes':q.stat().st_size,'sha256':hashlib.sha256(q.read_bytes()).hexdigest()} for q in r.iterdir() if q.is_file()}
(r/'finished.json').write_text(json.dumps(result,indent=2))
