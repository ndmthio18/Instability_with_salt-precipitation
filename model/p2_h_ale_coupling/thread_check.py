"""Operational thread check: exact same 120 s case, immutable new directory."""
import pathlib,subprocess,json,os,time,hashlib
root=pathlib.Path.cwd();path=root/'thread_check';path.mkdir(exist_ok=False)
ref=root/'smoke/smoke_n320';(path/'params.input').write_bytes((ref/'params.input').read_bytes())
args=json.loads((ref/'command.json').read_text());(path/'command.json').write_text(json.dumps(args,indent=2))
t=time.monotonic()
with (path/'solver.out').open('w') as f:
    p=subprocess.Popen(args,cwd=path,stdout=f,stderr=subprocess.STDOUT)
    time.sleep(.1)
    status=pathlib.Path(f'/proc/{p.pid}/status')
    (path/'process_status.txt').write_text(status.read_text() if status.exists() else 'already finished')
    code=p.wait(timeout=240)
elapsed=time.monotonic()-t
name='coupled_history.csv';same=(path/name).read_bytes()==(ref/name).read_bytes()
r=dict(exit_code=code,elapsed_s=elapsed,identical_history=same,sha256=hashlib.sha256((path/name).read_bytes()).hexdigest(),DUMUX_NUM_THREADS=os.environ.get('DUMUX_NUM_THREADS'))
(path/'result.json').write_text(json.dumps(r,indent=2));print(json.dumps(r));assert code==0 and same
