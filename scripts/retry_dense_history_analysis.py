"""Preserve failed analysis r01, repair JSON serialization, restart only postprocessing."""
import hashlib,json,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'analysis/dense_history_20261102'
if (OUT/'launch_r02.json').exists():raise SystemExit('Retry already recorded')
old=json.loads((OUT/'launch.json').read_text());plan=json.loads((OUT/'plan.json').read_text())
worker=(ROOT/'scripts/dense_history_worker.py').read_text()
revision=dict(reason='numpy.bool_ in audit checks was not JSON serializable; convert to builtin bool',
              original_worker_sha256=plan['worker_sha256'],actual_worker_sha256=hashlib.sha256(worker.encode()).hexdigest(),
              scientific_method_changed=False,simulations_started=False)
code='''
import os,signal,json,subprocess,sys,tarfile
from pathlib import Path
old=OLD;root_old=Path(old['remote_root']);pid=old['pid']
cmd=Path('/proc')/str(pid)/'cmdline'
if cmd.exists():
 args=cmd.read_bytes().split(b'\\0')
 assert str(root_old/'worker.py').encode() in args and str(root_old).encode() in args
 os.kill(pid,signal.SIGTERM)
(root_old/'interrupted_for_analysis_repair.json').write_text(json.dumps(REVISION,indent=2))
with tarfile.open(str(root_old)+'_failed.tar.gz','w:gz') as t:
 for p in root_old.iterdir():t.add(p,arcname=p.name)
root=Path('/temp/zhixin/paper2_history_analysis_20261102_r02');root.mkdir(exist_ok=False)
(root/'plan.json').write_text(json.dumps(PLAN,indent=2));(root/'worker.py').write_text(WORKER)
(root/'worker_revision.json').write_text(json.dumps(REVISION,indent=2))
f=(root/'worker.log').open('w')
p=subprocess.Popen(['nice','-n','10',sys.executable,str(root/'worker.py'),str(root)],stdout=f,stderr=subprocess.STDOUT,
 stdin=subprocess.DEVNULL,start_new_session=True,env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1'))
print(json.dumps(dict(pid=p.pid,remote_root=str(root),failed_archive=str(root_old)+'_failed.tar.gz',simulations_started=False)))
'''.replace('OLD',repr(old)).replace('PLAN',repr(plan)).replace('WORKER',repr(worker)).replace('REVISION',repr(revision))
bridge='/mnt/d/3OutputFile/1Research/Article/22Surrogate_model/paper/2/scripts/frodo_bridge.py'
p=subprocess.run(['wsl','python3',bridge,'run','python3 -'],input=code.encode(),capture_output=True,timeout=60)
if p.returncode:raise RuntimeError(p.stderr.decode(errors='replace'))
(OUT/'launch_r02.json').write_bytes(p.stdout);(OUT/'worker_revision.json').write_text(json.dumps(revision,indent=2))
print(p.stdout.decode())
