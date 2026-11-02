"""Freeze the historical 14-run selection and start read-only analysis on frodo."""
import csv
from datetime import datetime,timezone
import hashlib
import json
from pathlib import Path
import subprocess

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'analysis/dense_history_20261102';OUT.mkdir(parents=True,exist_ok=True)
if (OUT/'plan.json').exists():raise SystemExit('Existing analysis plan: preserve evidence; do not relaunch in place')
phase=ROOT/'analysis/phase0_20261102'
snap=json.loads((phase/'remote_frodo_inventory.json').read_text())
with (phase/'run_audit.csv').open(encoding='utf-8-sig',newline='') as f:audit={r['run_key']:r for r in csv.DictReader(f)}
runs=[]
for r in snap['runs']:
    if r.get('vtu_count')!=101:continue
    a=audit[r['run_key']]
    entry={k:r[k] for k in ['run_key','path','parameters','small_file_hashes','sigmatop_sha256']}
    entry.update(id=f"r{len(runs)+1:02d}_c{a['case_id']}_n{a['seeded_nx']}",case_id=int(a['case_id']),frozen_audit=a)
    runs.append(entry)
assert len(runs)==14 and len({r['case_id'] for r in runs})==10
worker=(ROOT/'scripts/dense_history_worker.py').read_text(encoding='utf-8')
plan=dict(created_utc=datetime.now(timezone.utc).isoformat(),selection='all 14 frozen run records with 101 VTUs in phase0 snapshot',
    source_snapshot_sha256=hashlib.sha256((phase/'remote_frodo_inventory.json').read_bytes()).hexdigest(),
    source_audit_sha256=hashlib.sha256((phase/'run_audit.csv').read_bytes()).hexdigest(),
    worker_sha256=hashlib.sha256(worker.encode()).hexdigest(),runs=runs,simulations_started=False,
    method={'storage':'phi Sw x / horizontal area; no density multiplier','cuts':[.9,.75,.5],
            'quadrature':'tensor BOX nodal volumes; nearest adjacent-row midpoint cut',
            'inferred_downward_cumulative':'I(0)-I(t), assuming zero other salt boundary exchange and neglecting unmeasured residual',
            'threshold':'first observed frame >=0.1; prior sampled frame bounds observed transition, hidden earlier excursions not excluded',
            'no_effective_transport_label':True,'no_VTK_velocity_flux':True})
(OUT/'plan.json').write_text(json.dumps(plan,indent=2),encoding='utf-8')
remote='/temp/zhixin/paper2_history_analysis_20261102_r01'
code='''
from pathlib import Path
import subprocess,sys,os,json
root=Path(REMOTE);root.mkdir(exist_ok=False)
(root/'plan.json').write_text(json.dumps(PLAN,indent=2));(root/'worker.py').write_text(WORKER)
f=(root/'worker.log').open('w')
p=subprocess.Popen(['nice','-n','10',sys.executable,str(root/'worker.py'),str(root)],stdout=f,stderr=subprocess.STDOUT,
 stdin=subprocess.DEVNULL,start_new_session=True,env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1'))
print(json.dumps(dict(pid=p.pid,remote_root=str(root),runs=14,frames=1414,simulations_started=False)))
'''.replace('REMOTE',repr(remote)).replace('PLAN',repr(plan)).replace('WORKER',repr(worker))
bridge='/mnt/d/3OutputFile/1Research/Article/22Surrogate_model/paper/2/scripts/frodo_bridge.py'
p=subprocess.run(['wsl','python3',bridge,'run','python3 -'],input=code.encode(),capture_output=True,timeout=60)
if p.returncode:raise RuntimeError(p.stderr.decode(errors='replace'))
(OUT/'launch.json').write_bytes(p.stdout);print(p.stdout.decode())
