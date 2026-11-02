"""One explicitly requested unperturbed run; frozen legacy grid and time window."""
import hashlib,json,subprocess
from pathlib import Path
from onset_semantics import FrozenOnsets, standard_fields
from datetime import datetime,timezone
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'analysis/unperturbed_36877_20261102'
OUT.mkdir(exist_ok=False)
c=json.loads((ROOT/'analysis/phase0_20261102/pilot_candidates.json').read_text())[0]
v=c['actual_input_parameters'].copy()
updates={'Problem.Name':'unperturbed_36877','Problem.PerturbationType':'Constant','Problem.Amplitude':'0','Problem.StandardDeviationSaltMoleFraction':'0','TimeLoop.OutputInterval':'50','Transport.StopAtSolubility':'true','Transport.SolubilityMoleFraction':'0.1'}
v.update(updates);v.pop('Problem.EigenprofilePath',None)
sections={}
for key,value in v.items():
 s,k=key.split('.',1);sections.setdefault(s,{})[k]=value
inp='\n\n'.join('['+s+']\n'+'\n'.join(k+' = '+v for k,v in vals.items()) for s,vals in sections.items())+'\n'
(OUT/'params.input').write_text(inp)
build=json.loads((ROOT/'analysis/transport_build_20261102/attempt_02/result.json').read_text())
for n,h in build['source_sha256'].items():assert hashlib.sha256((ROOT/'model/transport_v1/source'/n).read_bytes()).hexdigest()==h
plan={'created_utc':datetime.now(timezone.utc).isoformat(),'case_id':36877,'purpose':'observe whether numerical noise grows without imposed perturbation','parameters':v,'changes_from_legacy':updates,'removed_unused_key':'Problem.EigenprofilePath','source_run':c['preferred_existing_run_key'],'executable':build['executable'],'source_sha256':build['source_sha256'],'input_sha256':hashlib.sha256(inp.encode()).hexdigest(),'legacy_detect_s':50394.334,'horizon_s':float(v['TimeLoop.TEnd']),'claim_scope':'same legacy time window; absence of growth is not proof of stability','predeclared_detection':'surface sigma > 1e-8 for 5 consecutive accepted steps and above 1000 times initial noise floor; descriptive visibility criterion, not LSA onset','wall_timeout_s':43200}
plan['legacy_onset_reference']=standard_fields(FrozenOnsets().run_by_key(c['preferred_existing_run_key']))
(OUT/'plan.json').write_text(json.dumps(plan,indent=2))
remote='/temp/zhixin/paper2_transport_20261102/unperturbed_36877_20261102_r01'
worker='''import json,subprocess,os,time,hashlib,datetime
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
'''
(OUT/'worker.py').write_text(worker)
code="""import json,subprocess,sys
from pathlib import Path
r=Path(REMOTE);r.mkdir(exist_ok=False)
(r/'params.input').write_text(INPUT)
(r/'plan.json').write_text(json.dumps(PLAN,indent=2))
(r/'worker.py').write_text(WORKER)
with (r/'worker.log').open('w') as f:
 p=subprocess.Popen([sys.executable,str(r/'worker.py')],stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
print(json.dumps({'pid':p.pid,'root':str(r)}))
""".replace('REMOTE',repr(remote)).replace('INPUT',repr(inp)).replace('PLAN',repr(plan)).replace('WORKER',repr(worker))
bridge='/mnt/d/3OutputFile/1Research/Article/22Surrogate_model/paper/2/scripts/frodo_bridge.py'
p=subprocess.run(['wsl','python3',bridge,'run','python3 -'],input=code.encode(),capture_output=True,timeout=60)
if p.returncode:raise RuntimeError(p.stderr.decode())
(OUT/'launch.json').write_bytes(p.stdout);print(p.stdout.decode())
