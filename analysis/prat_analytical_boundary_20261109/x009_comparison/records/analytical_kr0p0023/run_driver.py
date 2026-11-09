"""User-requested x=0.09 comparison; same verified analytical binary, no recompilation."""
import csv,hashlib,json,os,pathlib,subprocess,time
root=pathlib.Path.cwd()
if not os.environ.get('SLURM_JOB_ID') or os.uname().nodename.startswith(('ibcln','login','test')): raise RuntimeError('compute allocation required')
exe=root/'build/soil_solver'
assert hashlib.sha256(exe.read_bytes()).hexdigest()=='69b2a5b7f5d580f66bc43a0ba7b619ae35d4906aaa43885faa72cdc65debbccf'
d=root/'runs/x009_analytical_n80_dt0p5'; d.mkdir(parents=True,exist_ok=False)
(d/'params.input').write_text((root/'repo/model/prat_analytical_boundary/column.input').read_text())
args=[str(exe),'params.input','-Problem.SaltMoleFraction','0.09']
(d/'command.json').write_text(json.dumps(args,indent=2))
(d/'case_spec.json').write_text(json.dumps({'initial_mole_fraction':.09,'reservoir_mole_fraction':.09,'soil_cells':80,
    'dt_max_s':.5,'requested_s':25000,'reaction_coefficient_m_s':.0023,'activation_m':3.5e-6,'porosity':.43,
    'internal_precipitation':False,'crust_grid':False,'binary_sha256':hashlib.sha256(exe.read_bytes()).hexdigest()},indent=2))
started=time.monotonic()
with (d/'solver.out').open('w') as f:
    try: code=subprocess.run(args,cwd=str(d),stdout=f,stderr=subprocess.STDOUT,timeout=16*3600).returncode
    except subprocess.TimeoutExpired: code=124
rows=list(csv.DictReader((d/'history.csv').open())); last=rows[-1]
result={'name':d.name,'exit_code':code,'job':os.environ['SLURM_JOB_ID'],'node':os.uname().nodename,
    'elapsed_s':time.monotonic()-started,'requested_s':25000,'last':last,
    'max_salt_residual':max(abs(float(x['salt_residual_rel'])) for x in rows),
    'max_water_residual':max(abs(float(x['water_residual_rel'])) for x in rows)}
result['completed']=code==0 and abs(float(last['time_s'])-25000)<1e-6
result['coupling_pass']=result['completed'] and last['active']=='1' and max(result['max_salt_residual'],result['max_water_residual'])<=1e-8
(d/'outcome.json').write_text(json.dumps(result,indent=2)); print(json.dumps(result),flush=True)
