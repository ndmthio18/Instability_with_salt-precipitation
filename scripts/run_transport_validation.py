"""Launch eight bounded background/diagnostic checks on frodo, never Pilot A/B.
New directories only; sequential serial runs with pinned threads and timeouts.
"""
import ast,base64,hashlib,json,subprocess
from datetime import datetime,timezone
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];PARENT=ROOT.parents[1]
OUT=ROOT/'analysis/transport_validation_20261102';OUT.mkdir(exist_ok=True,parents=True)
if (OUT/'plan.json').exists():raise SystemExit('Validation plan already exists; preserve it and use a new version for retries')
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
constants_path=PARENT/'code/compare_onset_vs_paper_2025.py'
names={'phi','K','D','vg_n','Swr','vg_L','p_B','alpha','E','H','ell','g','x0','mu_ref'}
c={}
for node in ast.parse(constants_path.read_text(encoding='utf-8')).body:
    if isinstance(node,ast.Assign) and len(node.targets)==1 and isinstance(node.targets[0],ast.Name) and node.targets[0].id in names:
        c[node.targets[0].id]=ast.literal_eval(node.value)
assert set(c)==names
source=json.loads((ROOT/'analysis/phase0_20261102/pilot_candidates.json').read_text(encoding='utf-8'))[0]['actual_input_parameters']
fluid={k:source['Problem.'+k] for k in ['PureWaterDensity','DensityLinearCoefficient','DensityQuadraticCoefficient','MolarMassH2O','MolarMassNaCl']}
x0=c['x0'];mw=float(fluid['MolarMassH2O']);ms=float(fluid['MolarMassNaCl']);w=x0*ms/((1-x0)*mw+x0*ms)
rho=float(fluid['PureWaterDensity'])+1000*w*(float(fluid['DensityLinearCoefficient'])+float(fluid['DensityQuadraticCoefficient'])*w)
base={'Grid':{'Positions0':f"0 {c['ell']}",'Positions1':f"0 {c['H']}",'Cells0':8,'Cells1':40,'Grading1':-0.98},
 'TimeLoop':{'TEnd':2000,'DtInitial':20,'MaxTimeStepSize':20,'OutputInterval':100000},
 'Problem':{'Name':'background','Permeability':c['K'],'Porosity':c['phi'],'Swr':c['Swr'],'PcN':c['vg_n'],
 'PcAlpha':c['alpha'],'SweL':c['vg_L'],'KrwMualemExponent':'true','WaterLevel':c['p_B']/(rho*c['g']),
 'EvaporationRate':c['E'],'SaltDiffusivity':c['D'],'ReferenceViscosity':c['mu_ref'],'SaltMoleFraction':c['x0'],
 'GravitationalAcceleration':c['g'],**fluid,'Amplitude':0,'WaveNumber':1,'StandardDeviationSaltMoleFraction':0,
 'PerturbationType':'Constant','RandomSeed':1,'EnableGravity':'true','PerturbationHeight':0},
 'Transport':{'SolubilityMoleFraction':0.1,'StopAtSolubility':'true'},'Vtk':{'Precision':'Float64','AddVelocity':'true'}}
specs=[
 ('closed_hydrostatic',{'Problem.EvaporationRate':0,'TimeLoop.TEnd':100,'TimeLoop.DtInitial':10,'TimeLoop.MaxTimeStepSize':10},'requested_horizon'),
 ('evap_nx8_dt20',{},'requested_horizon'),
 ('evap_nx16_dt20',{'Grid.Cells0':16},'requested_horizon'),
 ('evap_nx8_dt10',{'TimeLoop.DtInitial':10,'TimeLoop.MaxTimeStepSize':10},'requested_horizon'),
 ('threshold_dt2',{'Problem.SaltMoleFraction':.09999,'TimeLoop.DtInitial':2,'TimeLoop.MaxTimeStepSize':2},'surface_solubility_threshold'),
 ('threshold_dt1',{'Problem.SaltMoleFraction':.09999,'TimeLoop.DtInitial':1,'TimeLoop.MaxTimeStepSize':1},'surface_solubility_threshold'),
 ('threshold_initial',{'Problem.SaltMoleFraction':.1},'surface_solubility_threshold'),
 ('startup_invalid_time',{'TimeLoop.TEnd':-1},'startup_failure')]
runs=[]
for name,updates,reason in specs:
    v=json.loads(json.dumps(base));v['Problem']['Name']=name
    for key,value in updates.items():sec,k=key.split('.',1);v[sec][k]=value
    text='\n\n'.join('['+sec+']\n'+'\n'.join(f'{k} = {value}' for k,value in values.items()) for sec,values in v.items())+'\n'
    runs.append({'id':name,'input':text,'parameters':v,'expected_reason':reason,'expected_exit':1 if reason=='startup_failure' else 0,
                 'input_sha256':hashlib.sha256(text.encode()).hexdigest(),'wall_timeout_s':600})
build=json.loads((ROOT/'analysis/transport_build_20261102/attempt_02/result.json').read_text())
assert all(sha(ROOT/'model/transport_v1/source'/n)==h for n,h in build['source_sha256'].items())
plan={'created_utc':datetime.now(timezone.utc).isoformat(),'scope':'small uniform 2D background and diagnostic validation; no Pilot A/B',
 'constant_source':str(constants_path),'constant_source_sha256':sha(constants_path),'constants':c,'reference_density_kg_m3':rho,
 'fluid_source':'phase0 pilot_candidates.json actual input fluid constants only; no A/B case simulations',
 'fluid_source_sha256':sha(ROOT/'analysis/phase0_20261102/pilot_candidates.json'),
 'executable':build['executable'],'build_source_sha256':build['source_sha256'],'runs':runs,
 'criteria':{'relative_cumulative_inventory_residual_max':1e-6,'relative_residual_to_cumulative_absolute_flux_max':1e-3,
 'absolute_residual_floor_m':1e-12,'lateral_x_span_max':1e-9,'lateral_sw_span_max':1e-9,'lateral_pressure_span_Pa_max':1e-4,
 'horizontal_refinement_profile_x_relative_max':1e-5,'horizontal_refinement_inventory_relative_max':1e-6,
 'dt_halving_top_x_relative_max':5e-3,'vtk_inventory_relative_max':1e-10,'hydrostatic_flux_abs_m_s_max':1e-12},
 'event_tests_are_manufactured_high_initial_salinity':True,'max_concurrent_processes':1,'nice':10}
(OUT/'plan.json').write_text(json.dumps(plan,indent=2),encoding='utf-8')
inp=OUT/'inputs';inp.mkdir()
for r in runs:(inp/(r['id']+'.input')).write_text(r['input'],encoding='utf-8')
worker=r'''
import csv,json,os,subprocess,time,hashlib,datetime,xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
root=Path('/temp/zhixin/paper2_transport_20261102/validation_20261102_r01')
plan=json.loads((root/'plan.json').read_text());exe=Path(plan['executable']['path'])
assert hashlib.sha256(exe.read_bytes()).hexdigest()==plan['executable']['sha256']
def fieldcheck(folder,phi):
 meta=json.loads((folder/'transport_metadata.json').read_text());records=[];profiles=[]
 for p in sorted(folder.glob('*.vtu')):
  tree=ET.parse(p);arrays={}
  for a in tree.iter('DataArray'):
   if a.attrib.get('format')=='ascii':arrays[a.attrib.get('Name','unnamed')]=np.fromstring(a.text or '',sep=' ')
  xyz=arrays['Coordinates'].reshape((-1,3));xx=np.unique(xyz[:,0]);zz=np.unique(xyz[:,1]);width=xx[-1]-xx[0]
  def widths(a):return np.diff(np.r_[a[0],.5*(a[:-1]+a[1:]),a[-1]])
  weights=widths(xx)[np.searchsorted(xx,xyz[:,0])]*widths(zz)[np.searchsorted(zz,xyz[:,1])]/width
  water=phi*arrays['sw']*weights;salt=water*arrays['moleFraction']
  row={'file':p.name,'sha256':hashlib.sha256(p.read_bytes()).hexdigest(),'all_salt_m':float(salt.sum()),'all_water_m':float(water.sum()),'regions':{},'lateral_max_span':{}}
  for cut in meta['cuts']:
   mask=xyz[:,1]>zz[0]+cut['actual_z_over_H']*(zz[-1]-zz[0]);row['regions'][cut['name']]={'salt_m':float(salt[mask].sum()),'water_m':float(water[mask].sum())}
  for name in ['moleFraction','sw','pw']:
   spans=[]
   for z in zz:
    values=arrays[name][xyz[:,1]==z];spans.append(float(np.ptp(values)))
    profiles.append({'frame':p.name,'z_m':float(z),'field':name,'mean':float(values.mean()),'span':float(np.ptp(values))})
   row['lateral_max_span'][name]=max(spans)
  records.append(row)
 (folder/'field_checks.json').write_text(json.dumps(records,indent=2))
 if profiles:
  with (folder/'field_profiles.csv').open('w',newline='') as f:
   w=csv.DictWriter(f,fieldnames=list(profiles[0]));w.writeheader();w.writerows(profiles)
summary=[]
for spec in plan['runs']:
 folder=root/'runs'/spec['id'];folder.mkdir(parents=True,exist_ok=False)
 data=spec['input'].encode();assert hashlib.sha256(data).hexdigest()==spec['input_sha256'];(folder/'params.input').write_bytes(data)
 manifest={k:v for k,v in spec.items() if k!='input'};manifest.update(executable_sha256=plan['executable']['sha256'],build_source_sha256=plan['build_source_sha256'],purpose=plan['scope'],status='running',started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
 (folder/'run_manifest.json').write_text(json.dumps(manifest,indent=2));started=time.monotonic()
 result={'id':spec['id'],'expected_reason':spec['expected_reason'],'expected_exit':spec['expected_exit']}
 with (folder/'solver.log').open('w') as log:
  try:
   p=subprocess.run(['nice','-n','10',str(exe),'params.input'],cwd=folder,stdout=log,stderr=subprocess.STDOUT,timeout=spec['wall_timeout_s'],env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1'))
   result['returncode']=p.returncode
  except subprocess.TimeoutExpired:result.update(returncode=None,failure='wall_timeout')
 result['wall_s']=time.monotonic()-started
 if (folder/'stop_record.json').exists():result['stop']=json.loads((folder/'stop_record.json').read_text())
 if (folder/'transport_metadata.json').exists():
  try:fieldcheck(folder,float(spec['parameters']['Problem']['Porosity']))
  except Exception as e:result['field_check_error']=repr(e)
 result['files']={p.name:{'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in folder.iterdir() if p.is_file() and p.name!='run_manifest.json'}
 manifest.update(status='finished',result=result,finished_utc=datetime.datetime.now(datetime.timezone.utc).isoformat());(folder/'run_manifest.json').write_text(json.dumps(manifest,indent=2))
 summary.append(result);(root/'summary.json').write_text(json.dumps(summary,indent=2))
 if result.get('returncode')!=spec['expected_exit'] or result.get('stop',{}).get('reason')!=spec['expected_reason']:
  (root/'halted.json').write_text(json.dumps({'reason':'unexpected_exit_or_stop','run':result},indent=2));break
(root/'finished.json').write_text(json.dumps({'completed':len(summary),'planned':len(plan['runs']),'full_pilot_started':False}))
'''
(OUT/'remote_worker.py').write_text(worker,encoding='utf-8')
remote=r'''
from pathlib import Path
import json,subprocess,sys,hashlib
root=Path('/temp/zhixin/paper2_transport_20261102/validation_20261102_r01');root.mkdir(exist_ok=False)
plan=PLAN
exe=Path(plan['executable']['path']);assert hashlib.sha256(exe.read_bytes()).hexdigest()==plan['executable']['sha256']
(root/'plan.json').write_text(json.dumps(plan,indent=2));(root/'worker.py').write_text(WORKER)
f=(root/'worker.log').open('w');p=subprocess.Popen([sys.executable,str(root/'worker.py')],stdout=f,stderr=subprocess.STDOUT,stdin=subprocess.DEVNULL,start_new_session=True)
print(json.dumps({'pid':p.pid,'root':str(root),'planned_runs':len(plan['runs']),'concurrent':1}))
'''.replace('PLAN',repr(plan)).replace('WORKER',repr(worker))
bridge='/mnt/d/3OutputFile/1Research/Article/22Surrogate_model/paper/2/scripts/frodo_bridge.py'
p=subprocess.run(['wsl','python3',bridge,'run','python3 -'],input=remote.encode(),capture_output=True,timeout=60)
if p.returncode:raise RuntimeError(p.stderr.decode(errors='replace'))
(OUT/'launch.json').write_bytes(p.stdout);print(p.stdout.decode())
