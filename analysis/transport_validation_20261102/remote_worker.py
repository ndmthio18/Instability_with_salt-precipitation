
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
