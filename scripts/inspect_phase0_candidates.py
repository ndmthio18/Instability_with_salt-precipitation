"""Read a few existing ASCII VTU files remotely; never run the solver."""
import json,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'analysis/phase0_20261102'
snap=json.loads((OUT/'remote_frodo_inventory.json').read_text())
keys=['frodo:c1b/runs_600/shard-001/0017--case-36877--native',
      'frodo:c1/runs_36877_nx60/shard-000/0001--case-36877--native',
      'frodo:hf1/runs/shard-019/0079--case-38852--native',
      'frodo:paper1_wp1_paired_v2/runs/case-38852--both--nx-16']
targets=[r for r in snap['runs'] if r['run_key'] in keys]
code=r'''
import json,hashlib,xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
out=[]
for r in TARGETS:
 p=Path(r['path']);fs=sorted(p.glob('*.vtu'));q=fs[-1]
 tree=ET.parse(q);data={};arrays={}
 for x in tree.iter('DataArray'):
  if x.attrib.get('format')!='ascii':continue
  name=x.attrib.get('Name','unnamed');a=np.fromstring(x.text or '',sep=' ')
  arrays[name]=a
  data[name]={'count':a.size,'min':float(a.min()) if a.size else None,'max':float(a.max()) if a.size else None}
 x=arrays['moleFraction'];rho=arrays['rhoW']
 def density(r0,mw,ms,lin,quad):
  w=x*ms/((1-x)*mw+x*ms);return r0+1000*w*(lin+quad*w)
 rhoerr={'max_abs_error_explicit_input_law_kg_m3':float(np.max(np.abs(rho-density(998.203,.01801528,.05844,.670804,.373854)))),
         'max_abs_error_legacy_defaults_kg_m3':float(np.max(np.abs(rho-density(1000,.018,.0584428,.6708,.3739))))}
 xyz=arrays['Coordinates'].reshape((-1,3));top=np.isclose(xyz[:,1],.2,rtol=0,atol=1e-12)
 out.append({'run_key':r['run_key'],'file':str(q),'sha256':hashlib.sha256(q.read_bytes()).hexdigest(),'arrays':data,'density_check':rhoerr,
             'final_top_mean':float(x[top].mean()),'final_top_max':float(x[top].max()),
             'input_text':next(p.glob('*.input')).read_text(),'eigenprofile_sha256':hashlib.sha256(next(p.glob('*_eigenprofile.dat')).read_bytes()).hexdigest()})
print(json.dumps(out))
'''.replace('TARGETS',repr([{'run_key':r['run_key'],'path':r['path']} for r in targets]))
bridge='/mnt/d/3OutputFile/1Research/Article/22Surrogate_model/paper/1/analysis/onset_smoothness_audit_20261021/remote_bridge.py'
r=subprocess.run(['wsl','python3',bridge,'run','python3 -'],input=code.encode(),capture_output=True,timeout=100)
if r.returncode:raise RuntimeError(r.stderr.decode(errors='replace'))
j=json.loads(r.stdout)
(OUT/'candidate_field_evidence.json').write_text(json.dumps(j,ensure_ascii=False,indent=2),encoding='utf-8')
for x in j:print(x['run_key'],x['arrays'])
