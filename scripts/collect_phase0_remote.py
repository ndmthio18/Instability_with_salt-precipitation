"""Read existing frodo runs over SSH; no simulations or remote writes.

Bounded inventory: only frozen handoff paths and the 84 named WP1 run paths.
VTU header inspection reads at most 64 KiB per inspected file, not field arrays.
Uses the existing WSL SSH bridge, with strict host-key verification.
"""
import base64
import csv
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PAPER1 = ROOT.parents[1]
OUT = ROOT / 'analysis/phase0_20261102'
OUT.mkdir(parents=True, exist_ok=True)
handoff = ROOT / 'data/external/paper1/158395f3d25a/handoff'
runs = list(csv.DictReader((handoff / 'dumux_onsets.csv').open(encoding='utf-8-sig')))
targets = [{'run_key': r['run_key'], 'path': '/temp/zhixin/' + r['run_key'].split(':', 1)[1]}
           for r in runs if r['host'] == 'frodo']
for p in sorted((PAPER1 / 'paper/1/analysis/wp1_paired_v2_raw/runs').glob('*/case_manifest.json')):
    targets.append({'run_key': 'frodo:paper1_wp1_paired_v2/runs/' + p.parent.name,
                    'path': '/temp/zhixin/paper1_wp1_paired_v2/runs/' + p.parent.name})

REMOTE = r'''
import json,hashlib,re,csv,io,datetime,xml.etree.ElementTree as ET
from pathlib import Path
targets=TARGETS
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def inspect(t):
 p=Path(t['path']);r=dict(t,exists=p.is_dir())
 if not p.is_dir():return r
 fs=sorted(q for q in p.iterdir() if q.is_file());r['files']=[q.name for q in fs]
 r['small_file_hashes']={q.name:sha(q) for q in fs if q.suffix in ('.input','.json','.dat','.pvd') and q.stat().st_size<2000000}
 ins=[q for q in fs if q.suffix=='.input'];r['parameters']={}
 if ins:
  section=''
  for line in ins[0].read_text(errors='replace').splitlines():
   line=line.split('#')[0].strip()
   if line.startswith('['):section=line.strip('[]')
   elif '=' in line:
    k,v=line.split('=',1);r['parameters'][section+'.'+k.strip()]=v.strip()
 r['pvd_times_s']=[]
 for q in fs:
  if q.suffix=='.pvd':
   try:r['pvd_times_s'] += [float(x.attrib['timestep']) for x in ET.parse(q).iter('DataSet')]
   except Exception as e:r['pvd_error']=str(e)
 vtus=[q for q in fs if q.suffix=='.vtu'];r['vtu_count']=len(vtus)
 r['vtu_header_arrays']={}
 for q in ([vtus[0],vtus[-1]] if vtus else []):
  with q.open('rb') as h:b=h.read(65536).split(b'<AppendedData')[0]
  r['vtu_header_arrays'][q.name]=re.findall(r'<DataArray[^>]*Name="([^"]+)"',b.decode('utf-8',errors='replace'))
 sig=[q for q in fs if q.name.endswith('_sigmatop.csv')]
 if sig:
  q=sig[0];r['sigmatop_sha256']=sha(q)
  lines=q.read_text().splitlines();r['sigmatop_header']=[x for x in lines if x.startswith('#')][:3]
  vals=[]
  for line in lines:
   if not line.strip() or line.startswith('#'):continue
   try:vals.append([float(v) for v in line.replace(',',' ').split()])
   except ValueError:pass
  if vals:
   r['sigmatop_n']=len(vals);r['sigmatop_first']=vals[0];r['sigmatop_last']=vals[-1]
   r['sigmatop_column_min']=[min(row[i] for row in vals) for i in range(len(vals[0]))]
   r['sigmatop_column_max']=[max(row[i] for row in vals) for i in range(len(vals[0]))]
   # Confirm header; these diagnostics do NOT contain the local maximum.
   if any('time[s] x_top_mean[mol/mol] sigma_top[mol/mol]' in h for h in r['sigmatop_header']):
    r['first_top_mean_ge_0p1_s']=next((row[0] for row in vals if row[1]>=0.1),None)
   r['first_top_max_ge_0p1_s']=None
   r['top_max_status']='not_recorded_in_sigmatop'
 for name in ['case_manifest.json','run_meta.json']:
  q=p/name
  if q.exists():
   try:
    j=json.loads(q.read_text());r[name]=j if name=='run_meta.json' else {'dumux':j.get('dumux'),'dumux_numerics':j.get('dumux_numerics'),'constitutive_audit':j.get('constitutive_audit')}
   except Exception as e:r[name]={'error':str(e)}
 logs=[q for q in fs if q.suffix=='.log']
 if logs:
  with logs[0].open('rb') as h:
   h.seek(max(0,logs[0].stat().st_size-4000));r['log_tail']=h.read().decode(errors='replace')
 return r
base=Path('/temp/zhixin/surrogate/code/dumux_kiemle2023b_unsat/kiemle2023b-main/appl/simpleRichardsnc')
sources={str(base/n):{'sha256':sha(base/n),'text':(base/n).read_text()} for n in ['model.hh','main.cc']}
print(json.dumps({'host':'frodo','read_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'scope':'named paths only; read-only; no simulation','sources':sources,'runs':[inspect(t) for t in targets]}))
'''.replace('TARGETS', repr(targets))
command = 'python3 -'
bridge = '/mnt/d/3OutputFile/1Research/Article/22Surrogate_model/paper/1/analysis/onset_smoothness_audit_20261021/remote_bridge.py'
result = subprocess.run(['wsl','python3',bridge,'run',command],input=REMOTE.encode(),capture_output=True,timeout=120)
if result.returncode:
    raise RuntimeError(result.stderr.decode(errors='replace'))
data = json.loads(result.stdout)
(OUT / 'remote_frodo_inventory.json').write_text(json.dumps(data,ensure_ascii=False,indent=2),encoding='utf-8')
print(json.dumps({'targets':len(targets),'found':sum(r['exists'] for r in data['runs']),
                  'vtu_counts':{str(n):sum(r.get('vtu_count')==n for r in data['runs']) for n in sorted({r.get('vtu_count',0) for r in data['runs']})}},indent=2))
