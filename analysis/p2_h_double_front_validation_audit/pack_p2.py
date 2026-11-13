import csv,json,tarfile,hashlib,os
from pathlib import Path
configs=list(csv.DictReader(open('p2_cases.csv')));reports=[]
if os.getenv('CORE_P2_ONLY'):configs=[c for c in configs if c['id'].startswith('p2_')]
stem='p2_kinetic' if os.getenv('CORE_P2_ONLY') else 'p2'
Path('profiles_selected').mkdir(exist_ok=True)
for c in configs:
 name=c['id'];p=Path('raw')/(name+'_history.csv');rows=list(csv.DictReader(p.open()));last=rows[-1]
 status=dict(x.split('=',1) for x in (Path('raw')/(name+'_status.txt')).read_text().splitlines())
 residuals=['salt_res','water_res','solid_res','liquid_salt_res','crust_salt_res','liquid_water_res','crust_water_res','water_jump','contact_res']
 worst=max(abs(float(r[k])) for r in rows for k in residuals)
 jump=max(max(abs(float(r['Q_SPLUS'])-float(r['Q_SMINUS'])-float(r['Q_DB'])),abs(float(r['Q_MPLUS'])-float(r['Q_MMINUS'])-float(r['Q_DB']))) for r in rows)
 reports.append(dict(id=name,status=status['status'],end_s=float(last['t']),max_budget=worst,max_integrated_jump=jump,final_vt=float(last['vt']),final_vb=float(last['vb']),final_vnet=float(last['vnet']),final_hc=float(last['hc'])))
 with (Path('raw')/(name+'_profiles.csv')).open() as f,(Path('profiles_selected')/(name+'_profiles.csv')).open('w') as g:
  reader=csv.DictReader(f);writer=csv.DictWriter(g,fieldnames=reader.fieldnames);writer.writeheader()
  for r in reader:
   if float(r['t']) in [0,3600,43200,86400,float(last['t'])]:writer.writerow(r)
result=dict(cases=reports,passed=all(r['status']=='end' and r['max_budget']<2e-9 and r['max_integrated_jump']<2e-9 for r in reports))
Path((stem+'_gate.json')).write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
files=[Path((stem+'_gate.json'))]+list(Path('profiles_selected').glob('*'))
for c in configs:files+=list(Path('raw').glob(c['id']+'*history.csv'))+list(Path('raw').glob(c['id']+'*diagnostics.csv'))+list(Path('raw').glob(c['id']+'*status.txt'))+list(Path('raw').glob(c['id']+'*failure_state.csv'))
files+=list(Path('logs').glob('p2.*'))
with tarfile.open((stem+'_light.tar.gz'),'w:gz') as tar:
 for p in files:tar.add(p,arcname=str(p))
Path((stem+'_light.sha256')).write_text(hashlib.sha256(Path((stem+'_light.tar.gz')).read_bytes()).hexdigest()+'  '+stem+'_light.tar.gz\n')
