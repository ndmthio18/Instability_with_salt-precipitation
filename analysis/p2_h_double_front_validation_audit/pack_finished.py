import csv,json,tarfile,hashlib
from pathlib import Path
configs=list(csv.DictReader(open('final_cases.csv')))+list(csv.DictReader(open('fast_p01_cases.csv')))
Path('profiles_selected').mkdir(exist_ok=True);report=[];files=[]
for c in configs:
 name=c['id']
 if not (Path('raw')/(name+'_status.txt')).exists():continue
 history=Path('raw')/(name+'_history.csv');rows=list(csv.DictReader(history.open()));last=rows[-1]
 status_path=Path('raw')/(name+'_status.txt');status=dict(x.split('=',1) for x in status_path.read_text().splitlines())
 report.append(dict(id=name,status=status['status'],last_s=float(last['t']),target_s=float(c['end'])))
 for ending in ['history.csv','diagnostics.csv','status.txt','failure_state.csv']:
  p=Path('raw')/(name+'_'+ending)
  if p.exists():files.append(p)
 if name in {r['id'] for r in csv.DictReader(open('final_cases.csv'))}:
  p=Path('profiles_selected')/(name+'_profiles.csv')
  with (Path('raw')/(name+'_profiles.csv')).open() as f,p.open('w') as g:
   reader=csv.DictReader(f);writer=csv.DictWriter(g,fieldnames=reader.fieldnames);writer.writeheader()
   for r in reader:
    if float(r['t']) in [0,3600,43200,86400,float(last['t'])]:writer.writerow(r)
  files.append(p)
Path('partial_completion.json').write_text(json.dumps(report,indent=2));files.append(Path('partial_completion.json'))
files+=list(Path('logs').glob('*'))+list(Path('source_sqrt_v7').glob('*'))+list(Path('provenance').glob('*'))+list(Path('.').glob('*sha256.txt'))+list(Path('.').glob('*gate.json'))+list(Path('raw').glob('stress_repro*'))
with tarfile.open('partial_light.tar.gz','w:gz') as tar:
 for p in dict.fromkeys(files):
  if p.is_file():tar.add(p,arcname=str(p))
Path('partial_light.sha256').write_text(hashlib.sha256(Path('partial_light.tar.gz').read_bytes()).hexdigest()+'  partial_light.tar.gz\n')
print(json.dumps(report,indent=2))
