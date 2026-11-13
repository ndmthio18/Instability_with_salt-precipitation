from pathlib import Path
import pandas as pd
from analyze import ROOT,R
mapping={'diagnostic_cases.csv':'double_v1','p0_cases.csv':'longdouble_v2','early_cases.csv':'early_v3','retry_cases.csv':'setup_v6','p2_cases.csv':'setup_v6','production_cases.csv':'setup_v6','fast_p01_cases.csv':'sqrt_v7','final_cases.csv':'sqrt_v7','parallel_cases.csv':'parallel_v8'}
rows=[]
for filename,version in mapping.items():
 for r in pd.read_csv(ROOT/filename).to_dict('records'):
  r['implementation']=version;r['configuration_file']=filename
  status=ROOT/'raw'/f'{r["id"]}_status.txt'
  if not status.exists():status=ROOT/'provenance/run_status'/f'{r["id"]}_status.txt'
  if status.exists():r.update(dict(x.split('=',1) for x in status.read_text().splitlines()))
  else:r['status']='RAW_EVIDENCE_REMOTE_ONLY'
  r['local_history']=(ROOT/'raw'/f'{r["id"]}_history.csv').exists()
  r['remote_raw_file']='raw/'+r['id']+'_history.csv';rows.append(r)
for name,rt in [('legacy_r8',1e-8),('legacy_r9',1e-9),('legacy_r10',1e-10)]:
 r=dict(id=name,implementation='source_legacy',configuration_file='reproduce.sbatch',nc=128,nl=32,kappa=1,gap=1e-9,rtol=rt,u=.6,end=100)
 r.update(dict(x.split('=',1) for x in (ROOT/'raw'/f'{name}_status.txt').read_text().splitlines()));rows.append(r)
pd.DataFrame(rows).to_csv(R/'case_registry.csv',index=False)
