from pathlib import Path
import json,hashlib,subprocess
import pandas as pd
ROOT=Path(__file__).resolve().parent
checks=[]
def check(ok,message):
 if not ok:raise AssertionError(message)
 checks.append(message)
for name in ['README.md','NUMERICAL_FAILURE_DIAGNOSIS.md','INTERFACE_FLUX_CONVERGENCE.md','BOTTOM_KINETICS_AUDIT.md','TOP_AND_NET_GROWTH_VALIDATION.md','CORE_SOLID_FRACTION_AUDIT.md','FINAL_REVIEW.md']:
 check((ROOT/name).exists(),name)
status=json.loads((ROOT/'results/validation_status.json').read_text());check(len(status['gates'])==8,'eight separate gates');check(status['primary_completed']==51,'51 complete primary cases')
ids=list(pd.read_csv(ROOT/'fast_p01_cases.csv').id)+list(pd.read_csv(ROOT/'primary_cases.csv').id)
manifest={}
for line in (ROOT/'provenance/final_raw_sha256.txt').read_text().splitlines():
 digest,path=line.split(maxsplit=1);manifest[path.strip()]=digest
for name in ids:
 for suffix in ['history.csv','diagnostics.csv','status.txt']:
  path=ROOT/'raw'/f'{name}_{suffix}';key='raw/'+path.name
  check(hashlib.sha256(path.read_bytes()).hexdigest()==manifest[key],key+' matches final remote raw SHA-256')
 check((ROOT/'results/histories'/f'{name}.csv.gz').exists(),name+' lossless derived history')
for name in ['legacy_r9','legacy_r10']:
 check('solver_stop' in (ROOT/'raw'/f'{name}_status.txt').read_text(),name+' failed evidence retained')
for key in range(1,9):
 matches=list((ROOT/'figures').glob(f'figure{key}_*.pdf'));check(len(matches)==1,f'Figure {key} PDF')
 for ext in ['png','svg']:check(matches[0].with_suffix('.'+ext).exists(),f'Figure {key} {ext}')
qa=json.loads((ROOT/'results/figure_qa.json').read_text());check(len(qa)==11 and all(q['collision']=='PASS' and q['text_exit']==0 for q in qa),'final 11-figure rendered QA')
root=ROOT.parent.parent
changed=subprocess.check_output(['git','diff','04f9bb2d856512531ebcb6265e4b46ef35cc43b0','--name-only','--','analysis/p2_h_double_moving_boundary_validation'],cwd=root,text=True)
check(not changed.strip(),'Task27 preserved without changes')
(ROOT/'results/delivery_audit.json').write_text(json.dumps(dict(passed=len(checks),checks=checks),indent=2))
print('PASS',len(checks),'delivery and provenance assertions')
