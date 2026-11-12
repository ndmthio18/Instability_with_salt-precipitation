"""Run on a compute node after pilot array completion; fail closed."""
import csv,json
from pathlib import Path
R=Path.cwd();results={};valid=True
for case in csv.DictReader(open(R/'pilot_cases.csv')):
    name=case['id'];lines=list(csv.DictReader(open(R/'raw'/f'{name}_history.csv')))
    s=dict(x.split('=',1) for x in (R/'raw'/f'{name}_status.txt').read_text().splitlines())
    residual=max(abs(float(row[key])) for row in lines for key in ['salt_res','water_res','solid_res','liquid_salt_res','crust_salt_res','liquid_water_res','crust_water_res'])
    ok=s['status']=='end' and float(lines[-1]['t'])==3600 and residual<2e-9 and float(lines[-1]['zb'])>float(case['gap'])
    results[name]={'pass':ok,'budget':residual,'last':lines[-1]};valid &= ok
out={'pass':valid,'cases':results}
(R/'pilot_gate.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
assert valid,'P2 pilot gate failed; production prohibited'
