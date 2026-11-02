"""Pre-run tests of the resource gate and actual registered formal pair inputs."""
import ast
from copy import deepcopy
import json
from pathlib import Path
from analyze_matched_pair import validate_pair_inputs, validate_pair, rows, audit_series, compare_region, REGIONS
from pilot0_campaign import OUT, ROOT

node=next(n for n in ast.parse((ROOT/'scripts/pilot0_remote_worker.py').read_text(encoding='utf-8')).body if isinstance(n,ast.FunctionDef) and n.name=='eligible')
namespace={};exec(compile(ast.Module(body=[node],type_ignores=[]),'resource_gate','exec'),namespace)
eligible=namespace['eligible']
normal=dict(active_paper2_solvers=[],load_1m=1,load_5m=1,available_memory_bytes=8*1024**3,free_disk_bytes=100*1024**3)
checks={'free_resource_gate_accepts':eligible(normal)}
for name,key,value in [('existing_solver_blocks','active_paper2_solvers',[{'pid':123}]),
        ('high_load_blocks','load_1m',24),('sustained_load_blocks','load_5m',24),
        ('low_memory_blocks','available_memory_bytes',1024**3),('low_disk_blocks','free_disk_bytes',1024**3)]:
    case=deepcopy(normal);case[key]=value;checks[name]=not eligible(case)
v0,v1,plans=validate_pair_inputs(OUT/'runs/M0',OUT/'runs/M1')
checks['actual_registered_pair_inputs_pass_formal_interface']=True
try:validate_pair(OUT/'runs/M0',OUT/'runs/M1')
except FileNotFoundError:checks['unexecuted_M1_is_not_certified_as_complete']=True
else:checks['unexecuted_M1_is_not_certified_as_complete']=False
checks['same_grid_dt_end_and_threshold']=all(v0[k]==v1[k] for k in [
    'Grid.Cells0','Grid.Cells1','Grid.Grading1','TimeLoop.TEnd','TimeLoop.MaxTimeStepSize',
    'TimeLoop.DtInitial','Transport.StopAtSolubility','Transport.SolubilityMoleFraction'])
data=rows(OUT/'runs/M0/transport_budget.csv')
for region in REGIONS:
    series=[r for r in data if r['region']==region];audit_series(series)
    paired=compare_region(series,series)
    checks['M0_actual_budget_and_self_null/'+region]=len(paired)==601 and all(
        r['delta_C_m']==r['corrected_delta_inventory_m']==0 and
        r['M0_cumulative_total_down_m']==r['M1_cumulative_total_down_m'] for r in paired)
for path in ROOT.joinpath('scripts').glob('*.py'):ast.parse(path.read_text(encoding='utf-8'))
checks['all_python_sources_parse']=True
result=dict(passed=all(checks.values()),checks=checks,scope='pre-run guards, no claim of completed M1 or numerical validation')
(OUT/'execution_checks.json').write_text(json.dumps(result,indent=2)+'\n',encoding='utf-8')
print(json.dumps(result,indent=2))
raise SystemExit(0 if result['passed'] else 1)
