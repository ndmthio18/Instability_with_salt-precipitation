"""Test cumulative-budget comparisons on genuinely different accepted time grids."""
import json
from pathlib import Path
from copy import deepcopy
from analyze_matched_pair import compare_region


def series(times, initial, flux):
    return [dict(step=k,time_s=t,dt_s=t-times[k-1] if k else 0.,actual_z_over_H=.9,
        inventory_salt_m=initial-t*flux,inventory_initial_salt_m=initial,
        cut_salt_total_down_m_s=flux,other_salt_out_m_s=0.,
        cumulative_cut_salt_m=t*flux,cumulative_cut_adv_m=t*flux,
        cumulative_cut_diff_m=0.,salt_residual_m=0.,threshold_not_yet_observed=1.) for k,t in enumerate(times)]

a=series([0.,1.,2.,3.],.125,.0078125)
b=series([0.,.5,1.,1.5,2.,2.5,3.],.1875,.015625)
p=compare_region(a,b,common_accepted_only=True)
checks=dict(exact_shared_physical_times=[r['time_s'] for r in p]==[0.,1.,2.,3.],
    initial_offset_corrected=all(r['initial_delta_inventory_m']==.0625 for r in p),
    cumulative_flux_not_reintegrated=all(r['delta_C_m']==r['time_s']*.0078125 for r in p),
    identity_preserved=all(r['identity_error_m']==0 for r in p),
    differing_dt_retained=p[-1]['M0_accepted_dt_s']==1 and p[-1]['M1_accepted_dt_s']==.5)
try:compare_region(a,b)
except ValueError:checks['default_still_rejects_step_mismatch']=True
else:checks['default_still_rejects_step_mismatch']=False
bad=deepcopy(b)
for r in bad:r['actual_z_over_H']=.899
try:compare_region(a,bad,common_accepted_only=True)
except ValueError:checks['cut_shift_rejected']=True
else:checks['cut_shift_rejected']=False
try:compare_region(a,series([0.,.5,1.5,2.5],.1875,.015625),common_accepted_only=True)
except ValueError:checks['initial_only_overlap_rejected']=True
else:checks['initial_only_overlap_rejected']=False
bad=deepcopy(b)
for r in bad:
    if r['time_s']>=1.5:r['threshold_not_yet_observed']=0
checks['crossing_excluded_between_common_times']=[r['time_s'] for r in compare_region(a,bad,common_accepted_only=True)]==[0.,1.]
record=dict(passed=all(checks.values()),checks=checks,scope='exact shared accepted times only; no interpolation error estimate or convergence claim')
p=Path(__file__).resolve().parents[1]/'analysis/pilot0_integrated_20261102/common_time_checks.json'
p.write_text(json.dumps(record,indent=2)+'\n',encoding='utf-8')
print(json.dumps(record,indent=2));raise SystemExit(0 if record['passed'] else 1)
