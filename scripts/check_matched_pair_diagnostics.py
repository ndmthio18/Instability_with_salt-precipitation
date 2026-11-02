"""Scientific regression checks for sign, initial offsets, censoring and guards."""
from copy import deepcopy
import json
from pathlib import Path
from analyze_matched_pair import compare_region


def series(initial, flux):
    return [dict(step=k,time_s=float(k),dt_s=float(k>0),actual_z_over_H=.9,
                 inventory_salt_m=initial-k*flux, inventory_initial_salt_m=initial,
                 cut_salt_total_down_m_s=flux, other_salt_out_m_s=0.,
                 cumulative_cut_salt_m=k*flux,cumulative_cut_adv_m=k*flux,
                 cumulative_cut_diff_m=0.,salt_residual_m=0.,threshold_not_yet_observed=1.)
            for k in range(3)]


a,b=series(.125,.0078125),series(.1875,.015625)
p=compare_region(a,b)
checks={
    'nonzero_initial_offset_removed':p[-1]['initial_delta_inventory_m']==.0625,
    'positive_extra_down_flux_negative_corrected_inventory':p[-1]['delta_C_m']==.015625 and p[-1]['corrected_delta_inventory_m']==-.015625,
    'exact_conservation':p[-1]['identity_error_m']==0,
    'self_null':all(r['delta_C_m']==r['corrected_delta_inventory_m']==0 for r in compare_region(a,a))}
for name,key,value in [('time_mismatch_rejected','time_s',1.1),
                       ('cut_mismatch_rejected','actual_z_over_H',.899),
                       ('corrupt_cumulative_rejected','cumulative_cut_salt_m',.3)]:
    bad=deepcopy(b);bad[1][key]=value
    try:compare_region(a,bad)
    except ValueError:checks[name]=True
    else:checks[name]=False
c=deepcopy(b);c[-1]['threshold_not_yet_observed']=0
checks['first_crossing_step_excluded']=len(compare_region(a,c))==2
out=Path(__file__).resolve().parents[1]/'analysis/matched_pair_36877_20261102/diagnostic_checks.json'
out.write_text(json.dumps({'checks':checks,'all_passed':all(checks.values()),
    'scope':'synthetic algebra/guards plus separate real M0 self-null; not M1 physical validation'},indent=2))
assert all(checks.values()),checks
print('PASS',len(checks),'scientific diagnostic checks')
