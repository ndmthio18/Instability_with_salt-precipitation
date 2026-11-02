"""Conservative matched-pair CSV diagnostics. Never launches a solver.

Defaults to identical accepted time grids and actual BOX cuts. An explicit option
uses only exactly shared accepted times without interpolating cumulative budgets.
Diagnostic-only mode is for self-null checks, not physical pair certification.
"""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
from prepare_transport_run import parse

REGIONS = ('upper_0.900000', 'upper_0.750000', 'upper_0.500000')


def rows(path):
    with path.open(newline='') as f:
        result = [{k: v if k == 'region' else float(v) for k, v in r.items()}
                  for r in csv.DictReader(f)]
    if not result or not all(math.isfinite(v) for r in result for v in r.values()
                             if isinstance(v, float)):
        raise ValueError('Empty/nonfinite data: '+str(path))
    return result


def audit_series(series):
    if not series or series[0]['time_s'] != 0 or series[0]['dt_s'] != 0:
        raise ValueError('Initial row is required')
    initial = series[0]['inventory_salt_m']
    cumulative = 0.0
    previous = 0.0
    for k, r in enumerate(series):
        t, dt = r['time_s'], r['dt_s']
        if r['step'] != k or (k and (dt <= 0 or abs(t-previous-dt) > 1e-8)):
            raise ValueError('Invalid accepted-step history')
        if r['other_salt_out_m_s'] != 0:
            raise ValueError('Zero top/side salt flux required')
        cumulative += dt*r['cut_salt_total_down_m_s']
        residual = r['inventory_salt_m']-initial+cumulative
        for a, b in [(initial, r['inventory_initial_salt_m']),
                     (cumulative, r['cumulative_cut_salt_m']),
                     (residual, r['salt_residual_m']),
                     (r['cumulative_cut_adv_m']+r['cumulative_cut_diff_m'], cumulative)]:
            if abs(a-b) > 1e-18+1e-12*max(abs(a), abs(b)):
                raise ValueError('Inconsistent cumulative/initial budget fields')
        previous = t


def accepted_pairs(m0, m1, common_accepted_only=False):
    """Exact physical-time intersection; never interpolate or match by step index."""
    if not common_accepted_only:
        yield from zip(m0,m1)
        return
    i=j=0
    while i<len(m0) and j<len(m1):
        a,b=m0[i],m1[j]
        if a['time_s']==b['time_s']:
            yield a,b;i+=1;j+=1
        elif a['time_s']<b['time_s']:i+=1
        else:j+=1


def compare_region(m0, m1, common_accepted_only=False):
    audit_series(m0)
    audit_series(m1)
    initial_delta = m1[0]['inventory_salt_m']-m0[0]['inventory_salt_m']
    reference = m0[0]['inventory_salt_m']
    if reference <= 0:
        raise ValueError('Positive fixed M0 initial inventory required')
    result = []
    for a, b in accepted_pairs(m0, m1, common_accepted_only):
        required=[('time_s',0 if common_accepted_only else 1e-8),('actual_z_over_H',1e-12)]
        if not common_accepted_only:required.append(('dt_s',1e-8))
        if any(abs(a[k]-b[k]) > tol for k, tol in
               required):
            raise ValueError('Unmatched accepted time grids or actual cuts')
        # Exclude the first crossing step of either run, including interior events.
        if not (a['threshold_not_yet_observed'] and b['threshold_not_yet_observed']):
            break
        dc = b['cumulative_cut_salt_m']-a['cumulative_cut_salt_m']
        di = b['inventory_salt_m']-a['inventory_salt_m']-initial_delta
        dr = b['salt_residual_m']-a['salt_residual_m']
        algebra = di+dc-dr
        if abs(algebra) > 1e-18:
            raise ValueError('Pair conservation identity failed')
        result.append(dict(time_s=a['time_s'], actual_z_over_H=a['actual_z_over_H'],
            M0_time_s=a['time_s'],M1_time_s=b['time_s'],
            M0_accepted_dt_s=a['dt_s'],M1_accepted_dt_s=b['dt_s'],
            M0_cumulative_total_down_m=a['cumulative_cut_salt_m'],
            M1_cumulative_total_down_m=b['cumulative_cut_salt_m'],
            M0_cumulative_adv_down_m=a['cumulative_cut_adv_m'],
            M1_cumulative_adv_down_m=b['cumulative_cut_adv_m'],
            M0_cumulative_diff_down_m=a['cumulative_cut_diff_m'],
            M1_cumulative_diff_down_m=b['cumulative_cut_diff_m'],
            M0_inventory_salt_m=a['inventory_salt_m'],M1_inventory_salt_m=b['inventory_salt_m'],
            M0_total_down_m_s=a['cut_salt_total_down_m_s'],M1_total_down_m_s=b['cut_salt_total_down_m_s'],
            delta_C_m=dc, corrected_delta_inventory_m=di,
            initial_delta_inventory_m=initial_delta, reference_inventory_m=reference,
            delta_C_over_reference=dc/reference,
            delta_adv_C_m=b['cumulative_cut_adv_m']-a['cumulative_cut_adv_m'],
            delta_diff_C_m=b['cumulative_cut_diff_m']-a['cumulative_cut_diff_m'],
            pair_residual_m=dr, identity_error_m=algebra,
            budget_envelope_m=abs(a['salt_residual_m'])+abs(b['salt_residual_m'])))
    if not result:
        raise ValueError('No common pre-threshold interval')
    if common_accepted_only and not any(r['time_s']>0 for r in result):
        raise ValueError('No positive common accepted time; do not infer a time history')
    return result


def validate_pair_inputs(m0, m1):
    """Validate frozen pairing before launch; completion checks are separate."""
    v0, v1 = [parse((p/'params.input').read_text()) for p in (m0, m1)]
    allowed = {'Problem.Name', 'Problem.PerturbationType', 'Problem.Amplitude',
               'Problem.StandardDeviationSaltMoleFraction', 'Problem.EigenprofilePath'}
    for k in set(v0) | set(v1):
        if k not in allowed and v0.get(k) != v1.get(k):
            raise ValueError('Unmatched input: '+k)
    if (v0['Problem.PerturbationType'] != 'Constant' or float(v0['Problem.Amplitude']) != 0
            or float(v0['Problem.StandardDeviationSaltMoleFraction']) != 0
            or v1['Problem.PerturbationType'] != 'Eigenprofile'
            or float(v1['Problem.Amplitude']) != 1e-6 or float(v1['Problem.WaveNumber']) != 28):
        raise ValueError('Requires M0 uniform and M1 canonical nx28 amplitude 1e-6')
    plans = [json.loads((p/'plan.json').read_text()) for p in (m0,m1)]
    if any(int(plan['case_id']) != 36877 for plan in plans):
        raise ValueError('This frozen protocol is for case 36877')
    if plans[0]['executable']['sha256'] != plans[1]['executable']['sha256']:
        raise ValueError('Different executables')
    if plans[0]['source_sha256'] != plans[1]['source_sha256']:
        raise ValueError('Different sources')
    for p,plan in zip((m0,m1),plans):
        if hashlib.sha256((p/'params.input').read_bytes()).hexdigest()!=plan['input_sha256']:
            raise ValueError('Input hash mismatch')
    profile=m1/v1['Problem.EigenprofilePath']
    if hashlib.sha256(profile.read_bytes()).hexdigest()!=plans[1]['profile_sha256'] or plans[1]['profile_sha256']!='d721e605d126e555f92b43200aa73b999ed4e2b61397edd1e98c0d9fddb36575':
        raise ValueError('Profile hash mismatch')
    return v0,v1,plans


def validate_pair(m0, m1):
    v0,v1,plans=validate_pair_inputs(m0,m1)
    for p, plan in zip((m0,m1), plans):
        if hashlib.sha256((p/'params.input').read_bytes()).hexdigest() != plan['input_sha256']:
            raise ValueError('Input hash mismatch')
        stop = json.loads((p/'stop_record.json').read_text())
        finished = json.loads((p/'finished.json').read_text())
        if finished['returncode'] != 0 or finished.get('wrapper_stop_reason') or stop['reason'] not in (
                'requested_horizon','surface_solubility_threshold','interior_solubility_threshold'):
            raise ValueError('Incomplete/failed run; retain separately')
        for name in ('transport_budget.csv','surface_events.csv','stop_record.json'):
            if hashlib.sha256((p/name).read_bytes()).hexdigest() != finished['files'][name]['sha256']:
                raise ValueError('Output hash mismatch: '+name)
        events=rows(p/'surface_events.csv')
        if len(events)!=stop['accepted_steps']+1 or events[-1]['time_s']!=stop['last_accepted_time_s']:
            raise ValueError('Incomplete event history')
        budget=rows(p/'transport_budget.csv')
        for region in REGIONS:
            br=[r for r in budget if r['region']==region]
            if len(br)!=len(events) or any(r['time_s']!=e['time_s'] for r,e in zip(br,events)):
                raise ValueError('Incomplete or mismatched regional budget')
    profile = m1/v1['Problem.EigenprofilePath']
    if hashlib.sha256(profile.read_bytes()).hexdigest() != plans[1]['profile_sha256'] or plans[1]['profile_sha256'] != 'd721e605d126e555f92b43200aa73b999ed4e2b61397edd1e98c0d9fddb36575':
        raise ValueError('Profile hash mismatch')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    for name in ('m0','m1','out'):
        ap.add_argument('--'+name, type=Path, required=True)
    ap.add_argument('--diagnostic-only', action='store_true')
    ap.add_argument('--common-accepted-times', action='store_true',
                    help='Use exact common physical times only; no interpolation')
    args = ap.parse_args()
    if not args.diagnostic_only:
        validate_pair(args.m0,args.m1)
    data = [rows(p/'transport_budget.csv') for p in (args.m0,args.m1)]
    events = [{r['time_s']:r for r in rows(p/'surface_events.csv')} for p in (args.m0,args.m1)]
    output = []
    summary = {}
    for region in REGIONS:
        pair = compare_region(*[[r for r in d if r['region']==region] for d in data],
                              common_accepted_only=args.common_accepted_times)
        for r in pair:
            r['region'] = region
            for label, ev in zip(('M0','M1'),events):
                e = ev[r[label+'_time_s']]
                for key in ('x_top_area_mean','x_top_max','sigma_top_node','x_domain_max'):
                    r[label+'_'+key] = e[key]
        output.extend(pair)
        summary[region] = dict(samples=len(pair), last=pair[-1],
            max_abs_delta_C_m=max(abs(r['delta_C_m']) for r in pair),
            max_identity_error_m=max(abs(r['identity_error_m']) for r in pair),
            max_budget_envelope_m=max(r['budget_envelope_m'] for r in pair))
    args.out.mkdir(parents=True, exist_ok=True)
    with (args.out/'pair_timeseries.csv').open('w', newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(output[0]));w.writeheader();w.writerows(output)
    record = dict(status='diagnostic_only_not_physical_evidence' if args.diagnostic_only else 'paired_diagnostics_pending_sensitivity',
        sensitivity_error_available=False, t_eff_defined=False,
        comparison_times='exact_common_accepted' if args.common_accepted_times else 'matched_accepted_grid',
        interpolation_used=False,
        alignment_error_at_reported_exact_times_m=0 if args.common_accepted_times else None,
        alignment_scope='Exact cumulative values at shared times; no claim about skipped times or temporal discretization error',
        transport_t_eff_s=None, transport_t_eff_status='undefined_not_a_label', regions=summary,
        sources={str(p/name):hashlib.sha256((p/name).read_bytes()).hexdigest()
                 for p in (args.m0,args.m1) for name in ('transport_budget.csv','surface_events.csv')})
    (args.out/'summary.json').write_text(json.dumps(record,indent=2))
    print(record['status'])


if __name__ == '__main__':
    main()
