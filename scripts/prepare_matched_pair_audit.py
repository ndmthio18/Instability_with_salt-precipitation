"""Freeze a local M1 proposal and audit existing evidence; never execute it."""
import csv
import hashlib
import json
import shutil
from pathlib import Path
from prepare_transport_run import parse
from analyze_matched_pair import rows, audit_series, REGIONS

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT/'analysis/matched_pair_36877_20261102'
M0 = ROOT/'analysis/unperturbed_36877_20261102'


def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()


def main():
    remote = json.loads((OUT/'remote_audit.json').read_text())
    old = json.loads((M0/'plan.json').read_text())
    checks = []

    def check(name, passed):
        checks.append(dict(name=name,passed=bool(passed)))

    check('remote_executable_matches_M0',remote['executable']['sha256']==old['executable']['sha256'])
    for name,h in old['source_sha256'].items():
        check('source/'+name, sha(ROOT/'model/transport_v1/source'/name)==h==remote['source_hashes'][name])
    frozen = json.loads((ROOT/'analysis/phase0_20261102/pilot_candidates.json').read_text())[0]
    for name,h in frozen['existing_source_hashes'].items():
        current = next(v['sha256'] for p,v in remote['files'].items() if p.endswith('/'+name)) if not name.endswith('.pvd') else None
        if current:
            check('legacy_frozen/'+name,h==current)
    check('M0_actual_remote_input_hash',sha(OUT/'M0_actual_params.input')==old['input_sha256'])
    check('M0_local_input_only_line_ending_difference',
          (M0/'params.input').read_bytes().replace(b'\r\n',b'\n')==(OUT/'M0_actual_params.input').read_bytes())
    finished=json.loads((M0/'evidence/finished.json').read_text())
    for name in ('transport_budget.csv','surface_events.csv','transport_metadata.json','stop_record.json'):
        check('M0_output/'+name,sha(M0/'evidence'/name)==finished['files'][name]['sha256'])
    check('M0_normal_completion',finished['returncode']==0 and finished['stop']['accepted_steps']==600)
    check('legacy_initial_formula',remote['initial_field']['legacy_vs_formula_max_x_error']<1e-16)
    check('same_initial_pressure',remote['initial_field']['legacy_vs_M0_max_pressure_error_Pa']==0)
    check('same_initial_saturation',remote['initial_field']['legacy_vs_M0_max_Sw_error']==0)
    check('profile_normalized',remote['profile']['A_max']==1 and remote['profile']['z_min']==0 and remote['profile']['z_max']==.2)

    v0=parse((M0/'params.input').read_text())
    legacy=parse(next((OUT/'legacy_nx28').glob('*.input')).read_text())
    v1=v0.copy()
    v1.update({'Problem.Name':'matched_36877_nx28', 'Problem.PerturbationType':'Eigenprofile',
               'Problem.Amplitude':legacy['Problem.Amplitude'],
               'Problem.StandardDeviationSaltMoleFraction':legacy['Problem.StandardDeviationSaltMoleFraction'],
               'Problem.EigenprofilePath':'eigenprofile.dat'})
    check('canonical_mode',v1['Problem.WaveNumber']=='28' and float(v1['Problem.Amplitude'])==1e-6)
    target=OUT/'planned_M1';target.mkdir(exist_ok=True)
    sections={}
    for key,val in v1.items():
        s,k=key.split('.',1);sections.setdefault(s,{})[k]=val
    inp='\n\n'.join('['+s+']\n'+'\n'.join(k+' = '+v for k,v in d.items()) for s,d in sections.items())+'\n'
    (target/'params.input').write_bytes(inp.encode())
    profile=next((OUT/'legacy_nx28').glob('*_eigenprofile.dat'))
    shutil.copyfile(profile,target/'eigenprofile.dat')
    proposed=dict(status='prepared_not_authorized_not_executed', solver_executed=False,
        case_id=36877,role='M1',executable=old['executable'],source_sha256=old['source_sha256'],
        parameters=v1,profile_sha256=sha(profile),input_sha256=sha(target/'params.input'),
        M0_plan_sha256=sha(M0/'plan.json'), M0_input_sha256=sha(M0/'params.input'),
        initialization_time_s=0, wall_timeout_s_proposed=43200,
        expected_accepted_steps_if_no_cuts=600,
        launch_gate='User approval and remote executable/source/input hashes; retain actual accepted dt and initial field')
    (target/'plan.json').write_text(json.dumps(proposed,indent=2),encoding='utf-8')
    with (OUT/'input_comparison.csv').open('w',newline='',encoding='utf-8-sig') as f:
        w=csv.writer(f);w.writerow(['parameter','M0','planned_M1','legacy_nx28','pair_status'])
        for k in sorted(set(v0)|set(v1)|set(legacy)):
            w.writerow([k,v0.get(k,''),v1.get(k,''),legacy.get(k,''),
                        'same' if v0.get(k)==v1.get(k) else 'intentional_initialization_or_name'])
    budget=rows(M0/'evidence/transport_budget.csv')
    summary=[]
    for region in REGIONS:
        series=[r for r in budget if r['region']==region]
        audit_series(series)
        r0=series[0];last=series[-1]
        summary.append(dict(region=region,actual_z_over_H=r0['actual_z_over_H'],
            initial_inventory_m=r0['inventory_salt_m'],final_inventory_m=last['inventory_salt_m'],
            final_cumulative_cut_m=last['cumulative_cut_salt_m'],
            max_abs_budget_residual_m=max(abs(r['salt_residual_m']) for r in series),
            max_relative_budget_residual=max(abs(r['salt_residual_m'])/r0['inventory_salt_m'] for r in series)))
    check('M0_accepted_time_steps',len({r['time_s'] for r in budget})==601)
    check('M0_no_observed_threshold',all(r['threshold_not_yet_observed']==1 for r in budget))
    check('M0_initial_inventory_reconstruction',max(abs(r['M0_csv_reconstruction_error_m']) for r in remote['initial_checks'])<2e-16)
    result=dict(checks=checks,all_passed=all(r['passed'] for r in checks),
        M0_local_raw_input_sha256=sha(M0/'params.input'),
        M0_actual_remote_input_sha256=sha(OUT/'M0_actual_params.input'),
        M0_local_input_caveat='Local initial writer used CRLF; remote submitted LF. Actual remote bytes separately archived; parameters identical.',
        M0_budgets=summary, M1_available=False, M1_initial_audit='predicted_from_formula_plus_legacy_t0_not_new_M1_execution',
        readiness='conditional_GO_for_one_approved_M1_only',physical_feasibility='NO_GO_pending_pair_and_sensitivity',
        current_dt_min_s=min(r['dt_s'] for r in budget if r['step']>0),
        current_dt_max_s=max(r['dt_s'] for r in budget),
        cost=dict(M0_observed_wall_s=finished['wall_s'],legacy_M1_observed_wall_s=19122,
                  proposed_M1_wall_hours_range=[1,6],proposed_cap_hours=12,
                  range_is_estimate_not_measured=True))
    (OUT/'audit_results.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    paths=[ROOT/'command/1.md',ROOT/'paper1_sources.lock.json',M0/'plan.json',M0/'params.input',
           M0/'evidence/finished.json',OUT/'remote_audit.json',OUT/'M0_actual_params.input',target/'params.input',target/'eigenprofile.dat',target/'plan.json']
    manifest={str(p.relative_to(ROOT)).replace('\\','/'):sha(p) for p in paths}
    (OUT/'source_manifest.json').write_text(json.dumps(manifest,indent=2),encoding='utf-8')
    print(json.dumps(result,indent=2))
    if not result['all_passed']:
        raise SystemExit('Audit failed')


if __name__=='__main__':
    main()
