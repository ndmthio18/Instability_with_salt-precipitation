"""Async-safe reanalysis of frozen 2T pair; writes only the command/8 review."""
import argparse,csv,hashlib,json
from pathlib import Path
import numpy as np
from wp6_paired_stopping import validate_leg,paired_times,budget_rows_at_events
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT/'analysis/wp6_execution_20261103';OUT=ROOT/'analysis/wp6_review_20261103'
def rows(p):
    with p.open() as f:return list(csv.DictReader(f))
def write(name,v):(OUT/name).write_text(json.dumps(v,indent=2,allow_nan=False)+'\n')
def savecsv(name,rr):
    with (OUT/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rr[0]));w.writeheader();w.writerows(rr)
def longest_duration(rr,predicate):
    start=None;longest=0
    for r in rr:
        if predicate(r):
            if start is None:start=r['time_s']
            longest=max(longest,r['time_s']-start)
        else:start=None
    return longest


def main(directories=None):
    legs={};budget={}
    directories=directories or {role:BASE/'runs'/('WP6-'+role+'-extend') for role in ['M0','M1']}
    canonical=all(directories[role].resolve()==(BASE/'runs'/('WP6-'+role+'-extend')).resolve() for role in directories)
    assert OUT.resolve().is_relative_to((ROOT/'analysis').resolve()) and OUT.resolve()!=BASE.resolve()
    OUT.mkdir(parents=True,exist_ok=True)
    for role in ['M0','M1']:
        d=directories[role];f=json.loads((d/'finished.json').read_text())
        legs[role]=validate_leg(rows(d/'surface_events.csv'),f['solver']['stop'],f['solver']['returncode'])
        assert json.loads((d/'post/independent_checks.json').read_text())['passed']
        budget[role]=rows(d/'transport_budget.csv')
        assert all(abs(float(r['salt_residual_m']))<1e-14 and abs(float(r['water_residual_m']))<1e-11 for r in budget[role])
        assert all(abs(float(r['whole_salt_residual_m']))<1e-14 and abs(float(r['whole_water_residual_m']))<1e-11 for r in budget[role])
    match=paired_times(legs['M0'],legs['M1']);write('async_stopping_checks.json',dict(legs=legs,pair=match))
    cuts=json.loads((directories['M0']/'transport_metadata.json').read_text())['cuts']
    cuts1=json.loads((directories['M1']/'transport_metadata.json').read_text())['cuts']
    assert [(c['name'],c['actual_z_over_H']) for c in cuts]==[(c['name'],c['actual_z_over_H']) for c in cuts1]
    pairs=[];summaries=[]
    for cut in cuts:
        reg=cut['name'];a=budget_rows_at_events(budget['M0'],legs['M0'],reg);b=budget_rows_at_events(budget['M1'],legs['M1'],reg)
        init=float(b[0]['inventory_salt_m'])-float(a[0]['inventory_salt_m']);rr=[]
        assert abs(init)<1e-15, 'initial fields not a matched pair'
        for i,j in match['indices']:
            x,y=a[i],b[j];dc=float(y['cumulative_cut_salt_m'])-float(x['cumulative_cut_salt_m']);di=float(y['inventory_salt_m'])-float(x['inventory_salt_m'])-init;dr=float(y['salt_residual_m'])-float(x['salt_residual_m']);env=abs(float(y['salt_residual_m']))+abs(float(x['salt_residual_m']))
            q=dict(region=reg,time_s=float(x['time_s']),initial_difference_m=init,corrected_inventory_difference_m=di,delta_C_m=dc,relative_delta_C_initial_inventory=dc/float(x['inventory_initial_salt_m']),delta_C_adv_m=float(y['cumulative_cut_adv_m'])-float(x['cumulative_cut_adv_m']),delta_C_diff_m=float(y['cumulative_cut_diff_m'])-float(x['cumulative_cut_diff_m']),M0_signed_cut_total=float(x['cut_salt_total_down_m_s']),M1_signed_cut_total=float(y['cut_salt_total_down_m_s']),residual_envelope_m=env,identity_error_m=di+dc-dr,candidate_threshold_m=max(1e-13,10*env),instantaneous_candidate=abs(dc)>max(1e-13,10*env))
            for role,r in [('M0',x),('M1',y)]:
                for component in ['adv','diff']:
                    q[role+'_signed_cut_'+component]=float(r['cut_salt_'+component+'_down_m_s'])
            rr.append(q);pairs.append(q)
        if len(rr)<2:
            summaries.append(dict(region=reg,comparison_status='insufficient_common_prethreshold_states',transport_effectiveness='unidentified'))
            continue
        assert max(abs(r['identity_error_m']) for r in rr)<1e-18
        trigger=False
        for q in rr:
            start=q['time_s']-10*227.3717729539139
            if start<0:continue
            window=[v for v in rr if v['time_s']>=start-1e-7 and v['time_s']<=q['time_s']]
            if window[-1]['time_s']-window[0]['time_s']<10*227.3717729539139-1e-7:continue
            trigger|=any(sum(v['instantaneous_candidate'] and sign*v['delta_C_m']>0 for v in window)/len(window)>=.8 for sign in [-1,1])
        final=rr[-1];maxenv=max(r['residual_envelope_m'] for r in rr)
        summaries.append(dict(region=reg,depth_mm=(1-cut['actual_z_over_H'])*200,
            delta_C_end_m=final['delta_C_m'],relative_to_initial_inventory=final['relative_delta_C_initial_inventory'],delta_C_adv_end_m=final['delta_C_adv_m'],delta_C_diff_end_m=final['delta_C_diff_m'],
            max_residual_envelope_m=maxenv,end_signal_over_max_closure=abs(final['delta_C_m'])/maxenv if maxenv>0 else None,
            max_abs_identity_error_m=max(abs(r['identity_error_m']) for r in rr),
            positive_exceeds_closure_longest_s=longest_duration(rr,lambda r:r['delta_C_m']>r['residual_envelope_m']),
            negative_exceeds_closure_longest_s=longest_duration(rr,lambda r:-r['delta_C_m']>r['residual_envelope_m']),
            positive_sign_longest_s=longest_duration(rr,lambda r:r['delta_C_m']>0),
            fraction_postinitial_above_closure=sum(abs(r['delta_C_m'])>r['residual_envelope_m'] for r in rr[1:])/(len(rr)-1),
            frozen_G2_trigger=trigger,transport_effectiveness='unresolved: discretization and domain errors unknown',
            **{k:final[k] for k in final if k.startswith('M0_signed') or k.startswith('M1_signed')}))
    if canonical:
        old=rows(BASE/'extend_paired_transport.csv');assert len(old)==len(pairs)
        differences={}
        for x,y in zip(old,pairs):
            assert x['region']==y['region']
            for key in x:
                if key=='region':continue
                if key=='instantaneous_candidate':assert (x[key]=='True')==y[key];continue
                differences[key]=max(differences.get(key,0),abs(float(x[key])-y[key]))
        assert all(v==0 for v in differences.values()),differences
        write('frozen_2T_parity.json',dict(passed=True,rows=len(old),shared_states=len(match['indices']),maximum_absolute_differences=differences,old_pair_sha256=hashlib.sha256((BASE/'extend_paired_transport.csv').read_bytes()).hexdigest(),frozen_files_modified=False))
    if pairs:savecsv('paired_transport.csv',pairs)
    savecsv('seven_plane_review.csv',summaries[:7])
    write('transport_evidence.json',dict(surfaces=summaries,unknown_errors=['paired_U_dt','U_z','U_x','U_domain'],t_eff=None,
        inference='Above closure does not establish convergence. Absence of G2 does not prove zero transport.',
        sign_convention='positive down; positive delta_C can mean less upward transport when both total fluxes are negative'))
    print(json.dumps(dict(pair_rows=len(pairs),shared_states=len(match['indices']),frozen_parity_passed=True if canonical else None,surface_event_difference=match['surface_time_difference'])))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--m0',type=Path);p.add_argument('--m1',type=Path);p.add_argument('--out',type=Path);a=p.parse_args()
    if any(v is not None for v in [a.m0,a.m1,a.out]):
        if not all(v is not None for v in [a.m0,a.m1,a.out]):p.error('supply --m0, --m1 and --out together')
        if a.out.exists():p.error('custom output directory must be new; preserve previous results')
        OUT=a.out;main(dict(M0=a.m0,M1=a.m1))
    else:main()
