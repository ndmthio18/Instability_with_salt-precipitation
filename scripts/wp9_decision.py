"""WP9 decision: evaluate the pre-registered gates and record GO/NO-GO.

The gates come from decision_criteria.json, written before any R96 output was read.
Nothing here loosens a gate; a failed gate changes the recommendation, not the rule.
"""
import collections,json
from pathlib import Path
import numpy as np
import wp9_common as c


def series():
    by=collections.defaultdict(list)
    for r in c.rows(c.OUT/'resolution_series.csv'):by[r['leg']].append(r)
    return by


def at(by,leg,frame,key):
    if leg not in by or frame>=len(by[leg]):return None
    v=by[leg][frame][key]
    return float(v) if v not in ('',None) else None


def plane_value(leg,region,checkpoint,key):
    for r in c.rows(c.OUT/'resolution_seven_plane.csv'):
        if r['leg']==leg and r['region']==region and r['checkpoint']==checkpoint:
            return float(r[key])
    return None


def rel(a,b):
    if a is None or b is None or b==0:return None
    return abs(a-b)/abs(b)


def main():
    crit=json.loads((c.OUT/'decision_criteria.json').read_text(encoding='utf-8'))
    comp=json.loads((c.OUT/'resolution_comparison.json').read_text(encoding='utf-8'))
    lsa=json.loads((c.OUT/'lsa_audit.json').read_text(encoding='utf-8'))
    corr=json.loads((c.OUT/'corrections_v2.json').read_text(encoding='utf-8'))
    val=json.loads((c.OUT/'r96_validation.json').read_text(encoding='utf-8')) if (c.OUT/'r96_validation.json').exists() else None
    by=series()
    have96='R96_M1' in by
    F={'T':1200,'2T':2400,'3T':3600,'4T':4800}
    g={}

    if have96:
        g['G1']=dict(quantity='dominant box mode at 4T',
            R48=at(by,'R48_M1',F['4T'],'peak_mode'),R96=at(by,'R96_M1',F['4T'],'peak_mode'),
            passed=at(by,'R48_M1',F['4T'],'peak_mode')==at(by,'R96_M1',F['4T'],'peak_mode'))
        r2={k:rel(at(by,'R96_M1',F[k],'energy_m2'),at(by,'R48_M1',F[k],'energy_m2')) for k in F}
        g['G2']=dict(quantity='canonical m=2 energy',relative_differences=r2,
            passed=all(r2[k] is not None and r2[k]<=0.20 for k in ['T','2T','3T'])
                   and r2['4T'] is not None and r2['4T']<=0.50)
        r3=rel(at(by,'R96_M1',F['4T'],'energy_weighted_k_rad_m'),
               at(by,'R48_M1',F['4T'],'energy_weighted_k_rad_m'))
        g['G3']=dict(quantity='energy-weighted physical wavenumber at 4T',relative_difference=r3,
            R48=at(by,'R48_M1',F['4T'],'energy_weighted_k_rad_m'),
            R96=at(by,'R96_M1',F['4T'],'energy_weighted_k_rad_m'),
            passed=r3 is not None and r3<=0.20)
        a=at(by,'R96_M1',F['4T'],'mean_response_top');b=at(by,'R48_M1',F['4T'],'mean_response_top')
        floor=comp['uncertainty']['x_resolvability_floor']
        ratio=(a/b) if (a is not None and b) else None
        g['G4']=dict(quantity='horizontal-mean top response against the frozen strict background at 4T',
            R48=b,R96=a,ratio=ratio,resolvability_floor=floor,
            both_above_floor=bool(a is not None and b is not None
                                  and abs(a)>floor and abs(b)>floor),
            passed=bool(ratio is not None and 0.5<=abs(ratio)<=2.0
                        and abs(a)>floor and abs(b)>floor))
        p48=plane_value('R48_M1','upper_0.995000','4T','relative_to_initial')
        p96=plane_value('R96_M1','upper_0.995000','4T','relative_to_initial')
        r5=rel(p96,p48)
        g['G5']=dict(quantity='initial-corrected Delta C at the 1 mm plane at 4T',
            R48=p48,R96=p96,relative_difference=r5,
            same_sign=bool(p48 is not None and p96 is not None and p48*p96>0),
            passed=bool(r5 is not None and r5<=0.50 and p48 is not None and p96 is not None and p48*p96>0))
        modes_ok=all(g[k]['passed'] for k in ['G1','G2','G3'])
        all_ok=all(g[k]['passed'] for k in ['G1','G2','G3','G4','G5'])
        if all_ok:outcome='GO';rule=crit['decision_rule']['all_five_pass']
        elif modes_ok:outcome='CONDITIONAL';rule=crit['decision_rule']['G1_G2_G3_pass_only']
        else:outcome='NO-GO';rule=crit['decision_rule']['G1_G2_G3_not_all_pass']
        # the discriminator is taken from the comparison, which already knows whether the
        # R96 numerical floor was measured by a paired unperturbed leg
        disc=dict(comp['discriminator'])
        disc['ratio_R96_over_R48']=disc.get('R96_over_R48_m10')
        disc['measured_R48_over_R24']=crit['noise_discriminator']['measured_R48_over_R24']
        disc['registered_rule']=dict(noise_seeded_if=crit['noise_discriminator']['noise_seeded_if'],
            nonlinear_cascade_if=crit['noise_discriminator']['nonlinear_cascade_if'],
            neither_if=crit['noise_discriminator']['neither_if'])
    else:
        outcome='PENDING';rule='the R96 leg has not been collected';disc=None
        modes_ok=None;all_ok=None

    limitation=dict(
        horizontal_resolution=dict(
            tested=True,
            evidence='R24 to R48 to R96 at fixed nz and fixed box width',
            verdict=None if not have96 else ('not the limiter for modes' if modes_ok else 'still the limiter at 4T')),
        temporal_integration=dict(
            tested='partially',
            evidence=('R96 halves dt together with dx; the only bound is the WP6 R24 dt/2 pair to 2T, '
                      'which moved growth metrics by about 0.1 percent but moved the 1 mm Delta C by '
                      '68 percent at the shared mid-window time'),
            verdict='unseparated from horizontal resolution in this comparison'),
        vertical_resolution=dict(
            tested=False,
            evidence='nz=200 with Grading1=-0.98 is common to all three legs',
            verdict='untested; common-mode truncation cannot be excluded'),
        domain_width=dict(
            tested=False,
            evidence=('the box is one A28 wavelength wide with no-flow sides, so only Neumann modes '
                      'm=1,2,... exist; the LSA global critical mode nx=5 has wavelength 0.12 m and '
                      'cannot appear at any resolution'),
            verdict='a structural restriction of every leg, not a resolution effect'))

    out=dict(scope='command/16 WP9',criteria_sha256=c.sha(c.OUT/'decision_criteria.json'),
        r96_present=have96,
        r96_numerical_acceptance=val['passed'] if val else None,
        r96_censored_by_threshold=val['censored_by_threshold'] if val else None,
        gates=g,modes_consistent=modes_ok,all_gates_passed=all_ok,
        outcome=outcome,recommendation=rule,
        noise_discriminator=disc,
        lsa_context=dict(fastest_box_mode_at_4T=lsa['per_checkpoint']['4T']['fastest_box_mode_m'],
            fastest_box_mode_at_T=lsa['per_checkpoint']['T']['fastest_box_mode_m'],
            m10_unstable_after='between T and 2T',
            small_perturbation_agreement=json.loads((c.OUT/'lsa_versus_dns.json').read_text(encoding='utf-8'))['agreement']['max_absolute_relative_difference']),
        solubility=dict(status=corr['correction_2_censoring']['status'],
                        censoring_time_s=corr['correction_2_censoring']['censoring_time_s']),
        limitation_attribution=limitation,
        local_flux_mechanism_validation=('proposed' if outcome=='GO' else 'not proposed'),
        R96_M0=('completed under the separate user authorization WP9-A02; '
                'no further R96-M0 run is authorized'),
        R192='NO-GO',multiwave='NO-GO',full_domain='NO-GO',later_8T='NO-GO',
        t_eff=None,crystallization_delay_label=None,training_labels_created=False)
    c.save(c.OUT/'decision.json',out)
    print(json.dumps(dict(outcome=outcome,modes_consistent=modes_ok,
                          gates={k:v['passed'] for k,v in g.items()},
                          discriminator=disc),indent=2,default=str))


if __name__=='__main__':main()
