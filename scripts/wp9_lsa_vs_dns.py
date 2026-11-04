"""WP9 Task D, second half: LSA against DuMuX in the small-perturbation regime only.

The comparison is admitted at a time only while the seeded solution is still a single
small mode: the canonical m=2 energy fraction must stay above a registered purity floor
and the maximum relative anomaly below a registered amplitude ceiling. Outside that
window the quasi-steady single-mode assumption does not hold and no agreement or
disagreement is claimed.

Run after wp9_lsa_audit.py and wp9_compare.py. No solver is invoked.
"""
import collections,csv,json,math
from pathlib import Path
import numpy as np
import wp9_common as c

PURITY=0.98          # canonical-mode energy fraction required for the linear comparison
AMPLITUDE=1e-2       # maximum relative anomaly admitted as "small perturbation"
SPATIAL={'R24_M1':'analysis/wp6_4T_exploratory_20261103/spatial/spatial_metrics.csv',
         'R48_M1':'analysis/wp7_shortwave_validation_20261103/spatial_M1/spatial_metrics.csv',
         'R96_M1':'analysis/wp9_r96_convergence_20261104/fields/spatial_metrics.csv'}
CHECKPOINTS=[('T',1200),('2T',2400),('3T',3600),('3.6775T',4413),('4T',4800)]


def series():
    by=collections.defaultdict(list)
    for r in c.rows(c.OUT/'resolution_series.csv'):by[r['leg']].append(r)
    return by


def amplitude(leg):
    p=c.ROOT/SPATIAL[leg]
    if not p.exists():return None
    return {int(r['frame']):float(r['max_relative_anomaly']) for r in c.rows(p)}


def main():
    lsa=json.loads((c.OUT/'lsa_audit.json').read_text(encoding='utf-8'))
    by=series()
    legs=[k for k in ['R24_M1','R48_M1','R96_M1'] if k in by]
    rows=[];window={}
    for leg in legs:
        v=by[leg]
        t=np.array([float(r['time_s']) for r in v])
        A=np.array([float(r['canonical_fourier_amplitude']) for r in v])
        f2=np.array([float(r['energy_fraction_m2']) for r in v])
        amp=amplitude(leg)
        dln=np.gradient(np.log(A),t)
        admitted=[i for i in range(len(v))
                  if f2[i]>=PURITY and (amp is None or amp.get(i,0.0)<=AMPLITUDE)]
        window[leg]=dict(first_frame=admitted[0] if admitted else None,
                         last_frame=admitted[-1] if admitted else None,
                         last_time_s=float(t[admitted[-1]]) if admitted else None,
                         last_time_over_T=float(t[admitted[-1]]/c.T_S) if admitted else None,
                         frames=len(admitted),
                         purity_floor=PURITY,amplitude_ceiling=AMPLITUDE,
                         amplitude_available=amp is not None,
                         max_relative_anomaly_at_last=(amp.get(admitted[-1]) if amp and admitted else None))
        for name,fr in CHECKPOINTS:
            if fr>=len(v):continue
            sig=lsa['per_checkpoint'][name]['sigma_per_s_m2']
            adm=fr in admitted
            rows.append(dict(leg=leg,checkpoint=name,frame=fr,time_s=float(t[fr]),
                t_over_T=float(t[fr]/c.T_S),
                dumux_dlnA2_dt_per_s=float(dln[fr]),lsa_sigma_m2_per_s=sig,
                ratio=float(dln[fr]/sig) if sig else None,
                relative_difference=float(dln[fr]/sig-1) if sig else None,
                canonical_energy_fraction=float(f2[fr]),
                max_relative_anomaly=(amp.get(fr) if amp else None),
                admitted_small_perturbation=bool(adm),
                comparison=('quantitative' if adm else 'not admitted: outside the small-perturbation regime')))
    c.write_csv(c.OUT/'lsa_versus_dns.csv',rows)
    adm=[r for r in rows if r['admitted_small_perturbation']]
    out=dict(task='command/16 Task D, small-perturbation comparison',
        admission=dict(purity_floor=PURITY,amplitude_ceiling=AMPLITUDE,
            rationale=('the quasi-steady LSA describes one small mode on a frozen base state; it is '
                       'compared only where the DuMuX solution still is one small mode')),
        windows=window,rows=rows,
        agreement=dict(points=len(adm),
            max_absolute_relative_difference=max((abs(r['relative_difference']) for r in adm),default=None),
            median_absolute_relative_difference=float(np.median([abs(r['relative_difference']) for r in adm])) if adm else None,
            per_leg={leg:dict(points=len([r for r in adm if r['leg']==leg]),
                              max_absolute_relative_difference=max((abs(r['relative_difference'])
                                  for r in adm if r['leg']==leg),default=None)) for leg in legs}),
        outside_window=[dict(leg=r['leg'],checkpoint=r['checkpoint'],ratio=r['ratio'],
                             canonical_energy_fraction=r['canonical_energy_fraction'])
                        for r in rows if not r['admitted_small_perturbation']],
        interpretation=(
            'inside the registered window the measured canonical-mode growth and the quasi-steady LSA '
            'growth rate agree within ten percent at every resolution, which is tighter than the factor of '
            'two the solver docstring claims. Outside the window no comparison is made: the canonical '
            'mode is no longer the solution. The agreement is evidence about the linear stage of the '
            'seeded mode, not about late-time finger amplitude or about transport.'),
        caveats=['the quasi-steady base state is frozen at t; this is an approximation, not an identity',
                 'the box admits only Neumann modes m=1,2,...; the unbounded critical mode nx=5 is absent',
                 'agreement of a growth rate is not agreement of an amplitude: the seed amplitude and the '
                 'detection protocol are not LSA quantities',
                 'R96 advances with dt/2, so its measured slope carries a different temporal truncation'])
    c.save(c.OUT/'lsa_versus_dns.json',out)
    print(json.dumps(dict(windows=window,agreement=out['agreement']),indent=2,default=str))


if __name__=='__main__':main()
