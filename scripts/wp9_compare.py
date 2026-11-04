"""WP9 Task C: compare R24, R48 and R96 at shared accepted physical times.

Every comparison is made on exactly shared accepted times. R24 and R48 advance with the
campaign step dt and write one frame per step; R96 advances with dt/2 and writes one
frame per two steps, so frame j is the same mathematical time in all three and the only
difference is the double accumulation, which is bounded by the registered budget.

Runs with whatever resolutions are present, so it can be executed before and after the
R96 leg. No solver is invoked.
"""
import json,math
from pathlib import Path
import numpy as np
import wp9_common as c

BGRUN=c.ROOT/'analysis/wp8_full4T_20261103/run'
BGPROF=c.ROOT/'tmp/wp8_full4T_20261103/run/post/profiles.npz'
ZERO=c.ROOT/'tmp/wp8_zero_mode_20261103'
R96RUN=c.OUT/'run'
CHECKPOINTS=[('T',1200),('2T',2400),('3T',3600),('3.6775T',4413),('4T',4800)]
MODES=(2,6,10)
# measured round-off floors in the salt mole fraction, from command/15 and from the
# cross-executable background bound recomputed in wp9_corrections.py
X_RESOLVABILITY=2.458970527197124e-15
DT_GROWTH_PRIOR=0.0011457   # |relative change| of the weighted surface sigma under dt/2 at R24 to 2T

LEGS={'R24_M0':dict(nx=24,run='analysis/wp6_4T_exploratory_20261103/runs/WP6-M0-explore4T',
                    npz=ZERO/'R24_M0/profiles_spectra.npz',stride=1),
      'R24_M1':dict(nx=24,run='analysis/wp6_4T_exploratory_20261103/runs/WP6-M1-explore4T',
                    npz=ZERO/'R24_M1/profiles_spectra.npz',stride=1),
      'R48_M0':dict(nx=48,run='analysis/wp7_shortwave_validation_20261103/runs/WP7-R48-M0-4T',
                    npz=ZERO/'R48_M0/profiles_spectra.npz',stride=1),
      'R48_M1':dict(nx=48,run='analysis/wp7_shortwave_validation_20261103/runs/WP7-R48-M1-4T',
                    npz=ZERO/'R48_M1/profiles_spectra.npz',stride=1),
      'R96_M0':dict(nx=96,run='analysis/wp9_r96_convergence_20261104/r96_m0/run',
                    npz=c.CACHE/'R96_M0_profiles_spectra.npz',stride=2),
      'R96_M1':dict(nx=96,run='analysis/wp9_r96_convergence_20261104/run',
                    npz=c.CACHE/'R96_M1_profiles_spectra.npz',stride=2)}


def num(v):
    """Empty diagnostic cells are genuinely undefined (e.g. the projection of an
    unperturbed leg); they are kept as None rather than coerced to zero."""
    return float(v) if v not in ('',None) else None


def available():
    out={}
    for k,v in LEGS.items():
        run=c.ROOT/v['run']
        if Path(v['npz']).exists() and (run/'surface_events.csv').exists():
            out[k]=dict(v,run=run)
    return out


def load(leg):
    # materialise once: indexing an NpzFile re-decompresses the whole array every time
    with np.load(leg['npz']) as f:
        z={k:np.asarray(f[k]) for k in ['time_s','energy','k_rad_m','top_coeff','x_mean']}
    ev=c.rows(leg['run']/'surface_events.csv')
    on=c.rows(leg['run']/'weighted_online.csv')
    s=leg['stride']
    ev=ev[::s];on=on[::s]
    n=min(len(ev),len(on),z['time_s'].shape[0])
    return dict(z=z,ev=ev[:n],on=on[:n],n=n,nx=leg['nx'],stride=s,run=leg['run'])


def budget_rows(path,stride):
    """Seven-plane budget indexed by frame, keeping the registered regions only."""
    import collections
    reg=collections.defaultdict(list)
    for r in c.rows(path):reg[r['region']].append(r)
    return {k:v[::stride] for k,v in reg.items()}


def main():
    legs={k:load(v) for k,v in available().items()}
    assert 'R24_M1' in legs and 'R48_M1' in legs,'the R24 and R48 legs are required'
    with np.load(BGPROF) as f:bg={'x_mean':np.asarray(f['x_mean'])}
    bgev=c.rows(BGRUN/'surface_events.csv')
    bgbud=budget_rows(BGRUN/'transport_budget.csv',1)
    nshared=min(min(v['n'] for v in legs.values()),bg['x_mean'].shape[0],len(bgev))

    # ---- time-axis evidence ---------------------------------------------------------
    ref=legs['R48_M1']['z']['time_s'][:nshared]
    axis={}
    for k,v in legs.items():
        d=np.abs(v['z']['time_s'][:nshared]-ref)
        axis[k]=dict(max_time_difference_s=float(d.max()),
                     bitwise_identical=int(np.sum(d==0.0)),states=int(nshared),
                     budget_s=9600*math.ulp(c.HORIZON_S),
                     within_budget=bool(d.max()<=9600*math.ulp(c.HORIZON_S)))
    bgd=np.abs(np.array([float(r['time_s']) for r in bgev[:nshared]])-ref)
    axis['background']=dict(max_time_difference_s=float(bgd.max()),
                            bitwise_identical=int(np.sum(bgd==0.0)),states=int(nshared),
                            budget_s=9600*math.ulp(c.HORIZON_S),
                            within_budget=bool(bgd.max()<=9600*math.ulp(c.HORIZON_S)))

    # ---- per-frame series -----------------------------------------------------------
    series=[];spectra_rows=[]
    for name,v in legs.items():
        z=v['z'];en=z['energy'][:nshared];k=z['k_rad_m'];W=c.WIDTH_M
        tot=en.sum(axis=1)
        for i in range(nshared):
            e=en[i];t=float(ref[i]);total=float(tot[i])
            frac=lambda m:float(e[m-1]/total) if total>0 and m<=len(e) else None
            peak=int(np.argmax(e))+1 if total>0 else None
            series.append(dict(leg=name,nx_cells=v['nx'],frame=i,time_s=t,t_over_T=t/c.T_S,
                total_mode_energy=total,
                **{f'energy_m{m}':float(e[m-1]) if m<=len(e) else None for m in MODES},
                **{f'energy_fraction_m{m}':frac(m) for m in MODES},
                peak_mode=peak,peak_k_rad_m=float(peak*math.pi/W) if peak else None,
                peak_cells_per_wavelength=float(2*v['nx']/peak) if peak else None,
                energy_weighted_k_rad_m=float(k@e/total) if total>0 else None,
                canonical_fourier_amplitude=abs(float(z['top_coeff'][i][1])),
                weighted_surface_sigma=num(v['on'][i]['weighted_top_sigma']),
                weighted_fourier=num(v['on'][i]['weighted_fourier']),
                a28_projection=num(v['on'][i]['weighted_projection']),
                surface_mean=num(v['ev'][i]['x_top_area_mean']),
                surface_local_max=num(v['ev'][i]['x_top_max']),
                domain_local_max=num(v['ev'][i]['x_domain_max']),
                mean_response_top=float(z['x_mean'][i][-1]-bg['x_mean'][i][-1]),
                mean_response_Linf=float(np.max(np.abs(z['x_mean'][i]-bg['x_mean'][i])))))
        for label,fr in CHECKPOINTS:
            if fr>=nshared:continue
            e=en[fr];total=float(e.sum())
            for m in range(1,len(e)+1):
                spectra_rows.append(dict(leg=name,nx_cells=v['nx'],checkpoint=label,frame=fr,
                    time_s=float(ref[fr]),mode_m=m,wavelength_m=2*W/m,k_rad_m=float(k[m-1]),
                    cells_per_wavelength=2*v['nx']/m,energy=float(e[m-1]),
                    energy_fraction=float(e[m-1]/total) if total>0 else None,
                    resolvable=2*v['nx']/m>=8))
    c.write_csv(c.OUT/'resolution_series.csv',series)
    c.write_csv(c.OUT/'resolution_spectra.csv',spectra_rows)

    # ---- checkpoint comparison with uncertainty -------------------------------------
    idx={(r['leg'],r['frame']):r for r in series}
    pairs=[('R48_M1','R24_M1'),('R96_M1','R48_M1'),('R96_M1','R24_M1')]
    fields=['total_mode_energy','energy_m2','energy_m6','energy_m10','canonical_fourier_amplitude',
            'weighted_surface_sigma','a28_projection','surface_mean','surface_local_max',
            'peak_k_rad_m','energy_weighted_k_rad_m','mean_response_top','mean_response_Linf']
    comp=[]
    for label,fr in CHECKPOINTS:
        if fr>=nshared:continue
        for hi,lo in pairs:
            if (hi,fr) not in idx or (lo,fr) not in idx:continue
            a,b=idx[(hi,fr)],idx[(lo,fr)]
            floor={}
            for leg0,nx in [(hi,a['nx_cells']),(lo,b['nx_cells'])]:
                mk=leg0.split('_')[0]+'_M0'
                floor[leg0]=idx[(mk,fr)]['total_mode_energy'] if (mk,fr) in idx else None
            for f in fields:
                x,y=a.get(f),b.get(f)
                if x is None or y is None:continue
                comp.append(dict(checkpoint=label,frame=fr,time_s=a['time_s'],
                    quantity=f,fine=hi,coarse=lo,fine_value=x,coarse_value=y,
                    absolute_difference=x-y,
                    relative_difference=(x-y)/y if y else None,
                    ratio=x/y if y else None,
                    fine_unperturbed_floor=floor[hi],coarse_unperturbed_floor=floor[lo],
                    above_unperturbed_floor=(None if floor[lo] is None or not f.startswith(('energy','total'))
                                             else bool(min(abs(x),abs(y))>floor[lo])),
                    x_resolvability_floor=(X_RESOLVABILITY if f.startswith(('mean_response','surface')) else None),
                    timestep_prior=(DT_GROWTH_PRIOR if f in ('weighted_surface_sigma','canonical_fourier_amplitude',
                                                             'a28_projection') else None)))
    c.write_csv(c.OUT/'resolution_checkpoint_comparison.csv',comp)

    # ---- seven-plane inventories and signed cumulative transport ---------------------
    plane=[]
    for name,v in legs.items():
        bud=budget_rows(v['run']/'transport_budget.csv',v['stride'])
        for region,rowsr in bud.items():
            if region not in bgbud:continue
            for label,fr in CHECKPOINTS:
                if fr>=min(nshared,len(rowsr)):continue
                r=rowsr[fr];r0=rowsr[0];bgr=bgbud[region][fr];bg0=bgbud[region][0]
                dI=float(r['inventory_salt_m'])-float(r0['inventory_salt_m'])
                dIb=float(bgr['inventory_salt_m'])-float(bg0['inventory_salt_m'])
                plane.append(dict(leg=name,nx_cells=v['nx'],region=region,checkpoint=label,frame=fr,
                    time_s=float(r['time_s']),
                    inventory_salt_m=float(r['inventory_salt_m']),
                    initial_corrected_delta_C_m=dI-dIb,
                    relative_to_initial=(dI-dIb)/float(r0['inventory_salt_m']),
                    cumulative_cut_salt_m=float(r['cumulative_cut_salt_m']),
                    cumulative_cut_adv_m=float(r['cumulative_cut_adv_m']),
                    cumulative_cut_diff_m=float(r['cumulative_cut_diff_m']),
                    cut_salt_total_down_m_s=float(r['cut_salt_total_down_m_s']),
                    salt_residual_m=float(r['salt_residual_m']),
                    water_residual_m=float(r['water_residual_m']),
                    background_delta_C_m=dIb))
    c.write_csv(c.OUT/'resolution_seven_plane.csv',plane)

    # ---- the registered discriminator ------------------------------------------------
    final=4800 if 4800<nshared else nshared-1
    g=lambda leg,f:idx[(leg,final)][f] if (leg,final) in idx else None
    disc=dict(frame=final,time_s=float(ref[final]),
        m10_energy={k:g(k,'energy_m10') for k in legs},
        m2_energy_fraction={k:g(k,'energy_fraction_m2') for k in legs},
        m10_energy_fraction={k:g(k,'energy_fraction_m10') for k in legs},
        unperturbed_total_energy={k:g(k,'total_mode_energy') for k in legs if k.endswith('M0')},
        R48_over_R24_m10=(g('R48_M1','energy_m10')/g('R24_M1','energy_m10')
                          if g('R24_M1','energy_m10') else None))
    if g('R96_M1','energy_m10') is not None:
        ratio=g('R96_M1','energy_m10')/g('R48_M1','energy_m10')
        paired=('R96_M0',final) in idx
        disc.update(R96_over_R48_m10=ratio,
            reading=('noise_seeded' if ratio<0.1 else
                     'nonlinear_cascade' if ratio>0.1 and ratio<10 else 'neither_expectation'),
            threshold='a factor of 10 separates the two registered readings',
            r96_has_a_paired_unperturbed_leg=paired,
            caveat=('the R96 numerical floor is measured by its own unperturbed leg, so the '
                    'comparison no longer rests on assuming the floor falls with refinement'
                    if paired else
                    'R96 has no matching unperturbed leg, so its own numerical floor is unmeasured; '
                    'the R48 unperturbed total energy is used as a conservative bound only under the '
                    'assumption that the floor decreases with refinement, which the R24-R48 pair '
                    'supports but does not prove'))
        if paired:
            floor96=g('R96_M0','total_mode_energy')
            disc.update(R96_unperturbed_total_energy=floor96,
                R96_m10_over_its_own_floor=(g('R96_M1','energy_m10')/floor96) if floor96 else None,
                unperturbed_floor_ratio_R96_over_R48=(
                    floor96/g('R48_M0','total_mode_energy') if g('R48_M0','total_mode_energy') else None))
    # ---- is the short-wave energy a function of the resolved canonical mode? ---------
    # A converged nonlinear cascade from m=2 would make E(m=10) a function of the
    # resolved state. At matched canonical-mode energy the resolutions can therefore be
    # compared directly: a large spread means the short waves are not a cascade product.
    cascade=[]
    for label,fr in CHECKPOINTS:
        if fr>=nshared:continue
        row=dict(checkpoint=label,frame=fr,time_s=float(ref[fr]))
        for leg in [k for k in legs if k.endswith('M1')]:
            r=idx.get((leg,fr))
            if not r:continue
            row[f'{leg}_energy_m2']=r['energy_m2']
            row[f'{leg}_energy_m10']=r['energy_m10']
            row[f'{leg}_m10_over_m2']=(r['energy_m10']/r['energy_m2']) if r['energy_m2'] else None
        pairs2=[('R48_M1','R24_M1'),('R96_M1','R48_M1')]
        for hi,lo in pairs2:
            if f'{hi}_energy_m2' in row and f'{lo}_energy_m2' in row:
                row[f'{hi}_over_{lo}_m2']=row[f'{hi}_energy_m2']/row[f'{lo}_energy_m2'] if row[f'{lo}_energy_m2'] else None
                row[f'{hi}_over_{lo}_m10']=row[f'{hi}_energy_m10']/row[f'{lo}_energy_m10'] if row[f'{lo}_energy_m10'] else None
        cascade.append(row)
    c.write_csv(c.OUT/'cascade_test.csv',cascade,
                fields=sorted({k for r in cascade for k in r},key=lambda k:(k!='checkpoint',k)))

    summary=dict(task='command/16 Task C',shared_states=int(nshared),
        legs={k:dict(nx_cells=v['nx'],states=v['n'],stride=v['stride'],run=str(v['run'].relative_to(c.ROOT)))
              for k,v in legs.items()},
        time_axis=axis,checkpoints=[dict(label=l,frame=f) for l,f in CHECKPOINTS if f<nshared],
        uncertainty=dict(x_resolvability_floor=X_RESOLVABILITY,
            x_resolvability_source='cross-executable one-cell background bound over 0 -> T (wp9_corrections)',
            modal_energy_floor='the unperturbed leg at the same resolution and time',
            timestep_prior_on_growth=DT_GROWTH_PRIOR,
            timestep_prior_source='WP6 R24 dt/2 pair to 2T; it does not extend to 3T-4T, to R96 or to transport',
            transport_timestep_prior='unbounded: the same WP6 study found a 68 percent mid-window change at 1 mm'),
        discriminator=disc,
        cascade_test=dict(rows=cascade,
            question='is E(m=10) a function of the resolved canonical mode',
            reading=('a converged nonlinear cascade would make E(m=10) follow E(m=2); a large '
                     'spread at matched E(m=2) means the short waves are set by the discretization')),
        independence_note=('R48 and R96 are compared directly and never through R24; R24 is reported '
                           'alongside but no R96 result is required to reproduce it'))
    c.save(c.OUT/'resolution_comparison.json',summary)
    print(json.dumps(dict(legs=list(legs),shared_states=int(nshared),
                          time_axis_ok=all(v['within_budget'] for v in axis.values()),
                          discriminator=disc),indent=2,default=str))


if __name__=='__main__':main()
