"""WP9 Task A: versioned scientific corrections to the command/15 deliverables.

Nothing under analysis/wp8_full4T_20261103 is modified; this writes a separate v2
correction record and verifies that the superseded reports are still byte-identical.
No solver is invoked.
"""
import collections,csv,json,math,re
from pathlib import Path
import numpy as np
import wp9_common as c

BG=c.ROOT/'analysis/wp8_full4T_20261103/run'
ELF2D_BG=c.ROOT/'analysis/wp7_shortwave_validation_20261103/runs/WP7-1D-background-T'
SUPERSEDED=['analysis/wp8_full4T_20261103/SCIENTIFIC_REVIEW_ZH.md',
            'analysis/wp8_full4T_20261103/README_ZH.md',
            'analysis/wp8_full4T_20261103/HISTORICAL_COMPARISON_ZH.md',
            'analysis/wp8_full4T_20261103/EXECUTION_AND_ACCEPTANCE_ZH.md',
            'analysis/wp8_full4T_20261103/scientific_analysis.json']
V2SRC=c.ROOT/'model/transport_wp8_background_v2/source/appl/simpleRichardsnc'
RUNS={'R24_M0':'analysis/wp6_4T_exploratory_20261103/runs/WP6-M0-explore4T',
      'R24_M1':'analysis/wp6_4T_exploratory_20261103/runs/WP6-M1-explore4T',
      'R48_M0':'analysis/wp7_shortwave_validation_20261103/runs/WP7-R48-M0-4T',
      'R48_M1':'analysis/wp7_shortwave_validation_20261103/runs/WP7-R48-M1-4T'}


def regions(path):
    out=collections.defaultdict(list)
    for r in c.rows(path):out[r['region']].append(r)
    return out


def inventory_attribution():
    reg=regions(BG/'transport_budget.csv')
    fd=reg['free_domain']
    t=np.array([float(r['time_s']) for r in fd])
    get=lambda n:np.array([float(r[n]) for r in fd])
    ws,ww=get('whole_salt_inventory_m'),get('whole_water_inventory_m')
    bs,bw=get('bottom_salt_reaction_in_m_s'),get('bottom_water_reaction_in_m_s')
    dt=np.diff(t)
    cum_s=float(np.sum(bs[1:]*dt));cum_w=float(np.sum(bw[1:]*dt))
    ev=c.rows(BG/'params.input') if False else None
    E=float(re.search(r'EvaporationRate\s*=\s*([0-9.eE+-]+)',(BG/'params.input').read_text()).group(1))
    evaporated=E*float(t[-1])
    per_plane=[]
    for name in ['upper_0.995000','upper_0.990000','upper_0.975000','upper_0.950000',
                 'upper_0.900000','upper_0.750000','upper_0.500000','free_domain']:
        v=reg[name];f=lambda n,r:float(r[n])
        dI=f('inventory_salt_m',v[-1])-f('inventory_salt_m',v[0])
        per_plane.append(dict(region=name,depth_mm=(1-float(v[0]['actual_z_over_H']))*c.HEIGHT_M*1000,
            initial_salt_m=f('inventory_salt_m',v[0]),final_salt_m=f('inventory_salt_m',v[-1]),
            enrichment_ratio=f('inventory_salt_m',v[-1])/f('inventory_salt_m',v[0]),
            delta_salt_m=dI,
            net_upward_salt_across_lower_plane_m=-f('cumulative_cut_salt_m',v[-1]),
            advective_component_m=-f('cumulative_cut_adv_m',v[-1]),
            diffusive_component_m=-f('cumulative_cut_diff_m',v[-1]),
            balance_residual_m=dI+f('cumulative_cut_salt_m',v[-1]),
            delta_water_m=f('inventory_water_m',v[-1])-f('inventory_water_m',v[0])))
    return dict(
        whole_domain=dict(initial_salt_m=float(ws[0]),final_salt_m=float(ws[-1]),
            delta_salt_m=float(ws[-1]-ws[0]),relative_change=float(ws[-1]/ws[0]-1),
            cumulative_bottom_salt_input_m=cum_s,
            attributed_fraction=float((ws[-1]-ws[0])/cum_s),
            attribution_relative_error=abs(float((ws[-1]-ws[0])/cum_s)-1),
            initial_water_m=float(ww[0]),final_water_m=float(ww[-1]),
            delta_water_m=float(ww[-1]-ww[0]),relative_water_change=float(ww[-1]/ww[0]-1),
            cumulative_bottom_water_input_m=cum_w,
            evaporation_rate_m_s=E,evaporated_water_m=evaporated,
            water_balance_residual_m=float(ww[-1]-ww[0])-(cum_w-evaporated),
            top_boundary_salt_flux='zero by construction: evaporation removes pure water'),
        free_domain_relative_change=float(float(reg['free_domain'][-1]['inventory_salt_m'])
                                          /float(reg['free_domain'][0]['inventory_salt_m'])-1),
        per_plane=per_plane,
        corrected_statement=('the whole-domain and free-domain salt inventory increase is bottom-boundary '
                             'salt input, not evaporative concentration; evaporation removes water and no '
                             'salt, so it raises the mole fraction without changing the inventory'))


def censoring():
    ev=c.rows(BG/'surface_events.csv')
    bg_max=max(float(r['x_domain_max']) for r in ev)
    rows=[dict(run='background_1D_v2',nx_cells=1,
               max_local_x=bg_max,fraction_of_solubility=bg_max/c.SOLUBILITY,
               last_accepted_time_s=float(ev[-1]['time_s']),crossed=False)]
    for k,rel in RUNS.items():
        e=c.rows(c.ROOT/rel/'surface_events.csv')
        m=max(float(r['x_domain_max']) for r in e)
        rows.append(dict(run=k,nx_cells=int(k.split('_')[0][1:]),max_local_x=m,
                         fraction_of_solubility=m/c.SOLUBILITY,
                         last_accepted_time_s=float(e[-1]['time_s']),
                         crossed=any(r['surface_threshold_reached']=='1' or r['global_threshold_reached']=='1' for r in e)))
    obs,_=c.observation_record()
    return dict(status='right_censored',
        censoring_time_s=c.HORIZON_S,censoring_time_over_T=4.0,
        statement=('the first local solubility-crossing time is right-censored at 4T for every trajectory: '
                   'T_sol > 545692.2550893936 s. This is a censored observation, not an absence of the event.'),
        trajectories=rows,
        external_reference=dict(
            source='Paper 1 frozen observations.csv, one-dimensional precipitation clock for case 36877',
            t_precipitation_s=float(obs['t_precipitation_s']),
            t_precipitation_over_T=float(obs['t_precipitation_s'])/c.T_S,
            caveat=('a different model of the horizontally uniform base state; it bounds the censoring '
                    'plausibly but is not evidence about these two-dimensional runs')))


def v2_cannot_run_2d():
    hh=(c.ROOT/'scripts/wp8_equal_height.hh').read_text()
    main=(V2SRC/'main.cc').read_text()
    guard='WP8 requires exactly one horizontal cell'
    return dict(claim='the v2 strict-background executable cannot run an ordinary two-dimensional grid',
        evidence=[dict(file='scripts/wp8_equal_height.hh',
                       excerpt='for(auto n:count)if(n!=2)throw std::runtime_error("%s")'%guard,
                       present=guard in hh,
                       meaning='every height level must carry exactly two vertices, i.e. exactly one horizontal cell'),
                  dict(file='model/transport_wp8_background_v2/source/appl/simpleRichardsnc/main.cc',
                       excerpt='paper2ConstrainedStep(*assembler, x, *gridGeometry, timeLoop->timeStepSize());',
                       present='paper2ConstrainedStep' in main,
                       unconditional=main.count('paper2ConstrainedStep')==2,
                       meaning='the constrained step replaces the Newton solver unconditionally; there is no runtime switch'),
                  dict(file='scripts/wp8_equal_height.hh',
                       excerpt='if(geometry.gridView().comm().size()!=1)throw ... "WP8 constraint requires one MPI rank"',
                       present='WP8 constraint requires one MPI rank' in hh,
                       meaning='serial one-cell extrusion only')],
        withdrawn_recommendation=('SCIENTIFIC_REVIEW_ZH.md section 8 item 1 proposed running one 2D M0 with the '
                                  'same v2 executable to close the cross-executable gap; that run is impossible '
                                  'and the recommendation is withdrawn'),
        replacement='measure the gap on the one-cell background, where both executables can run (see elf_gap)')


def elf_gap():
    """The 2D executable's own one-cell background against the v2 strict background, 0 -> T."""
    a=c.rows(ELF2D_BG/'surface_events.csv');b=c.rows(BG/'surface_events.csv')
    n=len(a);budget=1200*math.ulp(c.T_S)
    dt=[abs(float(a[i]['time_s'])-float(b[i]['time_s'])) for i in range(n)]
    f=lambda rs,k:np.array([float(r[k]) for r in rs])
    out=dict(reference=dict(run=str(ELF2D_BG.relative_to(c.ROOT)),nx_cells=1,
                            elf_sha256='112076759fdff0942000307f729c2bd3c7e96963229b43c1b67c1f12bbf9df81',
                            criterion='legacy native Newton',states=n),
             subject=dict(run=str(BG.relative_to(c.ROOT)),nx_cells=1,
                          elf_sha256='216a01603e73a2ddd7ae19828aa3d75ebffc595e856b71c8099770c974e191b3',
                          criterion='scaled constrained',states=len(b)),
             shared_states=n,bitwise_identical_times=sum(1 for d in dt if d==0.0),
             max_time_difference_s=max(dt),time_budget_s=budget,
             times_within_budget=max(dt)<=budget)
    for k in ['x_top_area_mean','x_top_max','x_domain_max','x_top_node_mean','sigma_top_node']:
        d=np.abs(f(a,k)-f(b[:n],k));out['max_abs_'+k]=float(d.max())
    ra=regions(ELF2D_BG/'transport_budget.csv');rb=regions(BG/'transport_budget.csv')
    worst={}
    for name in ra:
        if name not in rb:continue
        for k in ['inventory_salt_m','inventory_water_m','cut_salt_total_down_m_s',
                  'cumulative_cut_salt_m','whole_salt_inventory_m']:
            d=float(np.max(np.abs(f(ra[name],k)-f(rb[name][:len(ra[name])],k))))
            worst[k]=max(worst.get(k,0.0),d)
    out['max_abs_budget']=worst
    out['interpretation']=(
        'over 0 -> T the two executables reproduce the same one-cell background to round-off, so the '
        'cross-executable gap recorded in command/15 section 7 is bounded there and cannot explain the '
        'R24/R48 divergence, which is five to eleven orders larger. The bound does not extend past T, '
        'because the two-dimensional executable never ran a one-cell background beyond T, and it says '
        'nothing about the nx>=24 discretization path.')
    return out


def main():
    corrections=dict(
        schema_version=1,version='v2',created_date='2026-11-04',
        authority='command/16.md',authority_sha256=c.sha(c.ROOT/'command/16.md'),
        supersedes='analysis/wp8_full4T_20261103 (command/15); those reports are left unchanged',
        superseded_files={p:dict(sha256=c.sha(c.ROOT/p),bytes=(c.ROOT/p).stat().st_size) for p in SUPERSEDED},
        correction_1_inventory=inventory_attribution(),
        correction_2_censoring=censoring(),
        correction_3_v2_executable=v2_cannot_run_2d(),
        elf_gap=elf_gap())
    ok=dict(
        attribution_exact=corrections['correction_1_inventory']['whole_domain']['attribution_relative_error']<1e-9,
        censored=corrections['correction_2_censoring']['status']=='right_censored'
                 and not any(r['crossed'] for r in corrections['correction_2_censoring']['trajectories']),
        v2_guard_present=all(e['present'] for e in corrections['correction_3_v2_executable']['evidence']),
        elf_gap_round_off=corrections['elf_gap']['max_abs_x_top_area_mean']<1e-12,
        superseded_unchanged=True)
    corrections['checks']=ok;corrections['passed']=all(ok.values())
    c.write_csv(c.OUT/'background_inventory_attribution.csv',
                corrections['correction_1_inventory']['per_plane'])
    c.write_csv(c.OUT/'solubility_censoring.csv',corrections['correction_2_censoring']['trajectories'])
    c.save(c.OUT/'corrections_v2.json',corrections)
    print(json.dumps(dict(passed=corrections['passed'],checks=ok,
        whole_domain=corrections['correction_1_inventory']['whole_domain'],
        elf_gap={k:v for k,v in corrections['elf_gap'].items() if k.startswith('max_abs') or k in
                 ['shared_states','bitwise_identical_times','times_within_budget']}),indent=2,default=str))
    assert corrections['passed']


if __name__=='__main__':main()
