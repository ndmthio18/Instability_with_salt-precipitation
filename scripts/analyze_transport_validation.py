"""Audit the frozen small validation outputs. No solver or remote launch path."""
import csv
import hashlib
import json
import math
from pathlib import Path
import xml.etree.ElementTree as ET

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'analysis/transport_validation_20261102'
RAW = ROOT / 'data/generated/transport_validation_20261102/r01'
PLAN = json.loads((OUT / 'plan.json').read_text())
C = PLAN['criteria']
checks, summaries, regions = [], [], []


def check(name, condition, value=None, limit=None):
    checks.append(dict(name=name, passed=bool(condition), value=value, limit=limit))


def rows(path):
    with path.open(newline='') as f:
        return [{k: (v if k in ('region', 'frame', 'field') else float(v))
                 for k, v in r.items()} for r in csv.DictReader(f)]


def close(a, b, scale=1.0):
    return abs(a-b) <= 2e-14 * max(abs(a), abs(b), scale)


def field_audit(folder, budget, events, spec):
    """Reconstruct tensor-grid nodal BOX volumes from coordinates, not solver CVs."""
    phi = float(spec['parameters']['Problem']['Porosity'])
    datasets = ET.parse(next(folder.glob('*.pvd'))).getroot().iter('DataSet')
    inv_errors, surface_errors, spans = [], [], {'moleFraction': 0., 'sw': 0., 'pw': 0.}
    final_profiles = {}
    for ds in datasets:
        time = float(ds.attrib['timestep'])
        p = folder / ds.attrib['file']
        arrays = {a.attrib.get('Name'): list(map(float, (a.text or '').split()))
                  for a in ET.parse(p).getroot().iter('DataArray')
                  if a.attrib.get('format') == 'ascii'}
        coords = list(zip(*[iter(arrays['Coordinates'])]*3))
        xs, zs = sorted(set(v[0] for v in coords)), sorted(set(v[1] for v in coords))
        check(f'{spec["id"]}/{p.name}/tensor_grid', len(set(coords)) == len(xs)*len(zs))

        def nodal_widths(v):
            return {x: ((v[i+1] if i+1 < len(v) else x)
                        - (v[i-1] if i else x))/2 for i, x in enumerate(v)}

        wx, wz = nodal_widths(xs), nodal_widths(zs)
        width, height = xs[-1]-xs[0], zs[-1]-zs[0]
        water = [phi*s*wx[x]*wz[z]/width for (x, z, _), s in zip(coords, arrays['sw'])]
        salt = [v*x for v, x in zip(water, arrays['moleFraction'])]
        br = [r for r in budget if r['time_s'] == time]
        assert len(br) == 4
        for r in br:
            cut = zs[0]+height*r['actual_z_over_H']
            mask = [i for i, xyz in enumerate(coords) if xyz[1] > cut]
            for species, values in [('salt', salt), ('water', water)]:
                calc = math.fsum(values[i] for i in mask)
                inv_errors.append(abs(calc-r[f'inventory_{species}_m'])/abs(calc))
                calc_all = math.fsum(values)
                inv_errors.append(abs(calc_all-r[f'whole_{species}_inventory_m'])/abs(calc_all))
        ev = next(e for e in events if e['time_s'] == time)
        top = [i for i, xyz in enumerate(coords) if xyz[1] == zs[-1]]
        mean = math.fsum(arrays['moleFraction'][i]*wx[coords[i][0]] for i in top)/width
        surface_errors.append(abs(mean-ev['x_top_area_mean']))
        check(f'{spec["id"]}/{p.name}/extrema',
              close(max(arrays['moleFraction']), ev['x_domain_max'])
              and close(min(arrays['moleFraction']), ev['x_domain_min'])
              and close(max(arrays['moleFraction'][i] for i in top), ev['x_top_max']))
        for name in spans:
            profile = []
            for z in zs:
                values = [arrays[name][i] for i, xyz in enumerate(coords) if xyz[1] == z]
                spans[name] = max(spans[name], max(values)-min(values))
                profile.append((z, math.fsum(values)/len(values)))
            final_profiles[name] = profile
    for field, key in [('moleFraction', 'lateral_x_span_max'), ('sw', 'lateral_sw_span_max'), ('pw', 'lateral_pressure_span_Pa_max')]:
        check(spec['id']+'/uniform_'+field, spans[field] <= C[key], spans[field], C[key])
    check(spec['id']+'/independent_VTU_inventory', max(inv_errors) <= C['vtk_inventory_relative_max'], max(inv_errors), C['vtk_inventory_relative_max'])
    check(spec['id']+'/independent_surface_mean', max(surface_errors) <= 1e-14, max(surface_errors), 1e-14)
    return spans, max(inv_errors), final_profiles


all_budget, all_events, profiles = {}, {}, {}
remote_summary = json.loads((RAW/'summary.json').read_text())
check('remote_plan_unchanged', json.loads((RAW/'plan.json').read_text()) == PLAN)
for spec in PLAN['runs']:
    name = spec['id']
    folder = RAW/'runs'/name
    result = next(r for r in remote_summary if r['id'] == name)
    for filename, record in result['files'].items():
        check(name+'/sha256/'+filename, hashlib.sha256((folder/filename).read_bytes()).hexdigest() == record['sha256'])
    check(name+'/input_sha256', hashlib.sha256((folder/'params.input').read_bytes()).hexdigest() == spec['input_sha256'])
    stop = json.loads((folder/'stop_record.json').read_text())
    check(name+'/exit', result['returncode'] == spec['expected_exit'])
    check(name+'/stop', stop['reason'] == spec['expected_reason'])
    summary = dict(id=name, reason=stop['reason'], accepted_steps=stop['accepted_steps'],
                   final_time_s=stop['last_accepted_time_s'], wall_s=result['wall_s'])
    summaries.append(summary)
    if name == 'startup_invalid_time':
        check(name+'/no_false_output', not (folder/'transport_budget.csv').exists()
              and not stop['initial_diagnostics_written'] and stop['last_accepted_time_s'] is None)
        continue
    budget, events = rows(folder/'transport_budget.csv'), rows(folder/'surface_events.csv')
    all_budget[name], all_events[name] = budget, events
    check(name+'/row_counts', len(budget) == 4*(stop['accepted_steps']+1) and len(events) == stop['accepted_steps']+1)
    check(name+'/finite', all(math.isfinite(v) for row in budget+events for v in row.values() if isinstance(v, float)))
    check(name+'/event_time_order', all(e['step'] == i and (i == 0 or close(e['time_s']-events[i-1]['time_s'], e['dt_s'])) for i, e in enumerate(events)))
    check(name+'/final_time', events[-1]['time_s'] == stop['last_accepted_time_s'])
    check(name+'/surface_uniform_every_step', max(e['sigma_top_node'] for e in events) <= C['lateral_x_span_max'])
    check(name+'/prescribed_evaporation', all(close(r['other_water_out_m_s'], float(spec['parameters']['Problem']['EvaporationRate']), 1e-10) for r in budget))
    for region in sorted(set(r['region'] for r in budget)):
        rr = [r for r in budget if r['region'] == region]
        entry = dict(id=name, region=region, actual_z_over_H=rr[0]['actual_z_over_H'])
        for species in ('salt', 'water'):
            initial = rr[0][f'inventory_{species}_m']
            previous, cumulative, abs_flux = initial, 0., 0.
            max_rel, max_abs, max_flux_rel, max_step, algebra_error = 0., 0., 0., 0., 0.
            flux_gate = True
            floor_rescued_rows = 0
            adv_cum, diff_cum = 0., 0.
            for i, r in enumerate(rr):
                dt = r['dt_s']
                flux = r['cut_salt_total_down_m_s']+r['other_salt_out_m_s'] if species == 'salt' else r['cut_water_down_m_s']+r['other_water_out_m_s']
                cumulative += dt*flux
                abs_flux += dt*abs(flux)
                inventory = r[f'inventory_{species}_m']
                residual = inventory-initial+cumulative
                recorded = r[f'{species}_residual_m']
                step_residual = inventory-previous+dt*flux
                previous = inventory
                col = 'cumulative_cut_salt_m' if species == 'salt' else 'cumulative_net_water_out_m'
                algebra_error = max(algebra_error, abs(cumulative-r[col]), abs(residual-recorded), abs(step_residual-r[f'step_{species}_residual_m']))
                max_rel = max(max_rel, abs(recorded)/initial)
                max_abs = max(max_abs, abs(recorded))
                max_step = max(max_step, abs(step_residual))
                if abs_flux > C['absolute_residual_floor_m']:
                    max_flux_rel = max(max_flux_rel, abs(recorded)/abs_flux)
                flux_gate &= abs(recorded) <= max(C['absolute_residual_floor_m'], C['relative_residual_to_cumulative_absolute_flux_max']*abs_flux)
                if abs(recorded) > C['relative_residual_to_cumulative_absolute_flux_max']*abs_flux and abs(recorded) <= C['absolute_residual_floor_m']:
                    floor_rescued_rows += 1
                if species == 'salt':
                    adv_cum += dt*r['cut_salt_adv_down_m_s']
                    diff_cum += dt*r['cut_salt_diff_down_m_s']
                    algebra_error = max(algebra_error, abs(adv_cum-r['cumulative_cut_adv_m']), abs(diff_cum-r['cumulative_cut_diff_m']), abs(r['cut_salt_total_down_m_s']-r['cut_salt_adv_down_m_s']-r['cut_salt_diff_down_m_s']))
            entry.update({species+'_residual_abs_max_m': max_abs, species+'_residual_rel_initial_max': max_rel,
                          species+'_residual_rel_abs_flux_max': max_flux_rel, species+'_step_residual_abs_max_m': max_step,
                          species+'_flux_gate_floor_rescued_rows': floor_rescued_rows,
                          species+'_final_inventory_change_m': rr[-1][f'inventory_{species}_m']-initial,
                          species+'_final_net_out_m': cumulative})
            check(f'{name}/{region}/{species}_inventory_balance', max_rel <= C['relative_cumulative_inventory_residual_max'], max_rel, C['relative_cumulative_inventory_residual_max'])
            check(f'{name}/{region}/{species}_flux_balance', flux_gate,
                  dict(max_relative_when_flux_above_floor=max_flux_rel, floor_rescued_rows=floor_rescued_rows),
                  dict(relative=C['relative_residual_to_cumulative_absolute_flux_max'], absolute_floor_m=C['absolute_residual_floor_m']))
            check(f'{name}/{region}/{species}_CSV_algebra', algebra_error <= 1e-15, algebra_error, 1e-15)
        regions.append(entry)
    # Whole-domain reaction is an inferred quantity: check its accounting, not independence.
    free = [r for r in budget if r['region'] == 'free_domain']
    for species in ('salt', 'water'):
        cumulative = 0.
        err = 0.
        for r in free:
            top = r[f'other_{species}_out_m_s']
            cumulative += r['dt_s']*(top-r[f'bottom_{species}_reaction_in_m_s'])
            residual = r[f'whole_{species}_inventory_m']-free[0][f'whole_{species}_inventory_m']+cumulative
            err = max(err, abs(residual-r[f'whole_{species}_residual_m']))
        check(f'{name}/whole_{species}_reaction_accounting', err <= 1e-15, err, 1e-15)
    spans, inv_error, profiles[name] = field_audit(folder, budget, events, spec)
    summary.update(lateral_x_span_max=spans['moleFraction'], lateral_sw_span_max=spans['sw'],
                   lateral_pressure_span_Pa_max=spans['pw'], independent_inventory_rel_error_max=inv_error,
                   surface_sigma_all_steps_max=max(e['sigma_top_node'] for e in events),
                   final_top_x=events[-1]['x_top_max'], threshold_low_s=stop['surface_threshold_low_s'], threshold_high_s=stop['surface_threshold_high_s'])
    for species in ('salt', 'water'):
        summary[species+'_balance_rel_max'] = max(r[species+'_residual_rel_initial_max'] for r in regions if r['id'] == name)
    if stop['reason'] == 'requested_horizon':
        check(name+'/no_false_threshold', all(e['surface_threshold_reached'] == e['global_threshold_reached'] == 0 and e['x_domain_max'] < .1 for e in events)
              and stop['surface_threshold_low_s'] is None and stop['global_threshold_low_s'] is None)
    else:
        last = events[-1]
        check(name+'/threshold_last_step', last['x_top_max'] >= .1 and last['surface_threshold_reached'] == 1
              and all(e['x_top_max'] < .1 and e['surface_threshold_reached'] == 0 for e in events[:-1]))
        check(name+'/threshold_bracket', stop['surface_threshold_high_s'] == last['time_s']
              and stop['surface_threshold_low_s'] == (events[-2]['time_s'] if len(events)>1 else 0)
              and close(stop['surface_threshold_high_s']-stop['surface_threshold_low_s'], last['dt_s']))
        check(name+'/global_bracket', stop['global_threshold_low_s'] == stop['surface_threshold_low_s']
              and stop['global_threshold_high_s'] == stop['surface_threshold_high_s'])

hydro = all_budget['closed_hydrostatic']
max_flux = max(abs(r[k]) for r in hydro for k in ['cut_water_down_m_s', 'cut_salt_total_down_m_s'])
check('closed_hydrostatic/flux', max_flux <= C['hydrostatic_flux_abs_m_s_max'], max_flux, C['hydrostatic_flux_abs_m_s_max'])
p8, p16 = profiles['evap_nx8_dt20']['moleFraction'], profiles['evap_nx16_dt20']['moleFraction']
assert [p[0] for p in p8] == [p[0] for p in p16]
profile_diff = max(abs(a[1]-b[1]) for a, b in zip(p8, p16))/max(abs(a[1]) for a in p16)
check('horizontal_refinement/x_profile', profile_diff <= C['horizontal_refinement_profile_x_relative_max'], profile_diff, C['horizontal_refinement_profile_x_relative_max'])
inv_diff = max(abs(a[f'inventory_{s}_m']-b[f'inventory_{s}_m'])/b[f'inventory_{s}_m'] for a, b in zip(all_budget['evap_nx8_dt20'][-4:], all_budget['evap_nx16_dt20'][-4:]) for s in ('salt', 'water'))
check('horizontal_refinement/inventory', inv_diff <= C['horizontal_refinement_inventory_relative_max'], inv_diff, C['horizontal_refinement_inventory_relative_max'])
dt_diff = abs(all_events['evap_nx8_dt20'][-1]['x_top_max']-all_events['evap_nx8_dt10'][-1]['x_top_max'])/all_events['evap_nx8_dt10'][-1]['x_top_max']
check('dt_halving/top_x', dt_diff <= C['dt_halving_top_x_relative_max'], dt_diff, C['dt_halving_top_x_relative_max'])
comparisons = dict(hydrostatic_flux_abs_max_m_s=max_flux, horizontal_profile_x_relative=profile_diff,
                   horizontal_inventory_relative=inv_diff, dt_halving_top_x_relative=dt_diff)
for filename, data in [('run_summary.csv', summaries), ('budget_summary.csv', regions)]:
    fields = list(dict.fromkeys(k for row in data for k in row))
    with (OUT/filename).open('w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader(); w.writerows(data)
report = dict(passed=all(c['passed'] for c in checks), check_count=len(checks), failed=[c for c in checks if not c['passed']],
              comparisons=comparisons, runs=summaries, checks=checks,
              limitations=['uniform 2D only; not independently coded 1D', 'no vertical mesh convergence',
                           'bottom reaction inferred, not an independent physical boundary flux',
                           'event stress tests use synthetic high initial salinity', 'no Pilot A/B or effective transport label'])
(OUT/'validation_checks.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
print(json.dumps({k: report[k] for k in ['passed', 'check_count', 'failed', 'comparisons', 'runs']}, indent=2))
raise SystemExit(0 if report['passed'] else 1)
