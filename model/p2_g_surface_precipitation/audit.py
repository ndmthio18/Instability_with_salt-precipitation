"""Audit complete accepted histories on a Slurm compute node (no new solves)."""
import csv
import hashlib
import json
import math
import pathlib

root = pathlib.Path.cwd()
out = root / 'results'
histories = {}
checks = {}
for path in sorted((root / 'runs').iterdir()):
    with (path / 'surface_precipitation.csv').open() as f:
        rows = [{k: float(v) for k, v in r.items()} for r in csv.DictReader(f)]
    histories[path.name] = rows
    first, last = rows[0], rows[-1]
    projected = [r for r in rows if r['precipitated_step_kg_m2'] > 1e-15]
    onset_index = rows.index(projected[0]) if projected else None
    salt_projection = max(abs(r['projection_salt_residual_kg']) for r in rows)
    water_projection = max(abs(r['projection_water_residual_kg']) for r in rows)
    accumulation = abs(math.fsum(r['precipitated_step_kg_m2'] for r in rows)-last['solid_kg_m2'])
    growth_error = 0.0
    for prev, r in zip(rows, rows[1:]):
        assert r['time_s'] > prev['time_s']
        assert abs(r['time_s']-prev['time_s']-r['dt_s']) < 1e-12
        assert r['solid_kg_m2'] >= prev['solid_kg_m2']
        if not prev['prat_active']:
            assert r['prat_demand_kg_m2'] == 0
            continue
        # Independent evaluation of the user-specified sqrt(Da*) form.
        eps, ceq = .43, .1*.05844/(.9*.01801528+.1*.05844)
        diffusion, kr, av = math.sqrt(eps)*1.3e-9, .0023, 36*(1-eps)/25e-6
        h = prev['h_prat_m']
        da = h*h*(1-ceq)*av*kr/(eps*diffusion)
        da_star = 4*(1-eps)**2*(1-ceq)*kr/(eps*diffusion*av)
        z = math.sqrt(da_star)*math.tanh(math.sqrt(da))
        predicted = ceq*r['evaporation_kg_m2_s']/(1-ceq)*z/(1+z)*r['dt_s']
        growth_error = max(growth_error, abs(predicted-r['prat_demand_kg_m2'])/predicted)
    checks[path.name] = {
        'accepted_steps': len(rows)-1,
        'onset_bracket_s': [rows[onset_index-1]['time_s'], rows[onset_index]['time_s']] if onset_index else None,
        'max_projection_salt_residual_kg': salt_projection,
        'max_projection_water_residual_kg': water_projection,
        'solid_accumulation_error_kg_m2': accumulation,
        'growth_formula_relative_error': growth_error,
        'finite_to_zero_reaction_demand_ratio': (
            math.fsum(r['prat_demand_kg_m2'] for r in rows)/
            math.fsum(r['prat_zero_reaction_demand_kg_m2'] for r in rows)
            if any(r['prat_active'] for r in rows) else None),
        'minimum_surface_sw': min(r['min_surface_sw'] for r in rows),
    }
    assert last['time_s'] == 1000
    assert salt_projection/first['soil_salt_kg'] < 1e-12
    assert water_projection/first['soil_water_kg'] < 1e-12
    assert accumulation < 1e-14 and growth_error < 1e-12

physics = ('time_s', 'dt_s', 'evaporation_kg_m2_s', 'precipitated_step_kg_m2',
           'solid_kg_m2', 'h_mass_m', 'top_mass_fraction', 'min_surface_sw',
           'soil_salt_kg', 'soil_water_kg', 'evaporated_water_kg')
for name in ('threshold_half', 'threshold_double'):
    assert len(histories[name]) == len(histories['n80'])
    assert all(a[k] == b[k] for a, b in zip(histories[name], histories['n80']) for k in physics)

control = histories['no_precip_control']
nominal = histories['n80']
prefix = [i for i, r in enumerate(nominal) if r['solid_kg_m2'] == 0]
assert all(nominal[i][k] == control[i][k] for i in prefix for k in physics)
summaries = {s['case']: s for s in json.loads((out / 'summary.json').read_text())}
comparisons = {}
for a, b, gate in [('n40', 'n80', .02), ('n80', 'n160', .02), ('n80', 'dt_half', .01)]:
    metrics = {k: abs(summaries[a][k]-summaries[b][k])/abs(summaries[b][k])
               for k in ('solid_kg_m2', 'h_prat_m', 'demand_over_actual', 'onset_s', 'continuum_s')}
    comparisons[a+'__'+b] = dict(relative_to=b, tolerance=gate,
                                  relative_differences=metrics,
                                  final_solid_pass=metrics['solid_kg_m2'] < gate)

result = dict(
    checks=checks, comparisons=comparisons,
    threshold_changes_leave_soil_histories_identical=True,
    no_precipitation_prefix_identical=True,
    no_precipitation_prefix_end_s=nominal[prefix[-1]]['time_s'],
    final_control_top_mass_fraction=control[-1]['top_mass_fraction'],
    final_precipitating_top_mass_fraction=nominal[-1]['top_mass_fraction'],
    final_evaporation_relative_change=nominal[-1]['evaporation_kg_m2_s']/control[-1]['evaporation_kg_m2_s']-1,
    compatibility_within_10_percent=abs(summaries['n80']['demand_over_actual']-1) < .1,
    note='Endpoint dt gate passes; grid and onset prediction not established. No 2D inference.')
(out / 'audit.json').write_text(json.dumps(result, indent=2))

# Per-run replay metadata and full accepted histories are small, unlike raw VTK.
bundle = out / 'replay'
for path in sorted((root / 'runs').iterdir()):
    dest = bundle / path.name
    dest.mkdir(parents=True, exist_ok=True)
    for name in ('params.input', 'command.json', 'exit_status.txt', 'accepted_steps.csv', 'surface_precipitation.csv'):
        (dest / name).write_bytes((path / name).read_bytes())
manifest = {str(p.relative_to(root)): hashlib.sha256(p.read_bytes()).hexdigest()
            for base in (root / 'runs', root / 'source') for p in sorted(base.rglob('*')) if p.is_file()}
manifest['build/solver'] = hashlib.sha256((root / 'build/solver').read_bytes()).hexdigest()
manifest['source.tar.gz'] = hashlib.sha256((root / 'source.tar.gz').read_bytes()).hexdigest()
(out / 'remote_manifest.json').write_text(json.dumps(manifest, indent=2))
print(json.dumps(result, indent=2))
