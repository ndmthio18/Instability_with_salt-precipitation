"""Numerical acceptance and provenance checks for the av sensitivity matrix."""
import csv,json,math,pathlib

root=pathlib.Path.cwd(); analysis=root/'analysis'; results=analysis/'results'
def read(path):
    def value(key, raw):
        if key in ('case', 'reference', 'stage', 'status'):
            return raw
        try:
            return float(raw)
        except ValueError:
            return raw
    with path.open(newline='') as handle:
        return [{key:value(key, raw) for key,raw in row.items()} for row in csv.DictReader(handle)]
summary=read(results/'summary.csv')
primary=[row for row in summary if row['stage']=='primary_exact']
halves=[row for row in summary if row['stage']=='dt_half_exact']
refined=[row for row in summary if row['stage'] in ('ale_refine_exact','ale256_exact')]
identity=read(results/'preactivation_identity.csv')
activation_identity=read(results/'activation_identity.csv')
by={(int(row['soil_cells']),row['av_multiplier']):row for row in primary}
by_half={(int(row['soil_cells']),row['av_multiplier']):row for row in halves}
checks={
    'primary_case_count':len(primary),
    'primary_all_reached_2500_s':len(primary)==6 and all(row['status']=='COMPLETED' and row['actual_end_s']>=2500 for row in primary),
    'primary_all_cases_accounted':len(primary)==6,
    'stopped_cases_have_rejected_bound_state':all(row['status']=='COMPLETED' or math.isfinite(row['first_rejected_time_s']) for row in primary),
    'preactivation_exact_identity':len(identity)==4 and all(row['maximum_absolute_difference']==0 and row['row_count_reference']==row['row_count_case'] for row in identity),
    'activation_exact_identity':len(activation_identity)==4 and all(
        row['activation_difference_s']==0 and row['initial_ale_max_abs']==0 and row['post_transfer_soil_max_abs']==0 and
        row['donor_profile_max_abs']==0 and row['ale_rows_reference']==row['ale_rows_case'] and
        row['soil_rows_reference']==row['soil_rows_case'] and row['donor_rows_reference']==row['donor_rows_case']
        for row in activation_identity),
    'zero_initial_solid':json.loads((analysis/'verification.json').read_text())['all_zero_initial_solid'],
    'global_mass_gate':max([row['salt_relative_residual'] for row in summary]+[row['water_relative_residual'] for row in summary])<=1e-8,
    'historical_m3_regression':(root/'stage_a_exact/verification.log').exists() and 'PASS' in (root/'stage_a_exact/verification.log').read_text(),
    'av_unit_test':(root/'unit_exact/av_parameter_test_r02.log').exists() and 'PASS' in (root/'unit_exact/av_parameter_test_r02.log').read_text(),
}
dt=[]
for key,row in by.items():
    if key not in by_half: continue
    half=by_half[key]
    differences={field:abs(row[field]-half[field])/max(abs(half[field]),1e-30) for field in
                 ('final_height_m','net_internal_solid_kg_m2','cumulative_positive_internal_kg_m2','cumulative_new_top_solid_kg_m2','epsilon_min')}
    dt.append(dict(soil_cells=key[0],av_multiplier=key[1],**differences,maximum_relative_difference=max(differences.values())))
checks['timestep_cases']=len(dt)
checks['timestep_maximum_relative_difference']=max((row['maximum_relative_difference'] for row in dt),default=math.nan)
checks['timestep_one_percent_gate']=len(dt)>=2 and checks['timestep_maximum_relative_difference']<=.01
with (results/'timestep_comparison.csv').open('w',newline='') as handle:
    if dt:
        writer=csv.DictWriter(handle,dt[0].keys()); writer.writeheader(); writer.writerows(dt)

ale=[]
coarse={(int(row['soil_cells']),row['av_multiplier']):row for row in primary}
for row in refined:
    key=(int(row['soil_cells']),row['av_multiplier'])
    if key not in coarse: continue
    base=coarse[key]
    ale.append(dict(soil_cells=key[0],av_multiplier=key[1],coarse_cells=int(base['ale_cells']),refined_cells=int(row['ale_cells']),
        height_relative=abs(base['final_height_m']-row['final_height_m'])/max(abs(row['final_height_m']),1e-30),
        internal_relative=abs(base['net_internal_solid_kg_m2']-row['net_internal_solid_kg_m2'])/max(abs(row['net_internal_solid_kg_m2']),1e-30),
        l95_shortest_relative=abs(base['l95_shortest_m']-row['l95_shortest_m'])/max(abs(row['l95_shortest_m']),1e-30),
        epsilon_absolute=abs(base['epsilon_min']-row['epsilon_min'])))
if ale:
    with (results/'ale_resolution_comparison.csv').open('w',newline='') as handle:
        writer=csv.DictWriter(handle,ale[0].keys()); writer.writeheader(); writer.writerows(ale)
checks['ale_refinement_cases']=len(ale)
checks['ale_resolution_gate']=bool(ale) and max(max(row['height_relative'],row['internal_relative'],row['l95_shortest_relative']) for row in ale)<=.02

observations={}
for n in (80,160):
    rows=[by[(n,av)] for av in (.1,1.,10.) if (n,av) in by]
    if len(rows)==3:
        observations[str(n)]={field:[row[field] for row in rows] for field in
            ('cumulative_positive_internal_kg_m2','net_internal_solid_kg_m2','cumulative_new_top_solid_kg_m2',
             'internal_fraction_net','l95_bottom_m','l95_shortest_m','epsilon_min','final_height_m')}
output=dict(checks=checks,numerical_observations_by_increasing_av=observations,
            note='Monotonicity is reported as an observation, not an acceptance condition.')
(analysis/'verification_final.json').write_text(json.dumps(output,indent=2,allow_nan=True))
print(json.dumps(output,indent=2,allow_nan=True))
if not all(checks[key] for key in ('primary_all_cases_accounted','stopped_cases_have_rejected_bound_state','preactivation_exact_identity','activation_exact_identity','zero_initial_solid','global_mass_gate','historical_m3_regression','av_unit_test')):
    raise SystemExit(1)
