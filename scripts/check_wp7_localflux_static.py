"""Static evidence that WP7 adds observation only and preserves integrated budgets."""
import hashlib,json,re
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
OLD=ROOT/'model/transport_wp6_v1/source/appl/simpleRichardsnc'
NEW=ROOT/'model/transport_wp7_localflux_v1/source/appl/simpleRichardsnc'
OUT=ROOT/'analysis/wp7_shortwave_validation_20261103'
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
checks={
    'main_byte_identical':sha(OLD/'main.cc')==sha(NEW/'main.cc'),
    'params_byte_identical':sha(OLD/'params.input')==sha(NEW/'params.input'),
    'integrated_budget_update_retained':all(token in (NEW/'transportdiagnostics.hh').read_text() for token in [
        'b[i].adv+=sign*flux.advective/width_',
        'b[i].diff+=sign*flux.diffusive/width_',
        'b[i].total+=sign*flux.total/width_',
        'b[i].waterOut+=sign*flux.water/width_']),
    'local_identity_columns':all(token in (NEW/'transportdiagnostics.hh').read_text() for token in [
        'sign*flux.water', 'sign*flux.advective', 'sign*flux.diffusive', 'sign*flux.total']),
    'local_output_opt_in':all(token in (NEW/'transportdiagnostics.hh').read_text() for token in [
        'Transport.WriteLocalFluxes', 'writeLocalFluxes_ && sign']),
    'upwind_value_is_observation_only':all(token in (NEW/'model.hh').read_text() for token in [
        'result[Indices::conti0EqIdx] = f.water;',
        'result[Indices::conti1EqIdx] = f.total;',
        'upVolVar.moleFraction()};']),
    'vtk_velocity_excluded_from_local_schema':'cell-centre Darcy reconstruction' in (NEW/'transportdiagnostics.hh').read_text(),
}
old_model=(OLD/'model.hh').read_text();new_model=(NEW/'model.hh').read_text()
old_return=re.search(r'return \{scvf\.area\(\) \* vn,.*?\};',old_model,re.S).group(0)
new_return=re.search(r'return \{scvf\.area\(\) \* vn,.*?\};',new_model,re.S).group(0)
checks['four_physical_flux_expressions_unchanged']=old_return.rstrip('};').strip()==new_return.rsplit(',',1)[0].strip()
result={'passed':all(checks.values()),'checks':checks,'old_main_sha256':sha(OLD/'main.cc'),'new_main_sha256':sha(NEW/'main.cc')}
(OUT/'local_flux_static_checks.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
assert result['passed']
