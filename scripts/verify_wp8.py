"""Hash, numerical and delivery audit, including honest native validation status."""
import hashlib,json,subprocess,sys,unittest,io
from pathlib import Path
from analyze_wp8 import ROOT,OUT,CACHE,RUNS,rows
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    checks={};result=json.loads((OUT/'results.json').read_text())
    for k,v in result['source_hashes'].items():checks['frozen/'+k]=sha(ROOT/k)==v
    for key,rel in RUNS.items():
        audit=json.loads((OUT/key/'audit.json').read_text());checks[key+'/fields']=audit['all_field_hashes_verified']
        checks[key+'/input']=audit['input_sha256']==sha(ROOT/rel/'params.input');checks[key+'/manifest']=audit['field_manifest_sha256']==sha(ROOT/rel/'post/field_hashes.json')
        checks[key+'/worker']=audit['script_sha256']==sha(ROOT/'scripts/wp8_fields.py')
        checks[key+'/cache']=sha(CACHE/key/'profiles_spectra.npz')==audit['files']['profiles_spectra.npz'];checks[key+'/csv']=sha(OUT/key/'timeseries.csv')==audit['files']['timeseries.csv']
    checks['full_shared_pair']=all(result['comparisons'][k]['shared_states']==4801 for k in ['R24_A','R48_A'])
    checks['background_window_only']=all(result['comparisons'][k]['end_s']<136423.0637723484 for k in ['R24_B','R24_C','R48_B','R48_C'])
    checks['failed_background_preserved']=not result['audit']['background_gate']['overall_preregistered_gate']
    checks['identity']=result['max_identity_error_m']<1e-18;checks['teff_unidentified']=result['t_eff'] is None
    from test_wp8 import TestWP8
    stream=io.StringIO();rr=unittest.TextTestRunner(stream=stream,verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(TestWP8));checks['numerical_tests']=rr.wasSuccessful()
    (OUT/'unit_tests.txt').write_text(stream.getvalue(),encoding='utf-8')
    reports=['METHOD','DATA_AUDIT','ZERO_MODE_RESPONSE','SPECTRAL_EVOLUTION','CONSERVATIVE_TRANSPORT','GRID_SENSITIVITY','BACKGROUND_VALIDATION_PLAN','SCIENTIFIC_DECISION']
    for n in reports:checks[n]=(OUT/(n+'_ZH.md')).exists()
    for c,name in zip('ABCDEFG',['mean_response','physical_wavenumbers','spectral_energy','conservative_transport','depth_profiles','early_background','own_mean_anomalies']):
        for ext in ['png','pdf']:checks['figure/'+c+'/'+ext]=(OUT/f'Figure_{c}_{name}.{ext}').exists()
    native={}
    for d in sorted(OUT.glob('background_smoke10*')):
        ex=json.loads((d/'execution.json').read_text());val=json.loads((d/'validation.json').read_text()) if (d/'validation.json').exists() else None
        native[d.name]=dict(returncode=ex.get('returncode'),validation=val,status='passed_short_bridge_only' if val and val['passed'] else 'failed_or_pending')
        receipt=json.loads((d/'transfer_receipt.json').read_text());cache=ROOT/'tmp/wp8_zero_mode_20261103'/d.name
        for name,v in receipt.items():
            p=(cache if name.endswith(('.vtu','.gz')) else d)/name;checks[d.name+'/'+name]=sha(p)==v['sha256']
    model=ROOT/'model/transport_wp8_background_v1';manifest=json.loads((model/'source_manifest.json').read_text())
    for name,v in manifest.items():checks['native_source/'+name]=sha(model/'source'/name)==v['sha256']
    build=json.loads((OUT/'build/r03_attempt_01/result.json').read_text());execution=json.loads((OUT/'background_smoke10_r03/execution.json').read_text())
    checks['native_build_exact_source']=build['source_sha256']=={k:v['sha256'] for k,v in manifest.items()}
    checks['native_run_exact_ELF']=execution['elf_sha256']==build['executable']['sha256']
    checks['native_short_bridge_passed']=native['background_smoke10_r03']['validation']['passed']
    previous=ROOT/'model/transport_wp7_localflux_v1/source'
    for p in previous.rglob('*'):
        if p.is_file() and p.name!='main.cc':checks['unchanged_native_kernel/'+p.relative_to(previous).as_posix()]=p.read_bytes()==(model/'source'/p.relative_to(previous)).read_bytes()
    audit=dict(passed=all(checks.values()),checks=checks,numerical_tests=rr.testsRun,native_validation=native,long_simulations_launched=False)
    (OUT/'verification.json').write_text(json.dumps(audit,indent=2)+'\n')
    files=[ROOT/'command/12.md']+list(OUT.rglob('*'))+list(model.rglob('*'))+[p for p in (ROOT/'scripts').iterdir() if p.is_file() and ('wp8' in p.name)]
    hashes={p.relative_to(ROOT).as_posix():dict(sha256=sha(p),bytes=p.stat().st_size) for p in files if p.is_file() and p.name!='provenance.json'}
    (OUT/'provenance.json').write_text(json.dumps(dict(reviewed_baseline='323fd639c81e47b201ecb0ae8b52ba8449892a59',files=hashes),indent=2)+'\n')
    print(json.dumps(dict(passed=audit['passed'],checks=len(checks),numerical_tests=rr.testsRun,failures=[k for k,v in checks.items() if not v],native_status={k:v['status'] for k,v in native.items()}),indent=2));assert audit['passed']
if __name__=='__main__':main()
