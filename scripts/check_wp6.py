"""Static/provenance checks and two bounded threshold-fixture checks."""
import argparse,csv,hashlib,json,subprocess
from pathlib import Path
from single_wavelength_post import params
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'analysis/wp6_execution_20261103'
def save(name,r):
    r['passed']=all(r['checks'].values());(OUT/name).write_text(json.dumps(r,indent=2)+'\n');assert r['passed'];print(json.dumps(r,indent=2))
def static():
    checks={};src=ROOT/'model/transport_wp6_v1/source';old=ROOT/'model/transport_v1/source'
    for p in old.rglob('*'):
        if p.is_file() and p.name!='transportdiagnostics.hh':checks['unchanged/'+str(p.relative_to(old))]=p.read_bytes()==(src/p.relative_to(old)).read_bytes()
    for d in (OUT/'runs').iterdir():
        plan=json.loads((d/'plan.json').read_text());p=params(d/'params.input');b=params(ROOT/'analysis/single_wavelength_validation_20261103/runs'/('R24-'+plan['role']+'-full/params.input'))
        changed={k for k in p if p[k]!=b.get(k)};allowed={'Problem.Name','TimeLoop.TEnd','Transport.CutCoordinates'}
        if plan['stage']=='smoke':allowed.add('Transport.SolubilityMoleFraction')
        checks[d.name+'/allowed_input_changes']=changed<=allowed
        checks[d.name+'/bounded_time']=float(p['TimeLoop.TEnd'])<=272846.1275446968
        checks[d.name+'/input_sha']=hashlib.sha256((d/'params.input').read_bytes()).hexdigest()==plan['input_sha256']
        checks[d.name+'/profile_sha']=hashlib.sha256((d/'eigenprofile.dat').read_bytes()).hexdigest()==plan['profile_sha256']
        checks[d.name+'/seven_cuts']=len(p['Transport.CutCoordinates'].split())==7
    before={str(p):p.stat().st_mtime_ns for p in (ROOT/'analysis/wp6_bounded_nonlinear_plan_20261103').iterdir() if p.is_file()}
    proc=subprocess.run(['python',str(ROOT/'scripts/prepare_wp6_plan.py'),'--verify'],cwd=ROOT,capture_output=True)
    checks['readonly_verify']=proc.returncode==0 and all(Path(p).stat().st_mtime_ns==t for p,t in before.items())
    proc=subprocess.run(['python',str(ROOT/'scripts/prepare_wp6_plan.py'),'--recreate-to',str(ROOT/'analysis/wp6_bounded_nonlinear_plan_20261103')],cwd=ROOT,capture_output=True)
    checks['frozen_recreation_rejected']=proc.returncode!=0 and all(Path(p).stat().st_mtime_ns==t for p,t in before.items())
    import prepare_wp6_plan
    try:prepare_wp6_plan.main();blocked=False
    except RuntimeError:blocked=True
    checks['direct_generator_frozen_write_rejected']=blocked and all(Path(p).stat().st_mtime_ns==t for p,t in before.items())
    checks['new_build']=json.loads((OUT/'build/attempt_01/result.json').read_text())['build_exit']==0
    save('static_checks.json',dict(checks=checks))
def smoke():
    checks={}
    for name,step in [('initial',0),('crossing',1)]:
        d=OUT/'runs'/('WP6-'+name+'-smoke');f=json.loads((d/'finished.json').read_text());s=f['solver']['stop'];checks[name+'/returncode']=f['solver']['returncode']==0;checks[name+'/first_crossing_stop']=s['reason']=='surface_solubility_threshold' and s['accepted_steps']==step
        checks[name+'/independent_budgets']=json.loads((d/'post/independent_checks.json').read_text())['passed']
        rows=list(csv.DictReader((d/'surface_events.csv').open()));threshold=float(params(d/'params.input')['Transport.SolubilityMoleFraction'])
        checks[name+'/threshold_bracket']=float(rows[-1]['x_domain_max'])>=threshold and (step==0 or float(rows[-2]['x_domain_max'])<threshold)
    save('smoke_checks.json',dict(checks=checks))
if __name__=='__main__':p=argparse.ArgumentParser();p.add_argument('action',choices=['static','smoke']);a=p.parse_args();globals()[a.action]()
