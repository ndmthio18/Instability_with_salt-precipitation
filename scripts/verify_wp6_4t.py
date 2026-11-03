"""Verify the exploratory pair, its sources, archives and unchanged history."""
import csv,hashlib,json,py_compile,subprocess,sys
from pathlib import Path
import fitz
from wp6_plan_integrity import verify
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'analysis/wp6_4T_exploratory_20261103';CACHE=ROOT/'tmp/wp6_4T_exploratory_20261103'
SCRIPTS=['prepare_wp6_4t.py','wp6_4t_manager.py','wp6_4t_campaign.py','wp6_4t_spatial.py','collect_wp6_4t_spatial.py','wp6_4t_metrics.py','test_wp6_4t.py','analyze_wp6_4t.py','plot_wp6_4t.py','verify_wp6_4t.py']
def load(p):return json.loads(p.read_text())
def sha(p):
    with p.open('rb') as f:return hashlib.file_digest(f,'sha256').hexdigest()
def main():
    verify();checks={}
    for name in ['single_wavelength_validation_20261103','wp6_execution_20261103','wp6_review_20261103','wp6_timestep_20261103']:
        for rel,r in load(ROOT/'analysis'/name/'artifact_manifest.json').items():assert sha(ROOT/rel)==r['sha256'],rel
        checks[name+'_unchanged']=True
    for name,s in load(ROOT/'model/transport_wp6_v1/source_manifest.json').items():assert sha(ROOT/'model/transport_wp6_v1/source'/name)==s
    checks['physics_sources_unchanged']=True
    assert len(list((OUT/'runs').iterdir()))==2
    state=load(OUT/'frodo_explore4T_status.json')['manager'];cpus=state['resources']['cpus'];assert len(set(cpus))==2
    done=state['finished'];assert len(done)==2 and all(v==0 for v in done.values()) and not state['state']['active']
    count=0;allocated=0.;cpu_s=0.;validation={}
    for d in (OUT/'runs').iterdir():
        f=load(d/'finished.json');a=load(d/'started.json');p=load(d/'plan.json')
        assert a['elf_sha256']==p['elf_sha256']==load(OUT/'amendment_20261103_D01.json')['elf_sha256']
        assert sha(d/'params.input')==p['input_sha256'] and sha(d/'eigenprofile.dat')==p['profile_sha256']
        assert sha(OUT/'amendment_20261103_D01.json')==p['preregistration_sha256']
        assert load(OUT/'amendment_20261103_D01.json')['created_utc']<a['utc']
        assert a['affinity']==[cpus[0 if p['role']=='M0' else 1]]
        for name,r in f['files'].items():
            path=d/name
            if not path.exists():path=CACHE/'runs'/d.name/name
            assert sha(path)==r['sha256'],str(path)
        ar=load(OUT/'explore4T_archive.json')[d.name];assert sha(Path(ar['local_path']))==ar['sha256'];assert ar['verified_vtu']==len(load(d/'post/field_hashes.json'));count+=ar['verified_vtu']
        lim=state['wall_limit_result'][d.name];assert lim['hard_limit_s']==3600 and lim['wall_s']<3600 and not lim['timed_out'];allocated+=lim['wall_s'];cpu_s+=f['solver']['cpu_s']
        validation[d.name]=dict(accepted_steps=f['solver']['stop']['accepted_steps'],accepted_time_s=f['solver']['stop']['last_accepted_time_s'],stop_reason=f['solver']['stop']['reason'],core=cpus[0 if p['role']=='M0' else 1])
    assert allocated<=7200 and cpu_s<=7200
    checks['two_independent_legs']=validation;checks['archive_verified_vtu']=count;checks['allocated_core_hours_including_post']=allocated/3600;checks['solver_cpu_hours']=cpu_s/3600
    assert allocated+2*900<=7200
    checks['including_additional_spatial_analysis_core_hour_upper_bound']=(allocated+2*900)/3600
    checks['additional_spatial_analysis_cpu_s']=None
    checks['additional_spatial_analysis_budget_note']='Two successful post-only calls, each limited to 900 s on one freed physical core; CPU usage not separately recorded'
    mapping={'wp6_manager.py':'wp6_4t_manager.py','wp6_worker.py':'wp6_timestep_worker.py','wp6_timestep_manager.py':'wp6_timestep_manager.py','wp6_post.py':'wp6_post.py','single_wavelength_post.py':'single_wavelength_post.py'}
    for name,s in load(OUT/'remote_execution_source_hashes.json').items():assert sha(ROOT/'scripts'/mapping[name])==s
    for folder in ['spatial','spatial_M0']:
        spatial=load(OUT/folder/'spatial_audit.json');assert spatial['script_sha256']==sha(ROOT/'scripts/wp6_4t_spatial.py') and spatial['metrics_helper_sha256']==sha(ROOT/'scripts/wp6_4t_metrics.py')
        for name,r in load(OUT/folder/'transfer.json').items():assert sha((CACHE if name.endswith('.npz') else OUT)/folder/name)==r['sha256']
    checks['spatial_source_and_products_verified']=True
    tests=[]
    for name,count in [('test_wp6_paired_stopping.py',22),('test_wp6_timestep.py',9),('test_wp6_4t.py',10)]:
        r=subprocess.run([sys.executable,str(ROOT/'scripts'/name),'-v'],capture_output=True,text=True,encoding='utf-8');assert r.returncode==0,r.stderr;tests.append(dict(script=name,count=count,passed=True,output=r.stdout+r.stderr))
    (OUT/'test_results.json').write_text(json.dumps(tests,indent=2));checks['tests_passed']=41
    for name in SCRIPTS:py_compile.compile(str(ROOT/'scripts'/name),doraise=True)
    acceptance=load(OUT/'numerical_acceptance.json');assert acceptance['checks']['core_numerical_acceptance']
    assert acceptance['historical_M0_absolute_screen']=='failed_retained' and acceptance['t_eff'] is None
    checks['numerical_acceptance_with_separate_symmetry_flags']=acceptance['checks']
    pair=load(OUT/'paired/async_stopping_checks.json');checks['shared_prethreshold_states']=len(pair['pair']['indices']);checks['event_difference']=pair['pair']['surface_time_difference']
    assert len(list(OUT.glob('Figure*.pdf')))==5 and len(list(OUT.glob('Figure*.png')))==5
    for p in OUT.glob('Figure*.pdf'):
        with fitz.open(p) as doc:assert len(doc)==1 and len(doc[0].get_text())>100
    checks['figure_QA']='Five PNG groups visually inspected; five PDF pages checked for selectable labels'
    checks['source_base_commit']='c9241f2c1f769b7e7fa263b1488d1fa16bc2db21';checks['execution_scope']='4T decision point only; no 8T, other cases or wavelength domains'
    (OUT/'delivery_checks.json').write_text(json.dumps(checks,indent=2)+'\n')
    files=[p for p in OUT.rglob('*') if p.is_file() and p.name!='artifact_manifest.json' and p.suffix not in ['.vtu','.pvd','.npz']]+[ROOT/'scripts'/name for name in SCRIPTS]+[ROOT/'command/10.md',ROOT/'model/transport_wp6_v1/source_manifest.json']
    (OUT/'artifact_manifest.json').write_text(json.dumps({p.relative_to(ROOT).as_posix():dict(sha256=sha(p),bytes=p.stat().st_size) for p in sorted(files)},indent=2)+'\n');print(json.dumps(checks,indent=2))
if __name__=='__main__':main()
