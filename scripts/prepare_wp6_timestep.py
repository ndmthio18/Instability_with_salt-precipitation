"""Register exactly one M1 half-step run; refuses to replace existing evidence."""
import datetime,hashlib,json,shutil
from pathlib import Path
from single_wavelength_post import params
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'analysis/wp6_timestep_20261103'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    OUT.mkdir(exist_ok=True)
    assert not (OUT/'frodo_launch_dt_half.json').exists(), 'cannot prepare after launch'
    base=ROOT/'analysis/wp6_execution_20261103/runs/WP6-M1-extend'
    m0=ROOT/'analysis/wp6_review_20261103/runs/WP6-M0-dt_half'
    reg=dict(id='20261103-C01',created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),authority='command/9.md',base_commit='7db1dd3f5b9ee98802d46c2d5fecfaf913c5d68d',allowed_new_solver_runs=1,role='M1',T_s=272846.1275446968,dt_s=113.68588647695695,wall_limit_s=3600,
        source_M1_input_sha256=sha(base/'params.input'),matched_M0_input_sha256=sha(m0/'params.input'),elf_sha256=json.loads((base/'plan.json').read_text())['elf_sha256'],
        original_endpoint_gate=dict(threshold=1e-12,baseline_status='failed_retained',changed=False),
        exploratory_A=dict(scope='growth/structure exploration only; not quantitative transport validation',registration_context='after known M0 noise audit, before new refined M1 outcome; engineering screening, not physical accuracy theorem',
        relative_sigma_fourier_amplification_tolerance=.05,vertical_sigma_profile_relative_L2_tolerance=.05,relative_surface_mean_max_tolerance=.001,
        detection_time_tolerance='max(2 coarse dt, 2% coarse detection time)',background_sigma_over_M1_sigma_limit=.001,
        justification='5% distinguishes robust exploratory growth/shape from protocol-dependent changes; 0.1% surface-state agreement is stricter for mean background; noise must be at least 1000 times smaller than resolved perturbation, following the existing relative symmetry margin. These criteria do not bound tiny paired flux error.',
        required=['independent hashes and budget checks','strict within-protocol pairing','M1 original relative symmetry screen','no extrapolation beyond solubility','retain absolute M0 gate as historical failed']),
        quantitative_B=dict(required=['paired timestep sensitivity does not cover signal','vertical and horizontal and domain errors bounded','solubility time intervals informative'],unknown_errors=['U_z','U_x','U_domain']),
        timing_rule='Use exactly shared coarse/fine accepted times where available. Endpoint gaps must be reported and bounded with local conservative flux; no interpolation of cumulative flux.',
        forbidden=['4T launch','8T launch','additional solver runs','amplitude or physics change'])
    registration=OUT/'amendment_20261103_C01.json'
    if registration.exists():
        old=json.loads(registration.read_text());reg['created_utc']=old['created_utc'];assert old==reg
    else:registration.write_text(json.dumps(reg,indent=2)+'\n')
    d=OUT/'runs/WP6-M1-dt_half';d.mkdir(parents=True,exist_ok=True)
    s=(base/'params.input').read_text().replace('Name = WP6-M1-extend','Name = WP6-M1-dt_half').replace('DtInitial = 227.3717729539139','DtInitial = 113.68588647695695').replace('MaxTimeStepSize = 227.3717729539139','MaxTimeStepSize = 113.68588647695695')
    (d/'params.input').write_text(s);shutil.copyfile(base/'eigenprofile.dat',d/'eigenprofile.dat')
    p=json.loads((base/'plan.json').read_text());p.update(id=d.name,stage='dt_half',input_sha256=sha(d/'params.input'),preregistration_sha256=sha(OUT/'amendment_20261103_C01.json'),hard_wall_limit_s=3600);(d/'plan.json').write_text(json.dumps(p,indent=2))
    a,b=params(d/'params.input'),params(base/'params.input');delta={k for k in a if a[k]!=b[k]};assert delta=={'Problem.Name','TimeLoop.DtInitial','TimeLoop.MaxTimeStepSize'}
    c=params(m0/'params.input');role_delta={k for k in a.keys()|c.keys() if a.get(k)!=c.get(k)};assert role_delta=={'Problem.Name','Problem.Amplitude','Problem.PerturbationType','Problem.EigenprofilePath','Problem.StandardDeviationSaltMoleFraction'}
    assert sha(d/'eigenprofile.dat')==sha(m0/'eigenprofile.dat')==p['profile_sha256']
    (OUT/'static_checks.json').write_text(json.dumps(dict(passed=True,coarse_to_fine_diff=sorted(delta),matched_M0_to_M1_diff=sorted(role_delta),same_ELF=True,same_profile=True,hard_wall_s=3600),indent=2))
    worker=(ROOT/'scripts/wp6_worker.py').read_text().replace('timeout=43100','timeout=3500')
    (ROOT/'scripts/wp6_timestep_worker.py').write_text(worker)
    campaign=(ROOT/'scripts/wp6_review_campaign.py').read_text().replace('wp6_review_20261103','wp6_timestep_20261103').replace('amendment_20261103_B01','amendment_20261103_C01')
    campaign=campaign.replace("(ROOT/'scripts'/name).read_bytes()","(ROOT/'scripts'/({'wp6_worker.py':'wp6_timestep_worker.py','wp6_manager.py':'wp6_timestep_manager.py'}.get(name,name))).read_bytes()")
    campaign=campaign.replace("['state','resources','finished']","['state','resources','finished','wall_limit_result']")
    campaign=campaign.replace("p=r/(STAGE+'_manager_'+n+'.json')","p=r/(STAGE+('_wall_limit_result.json' if n=='wall_limit_result' else '_manager_'+n+'.json'))")
    campaign=campaign.replace('"""Frodo migration after proven cancellation of the unstarted sal array."""','"""One authorized M1 half-step campaign on frodo; never repeat launch."""')
    (ROOT/'scripts/wp6_timestep_campaign.py').write_text(campaign)
    print('Registered C01 and one M1 dt/2 run; source physics unchanged.')
if __name__=='__main__':main()
