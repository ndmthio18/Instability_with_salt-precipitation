"""Prepare exactly two fresh exploratory runs from immutable verified inputs."""
import datetime,hashlib,json,shutil
from pathlib import Path
from single_wavelength_post import params
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'analysis/wp6_4T_exploratory_20261103'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
    assert not OUT.exists();OUT.mkdir()
    sources={'M0':ROOT/'analysis/wp6_review_20261103/runs/WP6-M0-dt_half','M1':ROOT/'analysis/wp6_timestep_20261103/runs/WP6-M1-dt_half'}
    reg=dict(id='20261103-D01',created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),authority='command/10.md',base_commit='c9241f2c1f769b7e7fa263b1488d1fa16bc2db21',scope='two dt/2 exploratory 4T legs from t=0 only',T_s=545692.2550893936,dt_s=113.68588647695695,wall_s_per_leg=3600,total_core_hours_cap=2,source_inputs={k:sha(p/'params.input') for k,p in sources.items()},elf_sha256=json.loads((sources['M0']/'plan.json').read_text())['elf_sha256'],
        retained_history=['original failed M0 absolute endpoint screen 1e-12','C01 engineering exploration separated from quantitative claims','unknown U_z U_x U_domain','no t_eff or surrogate labels'],
        analysis=dict(log_sigma_fit_windows_T=[[1,2],[2,3],[3,4]],T_base_s=136423.0637723484,energy_depth='90% and 99% of vertically integrated horizontal anomaly variance, measured from surface',harmonics='Neumann m=2 fundamental; m=4,6 are second/third harmonic; odd m diagnose symmetry',initial_structure_control='compare anomaly, normalized vertical profile and modal shape with t=0; inherited eigenprofile is not evidence of newly grown fingers',nonlinearity='descriptive relative amplitude/harmonic/shape evolution, not a universal nonlinear threshold'),forbidden=['8T','nx60','full domain','other cases','restart','source physics change'])
    (OUT/'amendment_20261103_D01.json').write_text(json.dumps(reg,indent=2)+'\n');checks={}
    for role,src in sources.items():
        d=OUT/'runs'/('WP6-'+role+'-explore4T');d.mkdir(parents=True)
        s=(src/'params.input').read_text().replace('Name = WP6-'+role+'-dt_half','Name = '+d.name).replace('TEnd = 272846.1275446968','TEnd = 545692.2550893936');(d/'params.input').write_text(s);shutil.copyfile(src/'eigenprofile.dat',d/'eigenprofile.dat')
        a,b=params(d/'params.input'),params(src/'params.input');delta={k for k in a if a[k]!=b[k]};assert delta=={'Problem.Name','TimeLoop.TEnd'}
        p=json.loads((src/'plan.json').read_text());p.update(id=d.name,stage='explore4T',T_s=reg['T_s'],input_sha256=sha(d/'params.input'),preregistration_sha256=sha(OUT/'amendment_20261103_D01.json'),hard_wall_limit_s=3600);(d/'plan.json').write_text(json.dumps(p,indent=2));checks[role]=dict(allowed_input_diff=sorted(delta),profile_unchanged=sha(d/'eigenprofile.dat')==p['profile_sha256'])
    (OUT/'static_checks.json').write_text(json.dumps(dict(passed=True,legs=checks,exactly_two_runs=True),indent=2))
    campaign=(ROOT/'scripts/wp6_timestep_campaign.py').read_text().replace('wp6_timestep_20261103','wp6_4T_exploratory_20261103').replace('amendment_20261103_C01','amendment_20261103_D01').replace("choices=['dt_half']","choices=['explore4T']")
    campaign=campaign.replace("'single_wavelength_post.py']","'single_wavelength_post.py','wp6_timestep_manager.py']").replace("'wp6_manager.py':'wp6_timestep_manager.py'","'wp6_manager.py':'wp6_4t_manager.py'")
    campaign=campaign.replace('One authorized M1 half-step campaign on frodo; never repeat launch.','Two authorized exploratory 4T legs on frodo; never repeat launch.')
    (ROOT/'scripts/wp6_4t_campaign.py').write_text(campaign)
    print('D01 registered; two inputs differ only in name and horizon from their verified fine-step references.')
if __name__=='__main__':main()
