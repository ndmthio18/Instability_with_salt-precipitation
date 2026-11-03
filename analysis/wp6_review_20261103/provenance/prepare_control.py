from pathlib import Path
import json,hashlib,datetime,shutil
root=Path.cwd();out=root/'analysis/wp6_review_20261103';out.mkdir(exist_ok=True)
base=root/'analysis/wp6_execution_20261103';src=base/'runs/WP6-M0-extend';d=out/'runs/WP6-M0-dt_half';d.mkdir(parents=True)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
amend=dict(id='20261103-B01',created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),authority='command/8.md Task1',base_commit='6e8e03428779dd1f64134cc5878c9f9ebb03e927',purpose='one bounded M0 dt/2 control to separate timestep dependence from unresolved background symmetry mechanisms',allowed_runs=1,T_s=272846.1275446968,dt_s=113.68588647695695,per_leg_wall_s=43200,baseline_threshold_unchanged=1e-12,new_symmetry_threshold=None,forbidden=['4T','8T','full domain','nx60','second diagnostic solver run'],baseline_input_sha256=sha(src/'params.input'),elf_sha256=json.loads((src/'plan.json').read_text())['elf_sha256'])
with (out/'amendment_20261103_B01.json').open('x') as f:json.dump(amend,f,indent=2)
s=(src/'params.input').read_text().replace('Name = WP6-M0-extend','Name = WP6-M0-dt_half').replace('DtInitial = 227.3717729539139','DtInitial = 113.68588647695695').replace('MaxTimeStepSize = 227.3717729539139','MaxTimeStepSize = 113.68588647695695')
(d/'params.input').write_text(s,encoding='utf-8');shutil.copyfile(src/'eigenprofile.dat',d/'eigenprofile.dat')
p=json.loads((src/'plan.json').read_text());p.update(id=d.name,stage='dt_half',input_sha256=sha(d/'params.input'),preregistration_sha256=sha(out/'amendment_20261103_B01.json'));(d/'plan.json').write_text(json.dumps(p,indent=2))
s=(root/'scripts/wp6_campaign.py').read_text().replace('analysis/wp6_execution_20261103','analysis/wp6_review_20261103').replace('tmp/wp6_execution_20261103','tmp/wp6_review_20261103').replace('/temp/zhixin/paper2_wp6_execution_20261103_r01','/temp/zhixin/paper2_wp6_review_20261103_r01').replace('amendment_20261103_A01.json','amendment_20261103_B01.json').replace("choices=['smoke','bridge','extend']","choices=['dt_half']")
(root/'scripts/wp6_review_campaign.py').write_text(s,encoding='utf-8')
import sys
sys.path.insert(0,str(root/'scripts'))
from single_wavelength_post import params
x=params(d/'params.input');y=params(src/'params.input');delta={k for k in x if x[k]!=y[k]};assert delta=={'Problem.Name','TimeLoop.DtInitial','TimeLoop.MaxTimeStepSize'}
(out/'static_checks.json').write_text(json.dumps(dict(passed=True,allowed_input_diff=sorted(delta),same_ELF=True,same_grid=True,same_threshold=True,only_one_run=True),indent=2))
print('Registered one M0-only half-step control; no threshold amendment.')
