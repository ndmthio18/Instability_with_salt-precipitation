from pathlib import Path
import json,hashlib,datetime,shutil
root=Path.cwd();out=root/'analysis/wp6_execution_20261103';out.mkdir(exist_ok=True)
reg=json.loads((root/'analysis/wp6_bounded_nonlinear_plan_20261103/preregistration.json').read_text());build=json.loads((out/'build/attempt_01/result.json').read_text())
amend=dict(id='20261103-A01',created_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),authority='command/7.md',prior_commit='314d75d9e96be24ee956b5677cb22c8e7bc0b989',frozen_preregistration_sha256=hashlib.sha256((root/'analysis/wp6_bounded_nonlinear_plan_20261103/preregistration.json').read_bytes()).hexdigest(),authorized_stages=['two micro event fixtures','T bridge M0/M1','2T M0/M1 if bridge passes'],forbidden=['4T','8T','full domain','multiwavelength execution','nx60 execution'],diagnostic_schema=2,cut_coordinates=[v['actual_z_m'] for v in reg['observation_planes']],elf_sha256=build['executable']['sha256'],per_leg_wall_s=43200,reference_restoration='hash-verified original field archives; regenerate ignored cache only',immutable_plan_policy='default verification; recreation only in fresh tmp directory; changes require dated amendment')
with (out/'amendment_20261103_A01.json').open('x') as f:json.dump(amend,f,indent=2)
for stage,mult,roles in [('smoke',2/600,['initial','crossing']),('bridge',1,['M0','M1']),('extend',2,['M0','M1'])]:
 for role in roles:
  actual='M0' if role in ['M0','initial'] else 'M1';name='WP6-'+role+'-'+stage;d=out/'runs'/name;d.mkdir(parents=True)
  src=root/'analysis/single_wavelength_validation_20261103/runs'/('R24-'+actual+'-full')
  s=(src/'params.input').read_text().replace('Name = R24-'+actual+'-full','Name = '+name).replace('TEnd = 136423.0637723484','TEnd = '+str(reg['T_s']*mult))
  s+='\nCutCoordinates = '+' '.join(format(v['actual_z_m'],'.17g') for v in reg['observation_planes'])+'\n'
  if role=='initial':s=s.replace('SolubilityMoleFraction = 0.1','SolubilityMoleFraction = 0.0035999')
  if role=='crossing':s=s.replace('SolubilityMoleFraction = 0.1','SolubilityMoleFraction = 0.00361')
  (d/'params.input').write_text(s,encoding='utf-8');shutil.copyfile(src/'eigenprofile.dat',d/'eigenprofile.dat')
  sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
  plan=dict(id=name,role=actual,stage=stage,T_s=reg['T_s']*mult,input_sha256=sha(d/'params.input'),profile_sha256=sha(d/'eigenprofile.dat'),elf_sha256=build['executable']['sha256'],preregistration_sha256=sha(out/'amendment_20261103_A01.json'),nx_cells=24,nz_cells=200,DOF_vertices=5025,scalar_unknowns=10050)
  (d/'plan.json').write_text(json.dumps(plan,indent=2))
s=(root/'scripts/single_wavelength_frodo_worker.py').read_text().replace('WP5','WP6').replace("root/'single_wavelength_post.py'","root/'wp6_post.py'")
s=s.replace("if rc or not stop or stop['reason']!='requested_horizon':return 1","if rc or not stop or stop['reason'] not in ['requested_horizon','surface_solubility_threshold','interior_solubility_threshold']:return 1")
s=s.replace("assert abs(stop['last_accepted_time_s']-p['T_s'])<1e-7","assert stop['last_accepted_time_s']<=p['T_s']+1e-7\n if stop['reason']=='requested_horizon':assert abs(stop['last_accepted_time_s']-p['T_s'])<1e-7")
s=s.replace('>8*1024**3','>200*1024**3')
(root/'scripts/wp6_worker.py').write_text(s,encoding='utf-8')
s=(root/'scripts/single_wavelength_manager.py').read_text().replace('WP5','WP6').replace('single_wavelength_frodo_worker.py','wp6_worker.py');(root/'scripts/wp6_manager.py').write_text(s,encoding='utf-8')
print('Fresh inputs and immutable execution amendment written.')
