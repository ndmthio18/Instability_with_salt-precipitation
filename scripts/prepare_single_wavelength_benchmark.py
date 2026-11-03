"""Pre-register two 10-step full-width timing controls on the same frodo hardware."""
import json,hashlib,re,datetime
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'analysis/single_wavelength_validation_20261103';BASE=ROOT/'analysis/pilot0_integrated_20261102/sal_migration/runs'
def main():
 spec=dict(utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),purpose='Task5 required fair hardware comparison after sal had no free physical cores; two 10-step original-width controls only, no new physical case',steps=10,T_s=10*227.3717729539139,output_interval=1,comparison='Internal elapsed time for the same first 10 accepted steps, including VTU each step; same frodo CPU cores and original ELF. Excludes startup. Does not claim measured full-window M1 speedup.',authorized_scope='command/5.md necessary limited-scale simulations and Task5 performance verification')
 (OUT/'timing_control_preregistration.json').write_text(json.dumps(spec,indent=2))
 for role in ['M0','M1']:
  d=OUT/'runs'/('F600-'+role+'-bench');d.mkdir(exist_ok=True);s=(BASE/role/'params.input').read_text();s=re.sub(r'^TEnd = .*$',f'TEnd = {spec["T_s"]:.17g}',s,flags=re.M);s=re.sub(r'^OutputInterval = .*$','OutputInterval = 1',s,flags=re.M);s=re.sub(r'^Name = .*$',f'Name = F600-{role}-bench',s,flags=re.M);(d/'params.input').write_text(s);(d/'eigenprofile.dat').write_bytes((BASE/'M1/eigenprofile.dat').read_bytes())
  plan=dict(id=d.name,role=role,nx_cells=600,nz_cells=200,DOF_vertices=601*201,scalar_unknowns=2*601*201,input_sha256=hashlib.sha256((d/'params.input').read_bytes()).hexdigest(),profile_sha256=hashlib.sha256((d/'eigenprofile.dat').read_bytes()).hexdigest(),elf_sha256=json.load(open(BASE/'M1/plan.json'))['executable']['sha256'],preregistration_sha256=hashlib.sha256((OUT/'preregistration.json').read_bytes()).hexdigest(),T_s=spec['T_s'],stage='bench',output_interval=1,additional_scope=spec['authorized_scope']);(d/'plan.json').write_text(json.dumps(plan,indent=2))
 print('registered two 10-step timing controls')
if __name__=='__main__':main()
