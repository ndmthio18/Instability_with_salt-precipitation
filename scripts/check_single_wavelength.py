"""Analytic quadrature/phase and permitted-input-change checks for WP5."""
import json
from pathlib import Path
import numpy as np
from single_wavelength_post import params,widths
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'analysis/single_wavelength_validation_20261103'
def main():
 checks={}
 for d in (OUT/'runs').iterdir():
  p=params(d/'params.input');role=d.name.split('-')[1];b=params(ROOT/'analysis/pilot0_integrated_20261102/sal_migration/runs'/role/'params.input');changed={k for k in p if p[k]!=b.get(k)}
  checks[d.name+'/only_allowed_changes']=changed<={'Grid.Positions0','Grid.Cells0','TimeLoop.TEnd','TimeLoop.OutputInterval','Problem.Name','Problem.WaveNumber'}
  W=float(p['Grid.Positions0'].split()[-1]);n=int(p['Grid.Cells0']);checks[d.name+'/wave_number']=abs(float(p['Problem.WaveNumber'])/W-28/float(b['Grid.Positions0'].split()[-1]))<1e-12
  x=np.linspace(0,W,n+1);wx=widths(x);v=np.cos(2*np.pi*float(p['Problem.WaveNumber'])*x/W);checks[d.name+'/zero_weighted_mean']=abs(v@wx/W)<1e-13;checks[d.name+'/weighted_variance']=abs(v**2@wx/W-.5)<1e-13;checks[d.name+'/reflection']=max(abs(v-v[::-1]))<1e-12;checks[d.name+'/initial_salt_inventory_cancel']=abs(np.sum(1e-6*v*wx)/W)<1e-19
 checks['original_ELF_unchanged']=json.load(open(OUT/'runs/R24-M1-short/plan.json'))['elf_sha256']=='d78ce5122cf69e4dee5a4d12bcc726779959f304b343246009cb45bb37e894cc'
 checks={k:bool(v) for k,v in checks.items()};(OUT/'configuration_checks.json').write_text(json.dumps(dict(checks=checks,passed=all(checks.values()),note='First serialization attempt before submission failed on numpy bool; recomputed and recorded after short submission. Preregistration and physical inputs were frozen before submission; no solver retry.'),indent=2));assert all(checks.values());print('PASS',len(checks),'configuration checks')
if __name__=='__main__':main()
