"""WP5 analytic and cross-output scientific checks, separate from calibration thresholds."""
import json,hashlib,unittest
from pathlib import Path
import numpy as np
from single_wavelength_post import widths,params
from analyze_post_onset import readcsv,num,savejson,writecsv
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'analysis/single_wavelength_validation_20261103';BASE=ROOT/'analysis/pilot0_integrated_20261102/sal_migration/runs'
class Checks(unittest.TestCase):
 def test_neumann_is_not_periodic(self):
  W=.6/28;x=np.linspace(0,W,49);v=np.cos(np.pi*x/W);self.assertGreater(abs(v[0]-v[-1]),1.9);self.assertLess(abs(np.sin(np.pi)),1e-14)
 def test_weighted_quadrature_invariant(self):
  for N,mode,W in [(600,28,.6000000000000001),(24,1,.6000000000000001/28),(48,1,.6000000000000001/28)]:
   x=np.linspace(0,W,N+1);w=widths(x);v=np.cos(2*np.pi*mode*x/W);self.assertAlmostEqual(v@w/W,0,places=13);self.assertAlmostEqual(v**2@w/W,.5,places=13);self.assertAlmostEqual(np.std(v)**2,(N/2+1)/(N+1)-1/(N+1)**2,places=13)
 def test_saved_initial_statistics(self):
  rows=[]
  for name in ['R24-M1-short','R48-M1-short']:
   d=OUT/'runs'/name;sig=np.loadtxt(next(d.glob('*sigmatop.csv')));w=readcsv(d/'post/weighted_diagnostics.csv')[0];N=int(params(d/'params.input')['Grid.Cells0']);self.assertAlmostEqual(sig[0,1],.0036+1e-6/(N+1),places=13);self.assertAlmostEqual(sig[0,2],1e-6*np.sqrt((N/2+1)/(N+1)-1/(N+1)**2),places=17);self.assertAlmostEqual(num(w['weighted_top_sigma']),1e-6/np.sqrt(2),places=17);self.assertAlmostEqual(num(w['weighted_projection']),1e-6,places=17)
   rows.append(dict(model=name,node_mean=sig[0,1],weighted_mean=num(w['weighted_top_mean']),node_sigma=sig[0,2],weighted_sigma=num(w['weighted_top_sigma']),node_fourier=sig[0,4],weighted_fourier=num(w['weighted_fourier']),node_projection=sig[0,5],weighted_projection=num(w['weighted_projection'])))
  writecsv(OUT/'initial_sampling_bias.csv',rows)
 def test_benchmark_trajectory_same_as_reference(self):
  for role in ['M0','M1']:
   sig=np.loadtxt(next((OUT/'runs'/f'F600-{role}-bench').glob('*sigmatop.csv')));ref=np.loadtxt(next((BASE/role).glob('*sigmatop.csv')))
   self.assertEqual(sig.shape[0],11);self.assertTrue(np.allclose(sig[:,:3],ref[:11,:3],rtol=1e-10,atol=1e-13))
 def test_input_and_elf_hashes(self):
  for d in (OUT/'runs').iterdir():
   p=json.load(open(d/'plan.json'));self.assertEqual(hashlib.sha256((d/'params.input').read_bytes()).hexdigest(),p['input_sha256']);self.assertEqual(p['elf_sha256'],'d78ce5122cf69e4dee5a4d12bcc726779959f304b343246009cb45bb37e894cc')
 def test_completed_no_reduced_dt(self):
  for d in (OUT/'runs').iterdir():
   p=json.load(open(d/'plan.json'));fin=json.load(open(d/'finished.json'));self.assertEqual(fin['solver']['returncode'],0);rows=readcsv(d/'surface_events.csv');self.assertTrue(all(abs(num(r['dt_s'])-227.3717729539139)<1e-7 for r in rows[1:]));self.assertLess(abs(num(rows[-1]['time_s'])-p['T_s']),1e-7)
 def test_pair_identity(self):
  for stage in ['short','full']:
   rows=readcsv(OUT/(stage+'_paired_transport.csv'));self.assertLess(max(abs(num(r['identity_error_m'])) for r in rows),1e-18);self.assertTrue(all(r['interpretation']=='unresolved' for r in rows))
 def test_reference_files_untouched(self):
  for name,h in json.load(open(OUT/'reference_hashes.json')).items():self.assertEqual(hashlib.sha256((ROOT/name).read_bytes()).hexdigest(),h)
 def test_gates(self):
  self.assertTrue(json.load(open(OUT/'short_validation.json'))['advance_to_full']);self.assertTrue(json.load(open(OUT/'full_validation.json'))['numerical_passed'])
if __name__=='__main__':
 result=unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromTestCase(Checks));savejson(OUT/'test_results.json',dict(tests=result.testsRun,failures=len(result.failures),errors=len(result.errors),passed=result.wasSuccessful()));raise SystemExit(0 if result.wasSuccessful() else 1)
