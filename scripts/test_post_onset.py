"""Independent analytic fixtures plus scientific-data acceptance gates for WP4."""
import unittest
from analyze_post_onset import *
from analyze_post_onset_fields import penetration,spectrum
class MetricsTests(unittest.TestCase):
 def test_irregular_exponential_slope(self):
  t=np.r_[0.,np.geomspace(1,100,50)];s=2e-7*np.exp(.003*t);f=fit(t,s,0,100);self.assertAlmostEqual(f['slope_s_inv'],.003,places=13);self.assertAlmostEqual(f['r2'],1,places=13)
 def test_transient_peak_not_hidden_by_endpoint(self):
  t=np.arange(101.);s=np.exp(-abs(t-60)/30)*1e-6;a=np.c_[t,np.full(101,.01),s];r=metrics(a,20,s[20]);self.assertGreater(r['G_max'],r['G_end']);self.assertEqual(r['tail_shape'],'rise_then_fall');self.assertFalse(r['sustained_growth_two_windows'])
 def test_missing_detection_stays_unresolved(self):
  a=np.c_[np.arange(20.),np.full(20,.01),np.full(20,1e-6)];r=metrics(a,np.nan,np.nan);self.assertNotIn('G_max',r);self.assertAlmostEqual(r['A_rel_end'],1e-4,places=15)
 def test_floor_sensitivity(self):
  a=np.c_[np.arange(20.),np.full(20,.01),np.geomspace(1e-16,1e-12,20)];r=metrics(a,1,1e-16);self.assertTrue(r['floating_floor_flag']);self.assertGreater(r['G_end'],r['G_end_floor_1e-10'])
 def test_invalid_time_rejected(self):
  with self.assertRaises(AssertionError):metrics(np.array([[1,.01,1e-6],[1,.01,2e-6]]),1,1e-6)
 def test_constant_depth_energy(self):
  z=np.linspace(0,.2,201);r=penetration(z,np.ones(201)*1e-6);self.assertAlmostEqual(r['depth_centroid_H'],.5);self.assertLess(abs(r['depth90_H']-.9),.006)
 def test_analytic_two_mode_spectrum(self):
  xx=np.linspace(0,.6,601);z=np.linspace(0,.2,21);x=.01+np.tile(1e-6*(np.cos(2*np.pi*3*xx/.6)+.1*np.cos(2*np.pi*6*xx/.6)),(21,1));r=spectrum(x,z,3);self.assertAlmostEqual(r['spectral_fundamental_energy_fraction'],1/1.01,places=10);self.assertAlmostEqual(r['spectral_second_harmonic_energy_fraction'],.01/1.01,places=10)
 def test_full_inventory_and_separate_denominators(self):
  r=readcsv(OUT/'run_metrics.csv');m=[x for x in r if x['source_group']=='main'];w=[x for x in r if x['source_group']=='wp1'];self.assertEqual((len(m),len(set(x['case_id'] for x in m)),len(w)),(436,253,84));self.assertEqual(sum(x['series_status']=='hash_verified' for x in m),320);self.assertTrue(all(x['dumux_t_detect_s']=='' for x in m if x['dumux_detection_class']!='resolved_onset'))
 def test_recomputed_detection_agrees(self):self.assertEqual(json.load(open(OUT/'metric_checks.json'))['class_disagreements'],[])
 def test_dense_time_and_source_reuse(self):
  r=readcsv(OUT/'dense_depth_timeseries.csv');self.assertEqual(len(r),1414)
  for rid in set(x['id'] for x in r):self.assertTrue(np.all(np.diff([num(x['time_s']) for x in r if x['id']==rid])>0))
  for q in json.load(open(OUT/'dense_reuse_hashes.json')):self.assertEqual(hashlib.sha256((ROOT/q['path']).read_bytes()).hexdigest(),q['sha256'])
 def test_quarantine_and_initial_only(self):
  r=readcsv(OUT/'sparse_field_metrics.csv');self.assertEqual(len(r),18);self.assertEqual(set(x['case_id'] for x in r if 'quarantined' in x['physical_field_status']),{'35733','19674','38852'});self.assertEqual(sum(x['case_id']=='30858' for x in r),1)
 def test_no_non_linear_candidate_A(self):
  r=readcsv(OUT/'candidate_shortlist.csv');self.assertEqual(len(set(x['case_id'] for x in r)),9);self.assertFalse(any(x['candidate_category']=='A' for x in r));self.assertTrue(all(x['nonlinear_evidence']=='not_established' for x in r))
if __name__=='__main__':
 suite=unittest.defaultTestLoader.loadTestsFromTestCase(MetricsTests);result=unittest.TextTestRunner(verbosity=2).run(suite);savejson(OUT/'test_results.json',dict(tests=result.testsRun,failures=len(result.failures),errors=len(result.errors),passed=result.wasSuccessful()));raise SystemExit(0 if result.wasSuccessful() else 1)
