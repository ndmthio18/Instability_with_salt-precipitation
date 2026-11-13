"""Scientific/data integrity checks; does not substitute for numerical convergence."""
import unittest,json
from pathlib import Path
import numpy as np
import pandas as pd
from analyze import ROOT,read,PHI,CR
class AuditTests(unittest.TestCase):
 def test_theoretical_constant(self):self.assertAlmostEqual(PHI,.8898899218857486,places=13)
 def test_no_theory_in_evolution(self):
  for name in ['two_front.hh','two_bdf.hh','initialize.hh']:
   s=(ROOT/name).read_text();self.assertNotIn('phi_core',s);self.assertNotIn('.889889',s)
 def test_histories_monotone_and_geometric(self):
  for p in (ROOT/'raw').glob('*_history.csv'):
   h=pd.read_csv(p)
   with self.subTest(case=p.name):
    self.assertTrue(np.all(np.diff(h.t)>0));self.assertLess(abs(h.hc-(h.zt-h.zb)).max(),1e-15)
    self.assertLess(abs(h.vnet-(h.vt-h.vb)).max(),1e-20)
 def test_error_decomposition(self):
  for lab in ['u06','u09','x008']:
   c=f'setup_p0_ld_{lab}_r1e-10'
   if not (ROOT/'raw'/f'{c}_history.csv').exists():c=c.removeprefix('setup_')
   h=read(c);self.assertLess(abs(h.delta_net-h.delta_top+h.delta_bottom).max(),1e-22)
   self.assertLess(abs(h.vnet-h.mass_term-h.structure_term).max(),1e-22)
   self.assertLess(abs(h.delta_net-h.supply_term-h.storage_term-h.density_term-h.structure_term).max(),1e-22)
 def test_bottom_jump_independent_integrals(self):
  for p in (ROOT/'raw').glob('*_history.csv'):
   h=pd.read_csv(p)
   with self.subTest(case=p.name):
    self.assertLess(abs(h.Q_SPLUS-h.Q_SMINUS-h.Q_DB).max(),2e-9)
    self.assertLess(abs(h.Q_MPLUS-h.Q_MMINUS-h.Q_DB).max(),2e-9)
    self.assertLess(abs(h.liquid_salt-h.liquid_salt.iloc[0]-h.Q_SOIL_SALT+h.Q_SMINUS).max(),2e-9)
 def test_original_failures_retained(self):
  for name in ['legacy_r9','legacy_r10']:
   status=(ROOT/'raw'/f'{name}_status.txt').read_text();self.assertIn('solver_stop',status)
   self.assertTrue((ROOT/'raw'/f'{name}_failure_state.csv').exists())
 def test_forcing_support(self):
  f=pd.read_csv(ROOT/'forcing.csv');self.assertAlmostEqual(float(f.iloc[-1,0]),160951.5)
  for p in (ROOT/'raw').glob('*x008*_history.csv'):self.assertLessEqual(pd.read_csv(p).t.max(),160951.5)
if __name__=='__main__':unittest.main(verbosity=2)
