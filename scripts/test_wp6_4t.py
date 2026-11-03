import unittest
import numpy as np
from wp6_4t_metrics import log_growth_fit,energy_depth,neumann_spectrum

class GrowthTests(unittest.TestCase):
    def test_exponential_rate(self):
        t=np.linspace(0,100,101);r=log_growth_fit(t,1e-6*np.exp(.01*t),0,100)
        self.assertAlmostEqual(r['slope_per_s'],.01);self.assertEqual(r['status'],'complete_window')
    def test_censored_window_not_completed(self):
        r=log_growth_fit([0,1,2,3],[1,2,4,8],0,10);self.assertEqual(r['status'],'partial_window')
    def test_insufficient(self):self.assertIsNone(log_growth_fit([0,1],[1,2],0,10)['slope_per_s'])
    def test_decaying_not_positive_growth(self):
        r=log_growth_fit([0,1,2,3],np.exp(-np.arange(4)),0,3);self.assertLess(r['slope_per_s'],0);self.assertIsNone(r['e_folding_s'])
    def test_zero_energy_unknown(self):self.assertIsNone(energy_depth([0,.1,.2],[0,0,0],.9))
    def test_surface_energy_zero_depth(self):self.assertEqual(energy_depth([0,.1,.2],[0,0,1],.9),0)
    def test_initial_deep_profile_not_new_penetration(self):
        z=np.linspace(0,.2,201);s=np.exp((z-.2)/.02)
        self.assertAlmostEqual(energy_depth(z,s,.9),energy_depth(z,s*100,.9))
    def test_invalid_grid(self):
        with self.assertRaises(ValueError):energy_depth([0,.2,.1],[1,1,1],.9)
    def test_fundamental_and_harmonic_are_m2_m4(self):
        x=np.linspace(0,1,25);z=np.linspace(0,1,4);d=np.tile(np.cos(2*np.pi*x)+.1*np.cos(4*np.pi*x),(4,1))
        c,e=neumann_spectrum(x,z,d);self.assertAlmostEqual(e[1],.5);self.assertAlmostEqual(e[3],.005);self.assertLess(sum(e[::2]),1e-28)
    def test_nyquist_energy_is_not_doubled(self):
        x=np.linspace(0,1,25);z=np.linspace(0,1,4);d=np.tile(np.cos(24*np.pi*x),(4,1))
        c,e=neumann_spectrum(x,z,d);self.assertAlmostEqual(e[-1],1.)
if __name__=='__main__':unittest.main()
