import unittest
from wp6_timestep_bounds import local_tail_bound,pair_tail_bound,sensitivity_interval

def row(t,c,adv=2.,diff=-1.):
    return dict(time_s=t,cumulative_cut_salt_m=c,cut_salt_adv_down_m_s=adv,cut_salt_diff_down_m_s=diff,cut_salt_total_down_m_s=adv+diff)

class BoundsTests(unittest.TestCase):
    def test_linear_flux_exact_tail_inside_bound(self):
        b=local_tail_bound(row(0,0),row(2,2),1.5)
        self.assertGreaterEqual(b['bound_m'],.5);self.assertAlmostEqual(b['time_flux_bound_m'],1.5)
    def test_no_extrapolation(self):
        for t in [-.1,2.1]:
            with self.assertRaises(ValueError):local_tail_bound(row(0,0),row(2,2),t)
    def test_cancellation_not_assumed(self):
        b=local_tail_bound(row(0,0),row(2,0,adv=1,diff=-1),1)
        self.assertEqual(b['time_flux_bound_m'],2)
    def test_roundoff_or_quadrature_defect_retained(self):
        b=local_tail_bound(row(0,0),row(2,2.01),2)
        self.assertAlmostEqual(b['step_quadrature_defect_m'],.01)
    def test_pair_uses_triangle_bound(self):
        b=pair_tail_bound(row(0,0),row(2,2),row(0,0),row(2,4,adv=3,diff=-1),1)
        self.assertGreaterEqual(b['bound_m'],7)
    def test_interval_contains_exact_time_difference(self):
        r=sensitivity_interval(3,4,.2);self.assertAlmostEqual(r['U_dt_lower_m'],.8);self.assertAlmostEqual(r['U_dt_upper_m'],1.2)
    def test_nonnegative_lower_limit(self):self.assertEqual(sensitivity_interval(1,1,.1)['U_dt_lower_m'],0)
    def test_nonfinite_rejected(self):
        with self.assertRaises(ValueError):sensitivity_interval(1,float('nan'),1)
    def test_unknown_cumulative_not_zero(self):
        with self.assertRaises(ValueError):local_tail_bound(row(0,0),row(2,float('nan')),1)
if __name__=='__main__':unittest.main()
