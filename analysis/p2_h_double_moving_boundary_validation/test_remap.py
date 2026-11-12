import unittest
from pathlib import Path
import numpy as np
from remap import overlap_remap,convert_connected_zero_prefix

class RemapTests(unittest.TestCase):
    def test_connected_zero_transfer(self):
        v=np.array([.1,.2,.3,.4]);s=np.array([0,0,.1,0]);q=np.array([.2,.3,.4,.5])
        a,b,c,d,e,n=convert_connected_zero_prefix([.2],[.1],v,s,q)
        self.assertEqual(n,2);self.assertEqual(len(d),2) # isolated interior zero stays
        self.assertAlmostEqual(a.sum()+c.sum(),1.2)
        self.assertAlmostEqual(b.sum()+e.sum(),1.5)
        self.assertAlmostEqual(d.sum(),s.sum())
        self.assertAlmostEqual((a.sum()+c.sum()-d.sum())*1200-b.sum()-e.sum(),
                               (1.2-s.sum())*1200-1.5)

    def test_no_threshold_deletion(self):
        result=convert_connected_zero_prefix([.1],[.01],[.2,.3],[1e-30,0],[.03,.04])
        self.assertEqual(result[-1],0)

    def test_overlap_mass_and_convergence(self):
        errors=[]
        for n in [32,64,128]:
            a=np.linspace(0,1,n+1);b=np.linspace(0,1,n+4)**1.3
            q=np.diff(a*a/2+a) # exact cell inventory of 1+z
            mapped=overlap_remap(a,q,b)
            self.assertAlmostEqual(mapped.sum(),q.sum(),places=14)
            self.assertTrue(np.all(mapped>=0))
            errors.append(np.sum(abs(mapped-np.diff(b*b/2+b))))
        self.assertLess(errors[1],errors[0]);self.assertLess(errors[2],errors[1])

    def test_kinematics_not_analytical_growth_law(self):
        source=(Path(__file__).parent/'two_front.hh').read_text()
        for token in ['phi_core','DaStar','sech','tanh','v_corr']:
            self.assertNotIn(token,source)

if __name__=='__main__':unittest.main()
