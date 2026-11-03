import unittest
import numpy as np
from analyze_wp8 import pairs,conservative
from wp8_constraint import EqualHeightConstraint
from wp6_4t_metrics import neumann_spectrum
from single_wavelength_post import widths

class TestWP8(unittest.TestCase):
    def test_neumann_energy_and_physical_k(self):
        for nx in [24,48]:
            x=np.linspace(0,.6/28,nx+1);z=np.linspace(0,.2,9);d=np.tile(3*np.cos(6*np.pi*x/x[-1]),(9,1));c,e=neumann_spectrum(x,z,d)
            self.assertAlmostEqual(e.sum(),.9);self.assertEqual(e.argmax()+1,6);self.assertAlmostEqual((np.arange(1,nx+1)*np.pi/x[-1])@e/e.sum(),6*np.pi/x[-1])
    def test_covariance_cannot_be_dropped(self):
        x=np.array([1.,2.,3.]);s=np.array([.1,.2,.3]);w=np.array([.25,.5,.25]);cov=((s-s@w)*(x-x@w))@w
        self.assertAlmostEqual((s*x)@w,(s@w)*(x@w)+cov);self.assertGreater(cov,0)
    def test_exact_time_no_interpolation(self):self.assertEqual(pairs([0.,1.,2.],[0.,1.00001,2.]),[(0,0),(2,2)])
    def test_initial_correction_and_identity(self):
        def r(t,I,C):return dict(time_s=t,actual_z_over_H=.9,inventory_salt_m=I,cumulative_cut_salt_m=C,salt_residual_m=0,cut_salt_adv_down_m_s=0,cut_salt_diff_down_m_s=0,cut_salt_total_down_m_s=0)
        rr=conservative([r(0,1,0),r(1,1,0)],[r(0,1.25,0),r(1,1.,.25)])
        self.assertEqual(rr[-1]['delta_I_corrected_m'],-.25);self.assertEqual(rr[-1]['identity_error_m'],0)
    def test_constraint_exact_uniform(self):
        c=EqualHeightConstraint([0,0,1,1]);u=np.arange(4.).reshape(2,2);v=c.lift(u);np.testing.assert_array_equal(v[::2],v[1::2]);np.testing.assert_array_equal(c.project_initial(v),u)
    def test_reduced_jacobian_and_conservation(self):
        c=EqualHeightConstraint([0,0,1,1]);J=np.arange(64.).reshape(8,8);u=np.arange(4.);f=lambda v:J@v.ravel();np.testing.assert_allclose(c.residual(u,f),c.jacobian(J)@u);self.assertEqual(c.residual(u,f).sum(),f(c.lift(u)).sum())
    def test_reject_posthoc_projection(self):
        c=EqualHeightConstraint([0,0]);
        with self.assertRaises(ValueError):c.project_initial([[1,2],[1,3]])
    def test_constraint_nonlinear_derivative(self):
        c=EqualHeightConstraint([0,0,1,1]);u=np.array([1.,2.,3.,4.]);f=lambda v:v.ravel()**2;J=np.diag(2*c.lift(u).ravel());h=1e-5
        fd=np.column_stack([(c.residual(u+np.eye(4)[j]*h,f)-c.residual(u-np.eye(4)[j]*h,f))/(2*h) for j in range(4)])
        np.testing.assert_allclose(fd,c.jacobian(J),rtol=1e-10)
    def test_conservative_nonlinear_step(self):
        # Integrated two-column nonlinear storage with cancelling internal flux.
        # Unequal lateral cell volumes still conserve total storage after reduction.
        c=EqualHeightConstraint([0,0,1,1],components=1);old=np.array([[.2],[.2],[.7],[.7]]);mass=np.array([.25,.75,.25,.75]);dt=.1
        L=np.array([[1,0,-1,0],[0,1,0,-1],[-1,0,1,0],[0,-1,0,1.]])
        def assembly(v):
            v=v.ravel();r=mass*(v+v*v-old.ravel()-old.ravel()**2)+dt*L@v;j=np.diag(mass*(1+2*v))+dt*L;return r,j
        v,info=c.solve(old,assembly,atol=1e-14);np.testing.assert_array_equal(v[::2],v[1::2]);self.assertLess(info['reduced_residual_Linf'],1e-14);self.assertAlmostEqual(mass@(v.ravel()+v.ravel()**2),mass@(old.ravel()+old.ravel()**2),places=13)
    def test_dirichlet_and_failed_newton(self):
        c=EqualHeightConstraint([0,0,1,1],components=1);a=lambda v:(v.ravel()**2-4,np.diag(2*v.ravel()))
        v,_=c.solve(np.ones((4,1)),a,fixed={0:3.});np.testing.assert_array_equal(v[:2],[[3.],[3.]]);np.testing.assert_allclose(v[2:],2.)
        with self.assertRaises(RuntimeError):c.solve(np.ones((4,1)),a,maxiter=0)
if __name__=='__main__':unittest.main()
