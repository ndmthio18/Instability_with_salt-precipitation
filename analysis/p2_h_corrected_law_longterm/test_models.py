import unittest
from pathlib import Path
import numpy as np
import pandas as pd
from analyze import frozen_bvp,law,birth,ceq,grain,ds,kr,rho,cr

class ModelChecks(unittest.TestCase):
    def test_frozen_bvp_constant_property_limit(self):
        # Independent BVP must tend to the verified closed-form law.
        for u in [0,.9]:
            h=3.5e-6;j=1.78e-5;B=u*ceq*j;n=1024
            g=pd.DataFrame({'eps':np.full(n,birth),'h':h,'width':np.full(n,1/n),'xi':(np.arange(n)+.5)/n})
            predicted=frozen_bvp(g,j,B);exact=law(h,j,B,birth)[0]
            self.assertLess(abs(predicted/exact-1),1e-5)

    def test_screening_limit(self):
        a=law(1e-3,1.78e-5,0,birth)[0]
        b=law(1e-3,1.78e-5,.9*ceq*1.78e-5,birth)[0]
        self.assertAlmostEqual(a,b,places=18)

    def test_material_definitions(self):
        result=law(3.5e-6,1.78e-5,0,birth)
        self.assertAlmostEqual(result[2],3.1750355738459297,places=12)

class ProductionEvidenceChecks(unittest.TestCase):
    root=Path(__file__).resolve().parent

    def test_all_registered_runs_accounted_for(self):
        tables=[pd.read_csv(self.root/f'cases{suffix}.csv') for suffix in ['', '_p1','_events','_refine']]
        expected=pd.concat(tables).id
        summary=pd.read_csv(self.root/'results/summary.csv')
        self.assertEqual(set(expected),set(summary.id))
        self.assertFalse(summary.id.duplicated().any())
        self.assertEqual(len(pd.read_csv(self.root/'results/pending_cases.csv')),0)
        # Original failed trajectories must remain in the record.
        self.assertEqual(summary.status.str.startswith('solver_stop').sum(),8)

    def test_reconstructed_inventory_budgets(self):
        h=pd.read_csv(self.root/'results/history.csv.gz')
        for name,g in h.groupby('id'):
            a=g.iloc[0]
            residuals=[g.solid-a.solid-g.Qtop-g.Qplus+g.Qminus,
                g.dissolved-a.dissolved-g.Qb+g.Qtop+g.Qplus-g.Qminus,
                g.water-a.water-g.Qliquid+g.Qb+g.Qevap]
            for residual in residuals:self.assertLess(residual.abs().max(),2e-9,name)

    def test_prescribed_diffusion_not_equilibrium_advection(self):
        h=pd.read_csv(self.root/'results/history.csv.gz')
        np.testing.assert_allclose(h.diffb,-h.B,atol=1e-13,rtol=1e-10)
        np.testing.assert_allclose(h.Jb,h.advb+h.diffb,atol=1e-13,rtol=1e-10)

    def test_thresholds_and_numerical_gates(self):
        s=pd.read_csv(self.root/'results/summary.csv')
        self.assertTrue((s.loc[s.status.str.startswith('solid_framework'),'final_epsmax']>=.999).all())
        self.assertTrue((s.loc[s.status.str.startswith('pore_filling'),'final_epsmin']<=.05).all())
        c=pd.read_csv(self.root/'results/convergence.csv')
        err=c.loc[c.id=='u0_n256_p2_r1e-08_dt30','velocity_max_relative'].iloc[0]
        self.assertLess(err,.01)
        t=pd.read_csv(self.root/'results/temporal_convergence.csv')
        self.assertLess(t.max_velocity_relative.max(),.01)

if __name__=='__main__':unittest.main()
