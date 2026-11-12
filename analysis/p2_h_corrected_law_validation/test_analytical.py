import unittest

import numpy as np

from analytical_model import (
    Parameters,
    analytical_profile,
    corrected_observables,
    solve_independent_bvp,
    thickness_for_da,
)


class AnalyticalLawTests(unittest.TestCase):
    def setUp(self):
        self.p = Parameters()

    def test_dimensionless_identity(self):
        expected = (
            self.p.beta**2
            * (1 - self.p.eps) ** 2
            * (1 - self.p.ceq)
            * self.p.kr
            / (self.p.eps * self.p.deff * self.p.av)
        )
        self.assertAlmostEqual(self.p.da_star, expected, places=14)

    def test_target_da_thickness(self):
        for da in (1, 4, 10, 20, 45, 100):
            h = thickness_for_da(da, self.p)
            self.assertAlmostEqual((h / self.p.ell) ** 2, da, places=12)

    def test_boundary_conditions(self):
        for da, u in ((1.0, 0.9), (20.0, 0.5), (100.0, 0.0)):
            h = thickness_for_da(da, self.p)
            c, dcdz = analytical_profile(np.array([0.0, h]), da, u, self.p)
            b = u * self.p.ceq * self.p.j
            self.assertAlmostEqual(
                self.p.rho * self.p.eps * self.p.deff * dcdz[0], b, places=15
            )
            top_residual = (
                self.p.s_ref
                - self.p.rho * self.p.mathcal_d * dcdz[-1]
                - self.p.rho * self.p.ks * c[-1]
            )
            self.assertAlmostEqual(top_residual, 0.0, places=15)

    def test_normalized_screening_identity(self):
        for da in (1.0, 4.0, 20.0, 100.0):
            v0 = corrected_observables(da, 0.0, self.p)["velocity_m_s"]
            for u in (0.25, 0.5, 0.9):
                ratio = corrected_observables(da, u, self.p)["velocity_m_s"] / v0
                self.assertAlmostEqual(
                    ratio, 1.0 - u / np.cosh(np.sqrt(da)), places=14
                )

    def test_independent_bvp(self):
        exact = corrected_observables(10.0, 0.75, self.p)
        bvp = solve_independent_bvp(10.0, 0.75, self.p)
        self.assertLess(
            abs(exact["velocity_m_s"] - bvp["velocity_m_s"])
            / exact["velocity_m_s"],
            1e-9,
        )


if __name__ == "__main__":
    unittest.main()
