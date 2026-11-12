from __future__ import annotations

import math
from pathlib import Path
import sys
import unittest

import numpy as np

REPO = Path(__file__).resolve().parents[2]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from model.p2_h_dual_fronts.analytical import (  # noqa: E402
    finite_difference_profile,
    modified_velocity_factor,
    prat_equivalent_u,
    prat_velocity_factor,
    velocity_ratio,
)
from model.p2_h_dual_fronts.dual_front import BoundaryKinematics, DualFrontModel  # noqa: E402


class AnalyticalTests(unittest.TestCase):
    def test_velocity_ratio_identity(self):
        for da in (0.1, 1.0, 20.0, 100.0):
            for x in (0.5, 2.0, 5.0):
                for u in (0.0, 0.5, 1.25):
                    ratio = modified_velocity_factor(da, x, u) / modified_velocity_factor(da, x, 0.0)
                    self.assertAlmostEqual(ratio, velocity_ratio(da, u), places=14)

    def test_independent_finite_difference_solution(self):
        xi, y = finite_difference_profile(20.0, 2.0, 0.75, cells=600)
        self.assertEqual(len(xi), 601)
        self.assertLess(abs(2.0 * y[-1] - modified_velocity_factor(20.0, 2.0, 0.75)), 2.0e-5)

    def test_prat_requires_state_dependent_flux(self):
        for da in (0.1, 1.0, 10.0, 45.0):
            for x in (0.5, 2.0, 5.0):
                u = prat_equivalent_u(da, x)
                self.assertAlmostEqual(
                    modified_velocity_factor(da, x, u),
                    prat_velocity_factor(da, x),
                    places=13,
                )


class ConservationTests(unittest.TestCase):
    def test_translation_conserves_total_salt_and_water(self):
        model = DualFrontModel(40, top_porosity=0.43, bottom_porosity=0.50)
        initial = model.initial_state(1.0e-4, 0.26494, porosity=0.465)
        final, _ = model.run(initial, BoundaryKinematics(1e-9, 1e-9, 0.4), 100.0, 0.5)
        before = model.inventory(initial)
        after = model.inventory(final)
        self.assertAlmostEqual(before["height_m"], after["height_m"], places=15)
        self.assertLess(abs(before["total_salt_with_external_kg_m2"] - after["total_salt_with_external_kg_m2"]), 1e-12)
        self.assertLess(abs(before["total_water_with_external_kg_m2"] - after["total_water_with_external_kg_m2"]), 1e-12)
        self.assertGreater(abs(before["solid_kg_m2"] - after["solid_kg_m2"]), 1e-8)

    def test_thinning_stops_at_safety_limit(self):
        model = DualFrontModel(20, minimum_height_m=9.5e-5)
        initial = model.initial_state(1.0e-4, 0.26494)
        with self.assertRaises(RuntimeError):
            model.run(initial, BoundaryKinematics(0.0, 1e-8), 1000.0, 1.0)


if __name__ == "__main__":
    unittest.main()
