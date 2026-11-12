import json
import unittest
from pathlib import Path

import pandas as pd


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"


class ValidationResultTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.summary = json.loads(
            (RESULTS / "validation_summary.json").read_text(encoding="utf-8")
        )
        cls.main = pd.read_csv(RESULTS / "main_validation.csv")
        cls.windows = pd.read_csv(RESULTS / "quasi_steady_windows.csv")

    def test_complete_controlled_matrix(self):
        self.assertEqual(len(self.main), 30)
        self.assertEqual(set(self.main["da"]), {1, 4, 10, 20, 45, 100})
        self.assertEqual(set(self.main["u"]), {0, 0.25, 0.5, 0.75, 0.9})
        self.assertEqual(self.summary["ale_completed_cases"], 135)
        statuses = pd.read_csv(RESULTS / "case_status.csv")
        self.assertEqual(len(statuses), 135)
        self.assertTrue((statuses["status"] == "forcing_end").all())

    def test_independent_bvp_passes(self):
        bvp = pd.read_csv(RESULTS / "bvp_verification.csv")
        self.assertEqual(len(bvp), 30)
        self.assertLess(bvp["velocity_relative_error"].max(), 1e-9)
        self.assertLess(bvp["profile_max_relative_error"].max(), 1e-8)

    def test_primary_growth_law_acceptance(self):
        self.assertTrue(self.main["quasi_steady_valid"].all())
        self.assertLess(self.main["ratio_absolute_error"].max(), 0.01)
        self.assertLess(self.main["velocity_relative_error"].max(), 0.01)

    def test_conservation_and_rate_identity(self):
        self.assertLess(self.summary["max_rate_balance_abs_residual_kg_m2_s"], 1e-15)
        self.assertLess(self.summary["max_salt_inventory_residual_kg_m2"], 1e-12)
        self.assertLess(self.summary["max_water_inventory_residual_kg_m2"], 1e-12)

    def test_numerical_convergence(self):
        self.assertLess(
            self.summary["spatial_n64_p1p5_max_relative_difference_from_n128"],
            0.01,
        )
        self.assertLess(self.summary["temporal_max_relative_difference"], 0.01)

    def test_boundary_is_not_silently_dirichlet(self):
        boundary = pd.read_csv(RESULTS / "boundary_audit.csv")
        d1 = boundary.loc[boundary["case_id"] == "boundary_da1_dirichlet"].iloc[0]
        n1 = self.main.loc[(self.main["da"] == 1) & (self.main["u"] == 0)].iloc[0]
        relative_difference = abs(d1["velocity_mean_m_s"] - n1["velocity_mean_m_s"]) / abs(
            n1["velocity_mean_m_s"]
        )
        self.assertGreater(relative_difference, 0.1)

    def test_x008_flux_mismatch_is_retained(self):
        x008 = self.summary["x008"]
        self.assertGreater(abs(x008["reconstructed_minus_soil_fraction"]), 0.01)


if __name__ == "__main__":
    unittest.main()
