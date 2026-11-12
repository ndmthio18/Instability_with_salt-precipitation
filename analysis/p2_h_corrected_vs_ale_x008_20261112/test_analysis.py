import json
import unittest
from pathlib import Path

import numpy as np
import pandas as pd

from analyze_compare import piecewise_linear_integral


HERE = Path(__file__).resolve().parent


class ComparisonAudit(unittest.TestCase):
    def test_piecewise_linear_integral(self):
        t = np.array([0.0, 2.0, 5.0])
        y = np.array([1.0, 3.0, 0.0])
        q = np.array([0.0, 1.0, 2.0, 3.0, 5.0])
        expected = np.array([0.0, 1.5, 4.0, 6.5, 8.5])
        np.testing.assert_allclose(piecewise_linear_integral(t, y, q), expected)

    def test_full_support_and_conservation(self):
        results = HERE / "results"
        ale = pd.read_csv(results / "ale_history.csv")
        forcing = pd.read_csv(results / "x008_forcing.csv")
        self.assertAlmostEqual(ale.time_s.iloc[-1], forcing.elapsed_time_s.iloc[-1], places=9)
        self.assertLess(ale.salt_residual_kg_m2.abs().max(), 1e-10)
        self.assertLess(ale.water_residual_kg_m2.abs().max(), 1e-10)
        self.assertTrue((ale.epsilon_bottom <= 1.0).all())
        self.assertAlmostEqual(ale.epsilon_bottom.iloc[-1], 1.0, places=12)

    def test_common_observable_agreement(self):
        summary = json.loads((HERE / "results" / "comparison_summary.json").read_text())
        self.assertEqual(summary["status"], "completed_full_forcing_support")
        self.assertLess(abs(summary["final_relative_difference_of_increment"]), 1e-4)
        self.assertGreater(summary["ale_final_geometric_top_m"], summary["ale_final_net_height_m"])


if __name__ == "__main__":
    unittest.main()
