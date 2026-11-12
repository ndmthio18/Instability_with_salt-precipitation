import csv
import unittest
from pathlib import Path


HERE = Path(__file__).resolve().parent


class ConfigurationTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        with (HERE / "configs" / "cases.csv").open(newline="", encoding="utf-8") as handle:
            cls.rows = list(csv.DictReader(handle))

    def test_case_count_and_ids(self):
        self.assertEqual(len(self.rows), 135)
        ids = [row["case_id"] for row in self.rows]
        self.assertEqual(len(ids), len(set(ids)))

    def test_main_matrix_complete(self):
        main = [row for row in self.rows if row["suite"] == "main"]
        self.assertEqual(len(main), 30)
        self.assertEqual({float(row["da"]) for row in main}, {1, 4, 10, 20, 45, 100})
        self.assertEqual({float(row["u"]) for row in main}, {0, 0.25, 0.5, 0.75, 0.9})
        self.assertEqual({row["boundary_mode"] for row in main}, {"neumann"})

    def test_da_star_not_varied(self):
        forbidden = {"kr", "av", "d0", "eps", "beta", "grain"}
        self.assertTrue(forbidden.isdisjoint(self.rows[0]))

    def test_convergence_matrix(self):
        spatial = [row for row in self.rows if row["suite"] == "spatial"]
        temporal = [row for row in self.rows if row["suite"] == "temporal"]
        self.assertEqual({int(row["cells"]) for row in spatial}, {32, 64, 128})
        self.assertEqual({float(row["power"]) for row in spatial}, {1.0, 1.5, 2.0, 3.0})
        self.assertEqual(len(temporal), 24)


if __name__ == "__main__":
    unittest.main()
