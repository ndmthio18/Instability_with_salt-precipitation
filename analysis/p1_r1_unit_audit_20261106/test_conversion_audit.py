#!/usr/bin/env python3
"""Regression tests for exact Case 36877 conversion inputs and identities."""

from __future__ import annotations

import tempfile
import unittest
from pathlib import Path

import conversion_audit as audit


class ConversionAuditTest(unittest.TestCase):
    def test_exact_frozen_parameters_are_loaded(self) -> None:
        self.assertEqual(audit.CASE_PARAMETERS, audit.EXPECTED_CASE_PARAMETERS)
        state = audit.properties(0.0036)
        self.assertAlmostEqual(state["Mbar_kg_per_mol"], 0.018160808992, places=15)
        self.assertAlmostEqual(state["w_NaCl"], 0.011584506, places=9)
        self.assertAlmostEqual(state["rho_kg_per_m3_liquid"], 1006.024104, places=6)
        self.assertAlmostEqual(state["c_mol_per_m3_liquid"], 55395.3354, places=4)

    def test_rounded_constants_are_rejected(self) -> None:
        rounded = """
[Problem]
MolarMassH2O = 0.018
MolarMassNaCl = 0.0584428
PureWaterDensity = 1000
DensityLinearCoefficient = 0.6708
DensityQuadraticCoefficient = 0.3739
EvaporationRate = 3.065544359e-9
SaltMoleFraction = 0.0036
""".strip()
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "rounded.input"
            path.write_text(rounded + "\n", encoding="utf-8")
            with self.assertRaisesRegex(ValueError, "inconsistent Case 36877"):
                audit.read_case_parameters(path)

    def test_mass_molar_identities(self) -> None:
        for x_nacl in (0.0, 0.0036, 0.008018127628, 0.1):
            state = audit.properties(x_nacl)
            salt_from_moles = (
                state["c_mol_per_m3_liquid"] * x_nacl * audit.M_NACL
            )
            water_from_moles = (
                state["c_mol_per_m3_liquid"]
                * (1.0 - x_nacl)
                * audit.M_H2O
            )
            self.assertAlmostEqual(
                state["salt_mass_concentration_kg_per_m3_liquid"],
                salt_from_moles,
                places=11,
            )
            self.assertAlmostEqual(
                state["water_mass_concentration_kg_per_m3_liquid"],
                water_from_moles,
                places=11,
            )


if __name__ == "__main__":
    unittest.main()
