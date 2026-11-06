#!/usr/bin/env python3
"""Pure numerical unit tests for the P1-R2 conversion and diffusion closure."""

from __future__ import annotations

import math
import unittest


M_WATER = 0.01801528
M_SALT = 0.05844
RHO_WATER = 998.203
LINEAR = 0.670804
QUADRATIC = 0.373854


def mixture_molar_mass(x_salt: float) -> float:
    return (1.0 - x_salt) * M_WATER + x_salt * M_SALT


def mass_fraction(x_salt: float) -> float:
    return x_salt * M_SALT / mixture_molar_mass(x_salt)


def mass_fraction_derivative(x_salt: float) -> float:
    return M_SALT * M_WATER / mixture_molar_mass(x_salt) ** 2


def density(x_salt: float) -> float:
    w_salt = mass_fraction(x_salt)
    return RHO_WATER + 1000.0 * w_salt * (LINEAR + QUADRATIC * w_salt)


class PhysicsConversionTest(unittest.TestCase):
    def test_mass_fraction_derivative(self) -> None:
        for x_salt in (0.0, 0.0036, 0.02, 0.08, 0.1):
            h = 1e-7
            numerical = (
                mass_fraction(x_salt + h) - mass_fraction(x_salt - h)
            ) / (2.0 * h)
            analytical = mass_fraction_derivative(x_salt)
            relative = abs(numerical - analytical) / analytical
            self.assertLess(relative, 1e-9)

    def test_component_mass_mole_identity(self) -> None:
        theta = 0.37
        for x_salt in (0.0, 0.0036, 0.02, 0.08, 0.1):
            mbar = mixture_molar_mass(x_salt)
            rho = density(x_salt)
            c = rho / mbar
            salt_mass = theta * rho * mass_fraction(x_salt)
            water_mass = theta * rho * (1.0 - mass_fraction(x_salt))
            self.assertTrue(
                math.isclose(
                    salt_mass,
                    M_SALT * theta * c * x_salt,
                    rel_tol=5e-13,
                    abs_tol=1e-14,
                )
            )
            self.assertTrue(
                math.isclose(
                    water_mass,
                    M_WATER * theta * c * (1.0 - x_salt),
                    rel_tol=5e-13,
                    abs_tol=1e-14,
                )
            )

    def test_normalized_conservation_is_not_physical_conservation(self) -> None:
        theta = 0.3
        before = (0.02, 0.08)
        after = (0.05, 0.05)
        self.assertAlmostEqual(sum(before), sum(after), places=15)

        def physical(values: tuple[float, float]) -> tuple[float, float]:
            total = sum(theta * density(value) for value in values)
            salt = sum(
                theta * density(value) * mass_fraction(value)
                for value in values
            )
            return total, salt

        total_before, salt_before = physical(before)
        total_after, salt_after = physical(after)
        self.assertGreater(
            abs(total_after - total_before) / total_before,
            1e-4,
        )
        self.assertGreater(
            abs(salt_after - salt_before) / salt_before,
            1e-4,
        )


if __name__ == "__main__":
    unittest.main()
