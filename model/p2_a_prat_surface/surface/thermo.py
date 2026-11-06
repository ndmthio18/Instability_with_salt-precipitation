"""Thermophysical correlations used by the P2-A surface subgrid.

Provenance of every correlation in this file (nothing is written from memory):

* Saturation pressure of pure water: IAPWS-IF97 Region 4 backward equation as
  implemented in DuMuX `dumux/material/components/iapws/region4.hh`
  (`Region4::saturationPressure`). The ten coefficients and the evaluation
  sequence are transcribed from that file.
* Binary gas diffusion coefficient of H2O in air: DuMuX
  `dumux/material/binarycoefficients/h2o_air.hh` (`H2O_Air::gasDiffCoeff`),
  i.e. ``D = 2.13e-5 * (1e5/p) * (T/273.15)**1.8`` after Vargaftik (1975) and
  Walker et al. (1981) as cited there.
* Universal gas constant: DuMuX `dumux/material/constants.hh`
  (``R = 8.314472``).
* Air molar mass: DuMuX `dumux/material/components/air.hh`
  (``molarMass() = 0.02896``).
* H2O and NaCl molar masses: the frozen Case 36877 input file
  (`Problem.MolarMassH2O = 0.01801528`, `Problem.MolarMassNaCl = 0.05844`).

The gas phase is treated as an **ideal binary mixture of water vapour and air**.
That is an explicit modelling assumption of this standalone module, taken so
that the film state can be evaluated without a DuMuX fluid system. It is NOT a
claim that the heterogeneity project's `FS::density` is an ideal-gas call.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

# --- source-traced constants -------------------------------------------------

R_GAS = 8.314472           # J mol^-1 K^-1   (DuMuX Constants::R)
M_H2O = 0.01801528         # kg mol^-1       (Case 36877 input)
M_NACL = 0.05844           # kg mol^-1       (Case 36877 input)
M_AIR = 0.02896            # kg mol^-1       (DuMuX Air::molarMass)

# IAPWS-IF97 Region 4 coefficients (DuMuX region4.hh)
_IAPWS_N = (
    0.11670521452767e4, -0.72421316703206e6, -0.17073846940092e2,
    0.12020824702470e5, -0.32325550322333e7, 0.14915108613530e2,
    -0.48232657361591e4, 0.40511340542057e6, -0.23855557567849,
    0.65017534844798e3,
)


def saturation_pressure(temperature: float) -> float:
    """Pure-water saturation pressure [Pa]; IAPWS-IF97 Region 4."""
    if not math.isfinite(temperature) or temperature <= 0.0:
        raise ValueError("saturation_pressure requires a positive finite temperature")
    n = _IAPWS_N
    sigma = temperature + n[8] / (temperature - n[9])
    a = (sigma + n[0]) * sigma + n[1]
    b = (n[2] * sigma + n[3]) * sigma + n[4]
    c = (n[5] * sigma + n[6]) * sigma + n[7]
    tmp = 2.0 * c / (math.sqrt(b * b - 4.0 * a * c) - b)
    tmp *= tmp
    tmp *= tmp
    return 1.0e6 * tmp


def vapour_diffusion_coefficient(temperature: float, pressure: float) -> float:
    """Binary H2O-air gas diffusion coefficient [m^2/s] (DuMuX H2O_Air)."""
    if not math.isfinite(temperature) or temperature <= 0.0:
        raise ValueError("vapour diffusion requires a positive finite temperature")
    if not math.isfinite(pressure) or pressure <= 0.0:
        raise ValueError("vapour diffusion requires a positive finite pressure")
    return 2.13e-5 * (1.0e5 / pressure) * (temperature / 273.15) ** 1.8


def gas_water_mass_fraction(y_water: float) -> float:
    """Mass fraction of water in an ideal H2O/air gas at mole fraction `y_water`.

    Mirrors `gasWaterMassFraction_` in the heterogeneity `problem.hh`.
    """
    if not math.isfinite(y_water) or y_water < 0.0 or y_water > 1.0:
        raise ValueError("gas water mole fraction must lie in [0, 1]")
    return M_H2O * y_water / (M_H2O * y_water + M_AIR * (1.0 - y_water))


def ideal_gas_density(temperature: float, pressure: float, y_water: float) -> float:
    """Ideal binary-gas mass density [kg/m^3] at the given vapour mole fraction."""
    molar_mass = y_water * M_H2O + (1.0 - y_water) * M_AIR
    return pressure * molar_mass / (R_GAS * temperature)


# --- liquid composition conversions (binary H2O-NaCl) ------------------------

def mass_fraction_from_mole_fraction(x_salt: float) -> float:
    """w_s from x_s, identical to the P1-R2 `massFraction` closure."""
    mixture_molar_mass = (1.0 - x_salt) * M_H2O + x_salt * M_NACL
    return x_salt * M_NACL / mixture_molar_mass


def mole_fraction_from_mass_fraction(w_salt: float) -> float:
    """Inverse of `mass_fraction_from_mole_fraction`."""
    if w_salt < 0.0 or w_salt >= 1.0:
        raise ValueError("salt mass fraction must lie in [0, 1)")
    n_salt = w_salt / M_NACL
    n_water = (1.0 - w_salt) / M_H2O
    return n_salt / (n_salt + n_water)


def molality_from_mole_fraction(x_salt: float) -> float:
    """m [mol NaCl / kg H2O]; mirrors the heterogeneity top-face expression
    ``m = x_salt / (x_water * M_H2O)``."""
    if not math.isfinite(x_salt) or x_salt < 0.0 or x_salt >= 1.0:
        raise ValueError("salt mole fraction must lie in [0, 1)")
    return x_salt / ((1.0 - x_salt) * M_H2O)


def molality_from_mass_fraction(w_salt: float) -> float:
    """m [mol NaCl / kg H2O] from the salt mass fraction."""
    if not math.isfinite(w_salt) or w_salt < 0.0 or w_salt >= 1.0:
        raise ValueError("salt mass fraction must lie in [0, 1)")
    return w_salt / (M_NACL * (1.0 - w_salt))


@dataclass(frozen=True)
class BrineDensityLaw:
    """The frozen P1 Richards-salt liquid density law, kept bit-for-bit.

    ``rho_l = rho_w0 + 1000 * w_s * (a + b * w_s)`` with the Case 36877
    coefficients. Reproduced from `main.cc::rhoW` of the P1-R2 isolated source.
    """

    pure_water_density: float = 998.203
    linear: float = 0.670804
    quadratic: float = 0.373854

    def density(self, w_salt: float) -> float:
        return self.pure_water_density + 1000.0 * w_salt * (
            self.linear + self.quadratic * w_salt
        )

    def density_from_mole_fraction(self, x_salt: float) -> float:
        return self.density(mass_fraction_from_mole_fraction(x_salt))
