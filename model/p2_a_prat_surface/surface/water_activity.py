"""NaCl water activity, reused from the heterogeneity project.

Source: `water_activity.hh::naclWaterActivity` of
`25 Heterogeneity/github_submission/heterogeniety/deliveries/
RH0_fixed_174h_20261103/model_revisions/20261102_cylinder_fixed_hydraulics_v12/
source/appl/water_activity.hh`.

That function implements the pure-NaCl Pitzer osmotic coefficient of
Rose et al. (2008), Appendix A, Eqs. A15-A16, with the Pitzer and Mayorga
(1973) parameters at 298.15 K, and returns

    a_w = exp(-2 * m * M_H2O * phi(m)).

The upstream header states explicitly that use at 293.15 K is a *stated
near-room-temperature approximation, not a fit*. That caveat is inherited here
unchanged: this module does not re-fit, re-calibrate or extend the correlation.

Deliberately NOT ported: `normalizedEvaporationWeights`. Task P2-A item 4
forbids normalising local evaporation weights to a fixed total rate, so the
normalisation helper is excluded on purpose rather than by oversight.
"""

from __future__ import annotations

import math

from .thermo import (
    M_H2O,
    molality_from_mass_fraction,
    molality_from_mole_fraction,
)

#: Upstream trial-state molality limit (`water_activity.hh`, v09 comment).
MOLALITY_LIMIT = 8.0

# Pitzer coefficients exactly as transcribed from `water_activity.hh`.
_A_PHI = 0.3915
_B_PITZER = 1.2
_BETA0 = 0.0765
_BETA1 = 0.2664
_ALPHA = 2.0
_C_PHI = 0.00127


def osmotic_coefficient(molality: float) -> float:
    """Pitzer osmotic coefficient phi(m) for aqueous NaCl."""
    _check_molality(molality)
    r = math.sqrt(molality)
    return (
        1.0
        - _A_PHI * r / (1.0 + _B_PITZER * r)
        + molality * (_BETA0 + _BETA1 * math.exp(-_ALPHA * r))
        + molality * molality * _C_PHI
    )


def water_activity(molality: float, water_molar_mass: float = M_H2O) -> float:
    """a_w for aqueous NaCl at the given molality [mol/kg H2O]."""
    if not math.isfinite(water_molar_mass) or water_molar_mass <= 0.0:
        raise ValueError("water activity requires a positive water molar mass")
    phi = osmotic_coefficient(molality)
    return math.exp(-2.0 * molality * water_molar_mass * phi)


def water_activity_from_mole_fraction(x_salt: float) -> float:
    """a_w from the liquid NaCl mole fraction (the Richards primary variable)."""
    return water_activity(molality_from_mole_fraction(x_salt))


def water_activity_from_mass_fraction(w_salt: float) -> float:
    """a_w from the liquid NaCl mass fraction (Prat's `C`)."""
    return water_activity(molality_from_mass_fraction(w_salt))


def _check_molality(molality: float) -> None:
    if not math.isfinite(molality) or molality < 0.0 or molality > MOLALITY_LIMIT:
        raise ValueError(
            "NaCl water activity requires trial molality in [0, "
            f"{MOLALITY_LIMIT}]; got {molality!r}. Accepted mineral equilibrium "
            "is checked separately and this range is not a validation of "
            "supersaturated states."
        )
