"""Exact and finite-difference solutions for the prescribed-base-flux BVP.

The coordinate ``xi`` points upward from the crust bottom to its top.  The
positive quantity ``B_soil`` is a *downward* diffusive salt flux, hence

    rho eps D d(C-Ceq)/dz |bottom = B_soil.

The dimensionless concentration returned here is scaled with
``Ceq*j/A_d``, where ``A_d=rho*eps*D*sqrt(Da)/h``.  In this scaling the lower
condition is ``y'(0)=u*sqrt(Da)`` and the upper Robin condition is
``y'(1)+sqrt(Da)*X*y(1)=sqrt(Da)``.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math
from typing import Iterable

import numpy as np


def mass_fraction_from_mole_fraction(x: float) -> float:
    m_water = 0.01801528
    m_nacl = 0.05844
    return x * m_nacl / ((1.0 - x) * m_water + x * m_nacl)


def brine_density(c: float) -> float:
    return 998.203 + 1000.0 * c * (0.670804 + 0.373854 * c)


def sech(value: float) -> float:
    a = abs(value)
    if a > 40.0:
        return 2.0 * math.exp(-a)
    return 1.0 / math.cosh(value)


@dataclass(frozen=True)
class AnalyticalParameters:
    """Actual P2-H parameters used for the 3.5 micrometre audit."""

    height_m: float = 3.5e-6
    porosity: float = 0.43
    crystal_density_kg_m3: float = 2165.0
    reaction_coefficient_m_s: float = 1.0e-3
    beta: float = 2.0
    molecular_diffusivity_m2_s: float = 1.3e-9
    grain_size_m: float = 25.0e-6
    av_multiplier: float = 1.0
    equilibrium_mole_fraction: float = 0.1
    evaporation_kg_m2_s: float = 1.7815e-5

    @property
    def equilibrium_mass_fraction(self) -> float:
        return mass_fraction_from_mole_fraction(self.equilibrium_mole_fraction)

    @property
    def solution_density_kg_m3(self) -> float:
        return brine_density(self.equilibrium_mass_fraction)

    @property
    def effective_diffusivity_m2_s(self) -> float:
        # Exact P2-H closure in ale.hh: d=sqrt(epsilon)*p.diffusion.
        return math.sqrt(self.porosity) * self.molecular_diffusivity_m2_s

    @property
    def specific_area_m_inv(self) -> float:
        return (
            self.av_multiplier
            * 36.0
            * (1.0 - self.porosity)
            / self.grain_size_m
        )

    @property
    def damkohler(self) -> float:
        c = self.equilibrium_mass_fraction
        return (
            (1.0 - c)
            * self.specific_area_m_inv
            * self.reaction_coefficient_m_s
            * self.height_m**2
            / (self.porosity * self.effective_diffusivity_m2_s)
        )

    @property
    def damkohler_star(self) -> float:
        c = self.equilibrium_mass_fraction
        return (
            self.beta**2
            * (1.0 - self.porosity) ** 2
            * (1.0 - c)
            * self.reaction_coefficient_m_s
            / (
                self.porosity
                * self.effective_diffusivity_m2_s
                * self.specific_area_m_inv
            )
        )

    @property
    def penetration_length_m(self) -> float:
        return self.height_m / math.sqrt(self.damkohler)

    def to_record(self) -> dict[str, float]:
        record = asdict(self)
        record.update(
            equilibrium_mass_fraction=self.equilibrium_mass_fraction,
            solution_density_kg_m3=self.solution_density_kg_m3,
            effective_diffusivity_m2_s=self.effective_diffusivity_m2_s,
            specific_area_m_inv=self.specific_area_m_inv,
            damkohler=self.damkohler,
            damkohler_star=self.damkohler_star,
            sqrt_damkohler_star=math.sqrt(self.damkohler_star),
            penetration_length_m=self.penetration_length_m,
        )
        return record


def normalized_top_concentration(da: float, x: float, u: float) -> float:
    """Return ``(Ctop-Ceq)/(Ceq*j/A_d)`` for the flux BVP."""

    lam = math.sqrt(da)
    return (1.0 - u * sech(lam)) / (x + math.tanh(lam))


def modified_velocity_factor(da: float, x: float, u: float) -> float:
    """Factor multiplying the ideal evaporation-limited velocity prefactor."""

    lam = math.sqrt(da)
    return x * (1.0 - u * sech(lam)) / (x + math.tanh(lam))


def prat_velocity_factor(da: float, x: float) -> float:
    lam = math.sqrt(da)
    tanh = math.tanh(lam)
    return x * tanh / (1.0 + x * tanh)


def velocity_ratio(da: float, u: float) -> float:
    return 1.0 - u * sech(math.sqrt(da))


def prat_equivalent_u(da: float, x: float) -> float:
    """State-dependent base flux that makes the flux BVP equal Prat's BVP."""

    lam = math.sqrt(da)
    return 1.0 / (math.cosh(lam) + x * math.sinh(lam))


def normalized_profile(
    xi: Iterable[float] | np.ndarray, da: float, x: float, u: float
) -> np.ndarray:
    """Exact normalized excess-concentration profile."""

    points = np.asarray(tuple(xi) if not isinstance(xi, np.ndarray) else xi, dtype=float)
    lam = math.sqrt(da)
    top = normalized_top_concentration(da, x, u)
    # The requested sweep has lam<=10, so direct hyperbolic functions are safe.
    return (
        top * np.cosh(lam * points) / math.cosh(lam)
        - u * np.sinh(lam * (1.0 - points)) / math.cosh(lam)
    )


def finite_difference_profile(
    da: float, x: float, u: float, cells: int = 400
) -> tuple[np.ndarray, np.ndarray]:
    """Independent second-order nodal solution of the same linear BVP."""

    if cells < 4:
        raise ValueError("at least four intervals are required")
    lam = math.sqrt(da)
    dx = 1.0 / cells
    matrix = np.zeros((cells + 1, cells + 1))
    rhs = np.zeros(cells + 1)

    # y'(0)=u*lambda, second-order one-sided derivative.
    matrix[0, 0:3] = (-3.0, 4.0, -1.0)
    rhs[0] = 2.0 * dx * u * lam
    for i in range(1, cells):
        matrix[i, i - 1 : i + 2] = (1.0, -2.0 - da * dx * dx, 1.0)
    # y'(1)+lambda*X*y(1)=lambda.
    matrix[cells, cells - 2 : cells + 1] = (1.0, -4.0, 3.0)
    matrix[cells, cells] += 2.0 * dx * lam * x
    rhs[cells] = 2.0 * dx * lam
    xi = np.linspace(0.0, 1.0, cells + 1)
    return xi, np.linalg.solve(matrix, rhs)
