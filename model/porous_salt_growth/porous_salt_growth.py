"""Standalone reduced model for a wet porous NaCl layer.

The module contains no soil saturation, soil diffusion, evaporation, or CT
observation model. A caller supplies the water evaporation mass flux ``j``
and either the net salt-supply ratio ``s`` or the downward salt flux ``B``.

``finite_layer_response`` solves the fixed-geometry reaction--diffusion
problem. ``thick_layer_response`` adds the conditional moving-bottom closure.
All fluxes and rates are per unit projected area. Positive bottom retreat is
upward, so ``dh/dt = top_velocity - bottom_retreat_velocity``.
"""

from __future__ import annotations

from dataclasses import dataclass
import math
from typing import Iterable, Optional


def _finite(name: str, value: float) -> float:
    value = float(value)
    if not math.isfinite(value):
        raise ValueError(f"{name} must be finite")
    return value


@dataclass(frozen=True)
class PorousSaltModel:
    """Reduced salt-layer parameters with ell and X=sqrt(Da*) explicit."""

    equilibrium_mass_fraction: float
    crystal_density_kg_m3: float
    top_porosity: float
    reaction_length_m: float
    sqrt_da_star: float

    def __post_init__(self) -> None:
        ceq = _finite("equilibrium_mass_fraction", self.equilibrium_mass_fraction)
        rho_cr = _finite("crystal_density_kg_m3", self.crystal_density_kg_m3)
        eps = _finite("top_porosity", self.top_porosity)
        ell = _finite("reaction_length_m", self.reaction_length_m)
        x = _finite("sqrt_da_star", self.sqrt_da_star)
        if not 0.0 < ceq < 1.0:
            raise ValueError("equilibrium_mass_fraction must lie in (0, 1)")
        if rho_cr <= 0.0 or ell <= 0.0 or x <= 0.0:
            raise ValueError("density, reaction length, and sqrt(Da*) must be positive")
        if not 0.0 < eps < 1.0:
            raise ValueError("top_porosity must lie in (0, 1)")

    @classmethod
    def from_transport_parameters(
        cls,
        *,
        equilibrium_mass_fraction: float,
        crystal_density_kg_m3: float,
        top_porosity: float,
        effective_diffusivity_m2_s: float,
        reaction_coefficient_m_s: float,
        specific_area_m_inv: float,
        beta: float,
    ) -> "PorousSaltModel":
        """Derive ell and X from one internally consistent parameter set."""
        ceq = _finite("equilibrium_mass_fraction", equilibrium_mass_fraction)
        eps = _finite("top_porosity", top_porosity)
        diffusivity = _finite("effective_diffusivity_m2_s", effective_diffusivity_m2_s)
        kr = _finite("reaction_coefficient_m_s", reaction_coefficient_m_s)
        av = _finite("specific_area_m_inv", specific_area_m_inv)
        beta = _finite("beta", beta)
        if diffusivity <= 0.0 or kr <= 0.0 or av <= 0.0 or beta <= 0.0:
            raise ValueError("transport and kinetic parameters must be positive")
        if not 0.0 < ceq < 1.0 or not 0.0 < eps < 1.0:
            raise ValueError("mass fraction and porosity must lie in (0, 1)")
        ell = math.sqrt(eps * diffusivity / ((1.0 - ceq) * av * kr))
        x = beta * (1.0 - eps) * (1.0 - ceq) * kr * ell / (eps * diffusivity)
        return cls(ceq, crystal_density_kg_m3, eps, ell, x)

    @property
    def top_skeleton_density_kg_m3(self) -> float:
        return self.crystal_density_kg_m3 * (1.0 - self.top_porosity)

    @property
    def wake_porosity(self) -> float:
        return self.top_porosity - (1.0 - self.top_porosity) / self.sqrt_da_star

    @property
    def thick_connected(self) -> bool:
        return self.wake_porosity > 0.0

    def reference_salt_flux(self, evaporation_flux_kg_m2_s: float) -> float:
        j = _finite("evaporation_flux_kg_m2_s", evaporation_flux_kg_m2_s)
        if j < 0.0:
            raise ValueError("evaporation flux must be non-negative")
        return self.equilibrium_mass_fraction * j / (1.0 - self.equilibrium_mass_fraction)

    def _forcing(
        self,
        evaporation_flux_kg_m2_s: float,
        *,
        supply_ratio: Optional[float],
        reverse_flux_kg_m2_s: Optional[float],
    ) -> tuple[float, float, float, float]:
        if (supply_ratio is None) == (reverse_flux_kg_m2_s is None):
            raise ValueError("provide exactly one of supply_ratio or reverse_flux_kg_m2_s")
        j = _finite("evaporation_flux_kg_m2_s", evaporation_flux_kg_m2_s)
        if j < 0.0:
            raise ValueError("evaporation flux must be non-negative")
        if supply_ratio is not None:
            s = _finite("supply_ratio", supply_ratio)
            if not 0.0 <= s <= 1.0:
                raise ValueError("the growth module requires 0 <= supply_ratio <= 1")
            u = 1.0 - s
            reverse_flux = u * self.equilibrium_mass_fraction * j
        else:
            reverse_flux = _finite("reverse_flux_kg_m2_s", reverse_flux_kg_m2_s)
            if reverse_flux < 0.0:
                raise ValueError("reverse flux must be non-negative")
            if j == 0.0:
                if reverse_flux != 0.0:
                    raise ValueError("non-zero reverse flux is incompatible with zero evaporation")
                u, s = 0.0, 1.0
            else:
                u = reverse_flux / (self.equilibrium_mass_fraction * j)
                s = 1.0 - u
                if not 0.0 <= s <= 1.0:
                    raise ValueError("reverse flux implies a supply ratio outside [0, 1]")
        return j, s, u, reverse_flux

    def prat_compatible_supply_ratio(self, height_m: float) -> float:
        """Supply ratio implied by Prat's fixed-equilibrium bottom condition."""
        height = _finite("height_m", height_m)
        if height < 0.0:
            raise ValueError("height must be non-negative")
        h = height / self.reaction_length_m
        x = self.sqrt_da_star
        if h > 50.0:
            u_critical = 2.0 * math.exp(-h) / (1.0 + x)
        else:
            u_critical = 1.0 / (math.cosh(h) + x * math.sinh(h))
        return 1.0 - u_critical

    def prat_finite_response(
        self, height_m: float, evaporation_flux_kg_m2_s: float
    ) -> "FiniteLayerResponse":
        return self.finite_layer_response(
            height_m,
            evaporation_flux_kg_m2_s,
            supply_ratio=self.prat_compatible_supply_ratio(height_m),
        )

    def finite_layer_response(
        self,
        height_m: float,
        evaporation_flux_kg_m2_s: float,
        *,
        supply_ratio: Optional[float] = None,
        reverse_flux_kg_m2_s: Optional[float] = None,
    ) -> "FiniteLayerResponse":
        """Return the fixed-geometry finite-layer analytical response.

        It solves ``y''-y=0``, ``y'(0)=u``, ``y'(H)+X*y(H)=1``.
        No bottom-motion law is assumed here.
        """
        height = _finite("height_m", height_m)
        if height < 0.0:
            raise ValueError("height must be non-negative")
        j, s, u, reverse_flux = self._forcing(
            evaporation_flux_kg_m2_s,
            supply_ratio=supply_ratio,
            reverse_flux_kg_m2_s=reverse_flux_kg_m2_s,
        )
        h = height / self.reaction_length_m
        x = self.sqrt_da_star
        e = math.exp(-h)
        denominator = (1.0 + x) + e * e * (x - 1.0)
        a = (1.0 + u * e * (x - 1.0)) / denominator
        b = a * e - u
        y_top = a + b * e
        eta_top = x * y_top
        reference = self.reference_salt_flux(j)
        top_solid = eta_top * reference
        top_velocity = top_solid / self.top_skeleton_density_kg_m3

        def primitive(position: float) -> float:
            return a * math.exp(position - h) - b * math.exp(-position)

        cuts = [0.0, h]
        if a * b < 0.0 and a != 0.0:
            root = 0.5 * (h + math.log(-b / a))
            if 0.0 < root < h:
                cuts.insert(1, root)
        positive_integral = 0.0
        negative_integral = 0.0
        for left, right in zip(cuts[:-1], cuts[1:]):
            integral = primitive(right) - primitive(left)
            middle = 0.5 * (left + right)
            y_middle = a * math.exp(middle - h) + b * math.exp(-middle)
            if y_middle >= 0.0:
                positive_integral += max(integral, 0.0)
            else:
                negative_integral += max(-integral, 0.0)
        internal_precipitation = reference * positive_integral
        internal_dissolution = reference * negative_integral
        bulk_net = internal_precipitation - internal_dissolution
        net_solid = top_solid + bulk_net
        expected_net = s * reference
        tolerance = 2.0e-12 * max(reference, 1.0e-30)
        if abs(net_solid - expected_net) > tolerance:
            raise ArithmeticError("finite-layer salt budget did not close")
        return FiniteLayerResponse(
            height_m=height,
            dimensionless_height=h,
            supply_ratio=s,
            normalized_reverse_flux=u,
            reverse_flux_kg_m2_s=reverse_flux,
            reference_salt_flux_kg_m2_s=reference,
            coefficient_a=a,
            coefficient_b=b,
            top_partition=eta_top,
            top_velocity_m_s=top_velocity,
            top_solid_rate_kg_m2_s=top_solid,
            internal_precipitation_kg_m2_s=internal_precipitation,
            internal_dissolution_kg_m2_s=internal_dissolution,
            internal_net_rate_kg_m2_s=bulk_net,
            net_solid_rate_kg_m2_s=net_solid,
        )

    def thick_layer_response(
        self,
        evaporation_flux_kg_m2_s: float,
        *,
        supply_ratio: Optional[float] = None,
        reverse_flux_kg_m2_s: Optional[float] = None,
    ) -> "ThickLayerResponse":
        """Return the thick wet-layer response with a moving dissolution front."""
        if not self.thick_connected:
            raise ValueError(
                "thick moving-front closure is invalid because wake_porosity <= 0"
            )
        j, s, u, reverse_flux = self._forcing(
            evaporation_flux_kg_m2_s,
            supply_ratio=supply_ratio,
            reverse_flux_kg_m2_s=reverse_flux_kg_m2_s,
        )
        reference = self.reference_salt_flux(j)
        x = self.sqrt_da_star
        top_partition = x / (1.0 + x)
        top_solid = top_partition * reference
        upper_precipitation = reference / (1.0 + x)
        bottom_dissolution = u * reference
        net_solid = s * reference
        top_velocity = top_solid / self.top_skeleton_density_kg_m3
        wake_density = self.crystal_density_kg_m3 * (1.0 - self.wake_porosity)
        bottom_retreat = bottom_dissolution / wake_density
        thickness_velocity = top_velocity - bottom_retreat
        tolerance = 2.0e-12 * max(reference, 1.0e-30)
        if abs(top_solid + upper_precipitation - bottom_dissolution - net_solid) > tolerance:
            raise ArithmeticError("thick-layer salt budget did not close")
        velocity_tolerance = 2.0e-12 * max(abs(top_velocity), 1.0e-30)
        if abs(thickness_velocity - s * top_velocity) > velocity_tolerance:
            raise ArithmeticError("thick-layer kinematic closure did not close")
        return ThickLayerResponse(
            supply_ratio=s,
            normalized_reverse_flux=u,
            reverse_flux_kg_m2_s=reverse_flux,
            reference_salt_flux_kg_m2_s=reference,
            top_partition=top_partition,
            wake_porosity=self.wake_porosity,
            top_velocity_m_s=top_velocity,
            bottom_retreat_velocity_m_s=bottom_retreat,
            thickness_velocity_m_s=thickness_velocity,
            top_solid_rate_kg_m2_s=top_solid,
            internal_precipitation_kg_m2_s=upper_precipitation,
            internal_dissolution_kg_m2_s=bottom_dissolution,
            net_solid_rate_kg_m2_s=net_solid,
        )

    def integrate_thick_layer(
        self,
        times_s: Iterable[float],
        evaporation_flux_kg_m2_s: Iterable[float],
        supply_ratio: Iterable[float],
        *,
        initial_top_m: float = 0.0,
        initial_base_m: float = 0.0,
    ) -> list["GrowthState"]:
        """Integrate the thick closure for external ``j(t), s(t)`` histories."""
        times = [_finite("time", value) for value in times_s]
        fluxes = [_finite("evaporation flux", value) for value in evaporation_flux_kg_m2_s]
        supplies = [_finite("supply ratio", value) for value in supply_ratio]
        if not (len(times) == len(fluxes) == len(supplies)) or len(times) < 2:
            raise ValueError("times, evaporation fluxes, and supply ratios need equal length >= 2")
        if any(right <= left for left, right in zip(times[:-1], times[1:])):
            raise ValueError("times must be strictly increasing")
        responses = [
            self.thick_layer_response(j, supply_ratio=s)
            for j, s in zip(fluxes, supplies)
        ]
        top = _finite("initial_top_m", initial_top_m)
        base = _finite("initial_base_m", initial_base_m)
        if top < base:
            raise ValueError("initial_top_m must be at or above initial_base_m")
        cumulative_top = 0.0
        cumulative_internal = 0.0
        cumulative_dissolution = 0.0
        cumulative_net = 0.0
        states = [GrowthState(times[0], top, base, top - base, 0.0, 0.0, 0.0, 0.0)]
        for index in range(1, len(times)):
            dt = times[index] - times[index - 1]
            previous, current = responses[index - 1], responses[index]

            def average(field: str) -> float:
                return 0.5 * (getattr(previous, field) + getattr(current, field))

            top += dt * average("top_velocity_m_s")
            base += dt * average("bottom_retreat_velocity_m_s")
            cumulative_top += dt * average("top_solid_rate_kg_m2_s")
            cumulative_internal += dt * average("internal_precipitation_kg_m2_s")
            cumulative_dissolution += dt * average("internal_dissolution_kg_m2_s")
            cumulative_net += dt * average("net_solid_rate_kg_m2_s")
            states.append(
                GrowthState(
                    times[index], top, base, top - base,
                    cumulative_top, cumulative_internal,
                    cumulative_dissolution, cumulative_net,
                )
            )
        return states


@dataclass(frozen=True)
class FiniteLayerResponse:
    height_m: float
    dimensionless_height: float
    supply_ratio: float
    normalized_reverse_flux: float
    reverse_flux_kg_m2_s: float
    reference_salt_flux_kg_m2_s: float
    coefficient_a: float
    coefficient_b: float
    top_partition: float
    top_velocity_m_s: float
    top_solid_rate_kg_m2_s: float
    internal_precipitation_kg_m2_s: float
    internal_dissolution_kg_m2_s: float
    internal_net_rate_kg_m2_s: float
    net_solid_rate_kg_m2_s: float

    def dimensionless_concentration(self, position_over_ell: float) -> float:
        position = _finite("position_over_ell", position_over_ell)
        if not 0.0 <= position <= self.dimensionless_height:
            raise ValueError("position must lie inside the layer")
        return (
            self.coefficient_a * math.exp(position - self.dimensionless_height)
            + self.coefficient_b * math.exp(-position)
        )


@dataclass(frozen=True)
class ThickLayerResponse:
    supply_ratio: float
    normalized_reverse_flux: float
    reverse_flux_kg_m2_s: float
    reference_salt_flux_kg_m2_s: float
    top_partition: float
    wake_porosity: float
    top_velocity_m_s: float
    bottom_retreat_velocity_m_s: float
    thickness_velocity_m_s: float
    top_solid_rate_kg_m2_s: float
    internal_precipitation_kg_m2_s: float
    internal_dissolution_kg_m2_s: float
    net_solid_rate_kg_m2_s: float


@dataclass(frozen=True)
class GrowthState:
    time_s: float
    top_position_m: float
    base_position_m: float
    thickness_m: float
    cumulative_top_solid_kg_m2: float
    cumulative_internal_precipitation_kg_m2: float
    cumulative_internal_dissolution_kg_m2: float
    cumulative_net_solid_kg_m2: float
