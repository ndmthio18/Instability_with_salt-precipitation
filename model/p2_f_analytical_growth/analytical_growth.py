"""Reaction-free analytical salt-crust surface growth.

This module is an isolated P2-F model.  It does not solve a crust transport
grid, does not change crust porosity, and has no internal precipitation.  The
growth law is derived in ``docs/P2_F_MODEL_FORMULATION.md``; it is not claimed
as an equation stated explicitly by Prat (2026).
"""

from __future__ import annotations

from dataclasses import dataclass, replace
import math


M_WATER = 0.01801528
M_NACL = 0.05844


class ModelDomainExit(RuntimeError):
    """The simplified closure has left its declared validity domain."""


def mass_fraction_from_mole_fraction(x: float) -> float:
    if not 0.0 <= x < 1.0:
        raise ValueError("mole fraction must lie in [0, 1)")
    mixture = (1.0 - x) * M_WATER + x * M_NACL
    return x * M_NACL / mixture


def brine_density(w: float) -> float:
    """Frozen P1-R2 density law [kg/m3], expressed in salt mass fraction."""
    if not 0.0 <= w < 1.0:
        raise ValueError("mass fraction must lie in [0, 1)")
    return 998.203 + 1000.0 * w * (0.670804 + 0.373854 * w)


@dataclass(frozen=True)
class Parameters:
    equilibrium_mass_fraction: float = mass_fraction_from_mole_fraction(0.1)
    porosity: float = 0.43
    crystal_density: float = 2165.0
    reaction_coefficient: float = 2.3e-3
    beta: float = 2.0
    molecular_diffusivity: float = 1.3e-9
    solution_density: float | None = None
    vapour_density: float = 0.0
    minimum_interface_saturation: float = 0.99
    interface_relative_tolerance: float = 0.02
    maximum_solution_peclet: float = 0.1
    maximum_moving_peclet: float = 0.1

    def __post_init__(self) -> None:
        if not 0.0 < self.equilibrium_mass_fraction < 1.0:
            raise ValueError("equilibrium mass fraction must lie in (0, 1)")
        if not 0.0 < self.porosity < 1.0:
            raise ValueError("porosity must lie in (0, 1)")
        for name in ("crystal_density", "reaction_coefficient", "beta",
                     "molecular_diffusivity"):
            value = getattr(self, name)
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be positive and finite")
        if self.vapour_density < 0.0:
            raise ValueError("vapour density must be nonnegative")

    @property
    def rho(self) -> float:
        return (self.solution_density if self.solution_density is not None
                else brine_density(self.equilibrium_mass_fraction))

    @property
    def effective_diffusivity(self) -> float:
        return math.sqrt(self.porosity) * self.molecular_diffusivity

    @property
    def height_coefficient(self) -> float:
        """K [1/m] such that Da_top = K h."""
        return (
            self.beta * (1.0 - self.porosity)
            * (1.0 - self.equilibrium_mass_fraction)
            * self.reaction_coefficient
            / (self.porosity * self.effective_diffusivity)
        )

    def damkohler_top(self, height: float) -> float:
        _require_height(height)
        return self.height_coefficient * height


def _require_height(height: float) -> None:
    if not math.isfinite(height) or height <= 0.0:
        raise ModelDomainExit("an active continuum crust requires finite h > 0")


def growth_prefactor(params: Parameters, evaporation_flux: float) -> float:
    if not math.isfinite(evaporation_flux) or evaporation_flux < 0.0:
        raise ValueError("evaporation flux must be finite and nonnegative")
    c = params.equilibrium_mass_fraction
    return c * evaporation_flux / (
        params.crystal_density * (1.0 - params.porosity) * (1.0 - c)
    )


def growth_velocity(params: Parameters, height: float,
                    evaporation_flux: float) -> float:
    """Reaction-free, negligible-advection surface-growth law [m/s]."""
    da_top = params.damkohler_top(height)
    return growth_prefactor(params, evaporation_flux) * da_top / (1.0 + da_top)


def large_damkohler_top_velocity(params: Parameters,
                                 evaporation_flux: float) -> float:
    """Da_top -> infinity reference, not the production default."""
    return growth_prefactor(params, evaporation_flux)


def _advective_gradient_factor(peclet: float) -> float:
    """Pe exp(Pe)/(exp(Pe)-1), with its regular Pe=0 limit."""
    if abs(peclet) < 1.0e-7:
        return 1.0 + 0.5 * peclet + peclet * peclet / 12.0
    if peclet > 700.0:
        return peclet
    return peclet / (-math.expm1(-peclet))


def solution_mass_flux(params: Parameters, evaporation_flux: float,
                       velocity: float) -> float:
    """rho U from Prat's moving-top total-mass balance [kg/m2/s]."""
    return evaporation_flux + (
        params.rho * params.porosity
        + params.crystal_density * (1.0 - params.porosity)
        - params.vapour_density
    ) * velocity


def solution_peclet(params: Parameters, height: float,
                     evaporation_flux: float, velocity: float) -> float:
    rho_u = solution_mass_flux(params, evaporation_flux, velocity)
    return rho_u * height / (
        params.rho * params.porosity * params.effective_diffusivity
    )


def moving_peclet(params: Parameters, height: float, velocity: float) -> float:
    return velocity * height / params.effective_diffusivity


def independent_advective_velocity(params: Parameters, height: float,
                                    evaporation_flux: float) -> float:
    """Independent reaction-free BVP reference retaining liquid advection.

    The exact constant-property profile solves
    ``U C_z = eps D C_zz``.  Its top-gradient factor is
    ``B(Pe)=Pe exp(Pe)/(exp(Pe)-1)``.  Because U also contains the moving-layer
    demand, the resulting scalar equation is solved by bisection.
    """
    _require_height(height)
    upper = large_damkohler_top_velocity(params, evaporation_flux)
    if upper == 0.0:
        return 0.0
    da_top = params.damkohler_top(height)

    def residual(v: float) -> float:
        pe = solution_peclet(params, height, evaporation_flux, v)
        rhs = upper * da_top / (da_top + _advective_gradient_factor(pe))
        return v - rhs

    low, high = 0.0, upper
    if residual(low) > 0.0 or residual(high) < 0.0:
        raise ModelDomainExit("failed to bracket the reaction-free advective reference")
    for _ in range(120):
        mid = 0.5 * (low + high)
        if residual(mid) > 0.0:
            high = mid
        else:
            low = mid
        if high - low <= 1.0e-14 * max(high, 1.0e-30):
            break
    value = 0.5 * (low + high)
    if abs(residual(value)) > 1.0e-12 * max(value, 1.0e-30):
        raise ModelDomainExit("reaction-free advective reference did not converge")
    return value


def top_mass_fraction(params: Parameters, velocity: float) -> float:
    """Outer kinetic relation corresponding to Prat equation (20)."""
    return params.equilibrium_mass_fraction + (
        params.crystal_density * velocity
        / (params.beta * params.rho * params.reaction_coefficient)
    )


@dataclass
class State:
    height: float = 0.0
    solid_salt: float = 0.0
    dissolved_salt: float = 0.0
    water: float = 0.0
    time: float = 0.0
    active: bool = False
    revision: int = 0
    cumulative_evaporation: float = 0.0

    @property
    def total_salt(self) -> float:
        return self.solid_salt + self.dissolved_salt


def inventories(params: Parameters, height: float) -> tuple[float, float, float]:
    if height < 0.0:
        raise ValueError("height must be nonnegative")
    solid = params.crystal_density * (1.0 - params.porosity) * height
    dissolved = (params.rho * params.porosity
                 * params.equilibrium_mass_fraction * height)
    water = (params.rho * params.porosity
             * (1.0 - params.equilibrium_mass_fraction) * height)
    return solid, dissolved, water


def activation_trial(params: Parameters, seed_height: float) -> State:
    _require_height(seed_height)
    solid, dissolved, water = inventories(params, seed_height)
    return State(height=seed_height, solid_salt=solid,
                 dissolved_salt=dissolved, water=water, active=True)


@dataclass
class Trial:
    base_revision: int
    proposed: State
    evaporation_flux: float
    velocity: float
    salt_flux: float
    water_flux: float
    total_flux: float
    top_mass_fraction: float
    solution_peclet: float
    moving_peclet: float
    committed: bool = False
    discarded: bool = False


def evaluate_trial(state: State, params: Parameters, evaporation_flux: float,
                   dt: float) -> Trial:
    if not state.active:
        raise ModelDomainExit("cannot grow an inactive crust")
    if not math.isfinite(dt) or dt <= 0.0:
        raise ValueError("dt must be positive and finite")
    velocity = growth_velocity(params, state.height, evaporation_flux)
    new_height = state.height + velocity * dt
    solid, dissolved, water = inventories(params, new_height)
    proposed = replace(
        state, height=new_height, solid_salt=solid, dissolved_salt=dissolved,
        water=water, time=state.time + dt,
        cumulative_evaporation=state.cumulative_evaporation + evaporation_flux * dt,
        revision=state.revision + 1,
    )
    salt_flux = ((solid + dissolved) - state.total_salt) / dt
    water_storage = (water - state.water) / dt
    water_flux = evaporation_flux - params.vapour_density * velocity + water_storage
    return Trial(
        base_revision=state.revision,
        proposed=proposed,
        evaporation_flux=evaporation_flux,
        velocity=velocity,
        salt_flux=salt_flux,
        water_flux=water_flux,
        total_flux=salt_flux + water_flux,
        top_mass_fraction=top_mass_fraction(params, velocity),
        solution_peclet=solution_peclet(
            params, state.height, evaporation_flux, velocity),
        moving_peclet=moving_peclet(params, state.height, velocity),
    )


def commit(state: State, trial: Trial) -> None:
    if trial.committed or trial.discarded:
        raise ModelDomainExit("trial is no longer committable")
    if state.revision != trial.base_revision:
        raise ModelDomainExit("stale trial")
    state.__dict__.update(trial.proposed.__dict__)
    trial.committed = True


def rollback(state: State, trial: Trial) -> None:
    if trial.committed:
        raise ModelDomainExit("a committed trial cannot be rolled back")
    if state.revision != trial.base_revision:
        raise ModelDomainExit("persistent state changed before rollback")
    trial.discarded = True


def validity(params: Parameters, state: State, trial: Trial,
             interface_mass_fraction: float,
             interface_saturation: float) -> dict[str, float | bool]:
    relative_interface_error = abs(
        interface_mass_fraction / params.equilibrium_mass_fraction - 1.0
    )
    return {
        "interface_relative_error": relative_interface_error,
        "interface_equilibrium_compatible": (
            relative_interface_error <= params.interface_relative_tolerance),
        "interface_saturation": interface_saturation,
        "interface_saturation_valid": (
            interface_saturation >= params.minimum_interface_saturation),
        "solution_peclet": abs(trial.solution_peclet),
        "solution_advection_negligible": (
            abs(trial.solution_peclet) <= params.maximum_solution_peclet),
        "moving_peclet": abs(trial.moving_peclet),
        "quasi_steady_moving_boundary": (
            abs(trial.moving_peclet) <= params.maximum_moving_peclet),
        "connected_saturated_crust_assumed": True,
        "valid": bool(
            relative_interface_error <= params.interface_relative_tolerance
            and interface_saturation >= params.minimum_interface_saturation
            and abs(trial.solution_peclet) <= params.maximum_solution_peclet
            and abs(trial.moving_peclet) <= params.maximum_moving_peclet
        ),
    }


def analytical_height(params: Parameters, initial_height: float,
                      evaporation_flux: float, elapsed: float) -> float:
    """Independent exact implicit solution for constant evaporation."""
    _require_height(initial_height)
    if elapsed < 0.0:
        raise ValueError("elapsed time must be nonnegative")
    g = growth_prefactor(params, evaporation_flux)
    if g == 0.0 or elapsed == 0.0:
        return initial_height
    k = params.height_coefficient

    def residual(h: float) -> float:
        return ((h - initial_height) + math.log(h / initial_height) / k
                - g * elapsed)

    low = initial_height
    high = initial_height + g * elapsed + 1.0 / k
    while residual(high) < 0.0:
        high *= 2.0
    for _ in range(160):
        mid = 0.5 * (low + high)
        if residual(mid) > 0.0:
            high = mid
        else:
            low = mid
    return 0.5 * (low + high)
