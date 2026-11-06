"""P2-A surface patch: one-way diagnostic chain and conservative base fluxes.

The chain required by P2-A item 7 is

    soil surface salinity -> water activity -> evaporation flux
                          -> Prat growth velocity -> crust height

`SurfacePatch.diagnose` evaluates exactly that chain and returns every
intermediate quantity. `SurfacePatch.integrate` advances the crust inventories
in time and additionally reconstructs the soil-side water and NaCl fluxes
(P2-A item 8) from Prat's moving-interface balances and the reduced crust
transport solution.

Evaporation-interface options (P2-A item 9)
-------------------------------------------
``EvaporationInterface.SOIL_TOP`` (option A)
    The water activity driving the gas film is evaluated at the **soil-top**
    (crust base) NaCl mass fraction. This is the simplified approximation: it
    ignores the supersaturation that the crust transport solution builds up
    between base and top.

``EvaporationInterface.CRUST_TOP`` (option B)
    The water activity is evaluated at the **crust-top** mass fraction
    ``C_top = C_eq + C*_top`` obtained from the reduced crust solution. Because
    ``C*_top`` itself depends on the evaporation flux through ``Pe``, option B
    is an implicit relation and is solved here by a bracketed bisection on the
    flux. No under-relaxation constant is tuned; the bracket is derived from
    the monotonicity of both branches.

Before an explicitly declared seed event there is no crust, both options
coincide, and the patch reports the ordinary soil evaporation boundary with
zero net salt flux, matching the P1-R2 top condition.

Conservation design
-------------------
The crust inventories ``M_s,solid``, ``M_s,d`` and ``M_w,c`` are **states**
advanced from the same rates that enter the reported interface fluxes:

    d h        / dt = v_s
    d M_s,solid/ dt = (1-eps_c) rho_cr v_s + R_int
    d M_s,d    / dt = eps_c S_c rho C_top v_s
    d M_w,c    / dt = eps_c S_c rho (1-C_top) v_s

so that

    F_s,soil = d(M_s,solid + M_s,d)/dt
    F_w,soil = j - rho_v v_s + d M_w,c/dt

hold exactly, up to the time-integration error, by construction. The quantity
that the quasi-steady reduction *drops* -- the change of the dissolved profile
shape -- is not hidden: `integrate` also evaluates the profile-diagnostic
inventory ``M_s,d^QS`` and reports the running difference as
``dissolved_storage_gap``. P2-A item 8's prohibition ("do not retain zero soil
salt flux while accumulating solid salt in the crust") is satisfied because
``F_s,soil`` is strictly positive whenever solid accumulates.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from enum import Enum

from . import prat_crust as pc
from .film_evaporation import GasFilm
from .water_activity import water_activity_from_mass_fraction


class EvaporationInterface(str, Enum):
    """Which concentration drives the evaporating-surface water activity."""

    SOIL_TOP = "soil_top"     # option A
    CRUST_TOP = "crust_top"   # option B


@dataclass(frozen=True)
class ChainDiagnostic:
    """One evaluation of the P2-A item 7 one-way chain."""

    time: float
    soil_mass_fraction: float
    soil_mole_fraction: float
    molality: float
    water_activity: float
    evaporation_interface: str
    interface_mass_fraction: float
    interface_activity: float
    driving_mass_fraction_difference: float
    evaporation_flux: float
    active: bool
    height: float
    growth_velocity: float
    top_supersaturation: float
    damkohler: float
    damkohler_star: float
    peclet: float
    efficiency: float
    fixed_point_iterations: int = 0
    fixed_point_residual: float = 0.0


@dataclass
class SurfacePatch:
    """One projected surface patch: gas film + declared seed + Prat crust."""

    film: GasFilm
    crust: pc.CrustParameters
    seed: pc.SeedEvent
    law: pc.GrowthLaw = pc.GrowthLaw.FINITE_DA
    interface: EvaporationInterface = EvaporationInterface.SOIL_TOP
    state: pc.CrustState = field(default_factory=pc.CrustState)
    fixed_point_tolerance: float = 1.0e-14
    fixed_point_max_iterations: int = 200

    # -- evaporation -------------------------------------------------------

    def _flux_soil_top(self, soil_mass_fraction: float) -> float:
        return self.film.flux_from_mass_fraction(soil_mass_fraction)

    def _flux_crust_top(self, height: float) -> tuple[float, float, int, float]:
        """Solve j = film(a_w(C_eq + C*_top(j))) by bracketed bisection.

        Returns ``(j, C_top, iterations, residual)``.

        Monotonicity: ``C*_top`` grows linearly with ``j`` (through ``Pe``) and
        the film flux decreases with ``C_top``, so ``g(j) = film(C_top(j)) - j``
        is strictly decreasing. ``g(0) >= 0`` whenever the film at equilibrium
        concentration evaporates, and ``g(j_eq) < 0`` for the equilibrium-
        concentration flux ``j_eq``, which therefore brackets the root.
        """
        c_eq = self.crust.equilibrium_mass_fraction

        def residual(j: float) -> tuple[float, float]:
            c_top = pc.top_mass_fraction(self.crust, height, j, self.law)
            if c_top >= 1.0:
                raise pc.CrustModelDomainExit(
                    "crust-top NaCl mass fraction reached unity; the dilute "
                    "linearised reaction closure is invalid"
                )
            try:
                return self.film.flux_from_mass_fraction(c_top) - j, c_top
            except ValueError as error:
                raise pc.CrustModelDomainExit(
                    "the crust-top supersaturation left the validity range of "
                    "the inherited NaCl activity correlation "
                    f"(C_top = {c_top!r}): {error}. Option B is reported as a "
                    "model-domain exit rather than extrapolated."
                ) from error

        lo = 0.0
        g_lo, c_lo = residual(lo)
        if g_lo <= 0.0:
            return lo, c_lo, 0, g_lo
        hi = self.film.flux_from_mass_fraction(c_eq)
        g_hi, c_hi = residual(hi)
        guard = 0
        while g_hi > 0.0 and guard < 60:
            hi *= 2.0
            g_hi, c_hi = residual(hi)
            guard += 1
        if g_hi > 0.0:
            raise pc.CrustModelDomainExit(
                "failed to bracket the option-B evaporation fixed point"
            )
        iterations = 0
        c_mid = c_hi
        mid = hi
        while iterations < self.fixed_point_max_iterations:
            mid = 0.5 * (lo + hi)
            g_mid, c_mid = residual(mid)
            iterations += 1
            if g_mid > 0.0:
                lo = mid
            else:
                hi = mid
            if hi - lo <= self.fixed_point_tolerance * max(hi, 1.0e-30):
                break
        g_mid, c_mid = residual(mid)
        return mid, c_mid, iterations, g_mid

    # -- chain -------------------------------------------------------------

    def diagnose(self, soil_mass_fraction: float, time: float | None = None) -> ChainDiagnostic:
        """Evaluate the one-way chain at the given soil-top NaCl mass fraction."""
        from .thermo import (
            molality_from_mass_fraction,
            mole_fraction_from_mass_fraction,
        )

        t = self.state.time if time is None else time
        molality = molality_from_mass_fraction(soil_mass_fraction)
        a_w_soil = water_activity_from_mass_fraction(soil_mass_fraction)

        iterations = 0
        residual = 0.0
        if not self.state.active:
            interface_c = soil_mass_fraction
            j = self._flux_soil_top(soil_mass_fraction)
        elif self.interface is EvaporationInterface.SOIL_TOP:
            interface_c = soil_mass_fraction
            j = self._flux_soil_top(soil_mass_fraction)
        else:
            j, interface_c, iterations, residual = self._flux_crust_top(
                self.state.height
            )

        if self.state.active:
            h = self.state.height
            v_s = pc.growth_velocity(self.crust, h, j, self.law)
            c_star_top = pc.top_supersaturation(self.crust, h, j, self.law)
            da = self.crust.damkohler(h)
            pe = self.crust.peclet(h, j)
            eff = (
                pc.efficiency_large_damkohler(self.crust.damkohler_star)
                if self.law is pc.GrowthLaw.LARGE_DA
                else pc.efficiency(da, self.crust.damkohler_star)
            )
        else:
            h = 0.0
            v_s = 0.0
            c_star_top = 0.0
            da = 0.0
            pe = 0.0
            eff = 0.0

        return ChainDiagnostic(
            time=t,
            soil_mass_fraction=soil_mass_fraction,
            soil_mole_fraction=mole_fraction_from_mass_fraction(soil_mass_fraction),
            molality=molality,
            water_activity=a_w_soil,
            evaporation_interface=(
                self.interface.value if self.state.active else "pre_activation"
            ),
            interface_mass_fraction=interface_c,
            interface_activity=water_activity_from_mass_fraction(interface_c),
            driving_mass_fraction_difference=(
                self.film.driving_mass_fraction_difference(
                    water_activity_from_mass_fraction(interface_c)
                )
            ),
            evaporation_flux=j,
            active=self.state.active,
            height=h,
            growth_velocity=v_s,
            top_supersaturation=c_star_top,
            damkohler=da,
            damkohler_star=self.crust.damkohler_star,
            peclet=pe,
            efficiency=eff,
            fixed_point_iterations=iterations,
            fixed_point_residual=residual,
        )

    # -- activation --------------------------------------------------------

    def consider_activation(self, soil_mass_fraction: float, time: float) -> bool:
        """Apply the declared seed rule. Returns True on the activating call.

        Eligibility (the solubility crossing) and activation are recorded
        separately so that a thermodynamic threshold is never reported as a
        nucleation prediction.
        """
        if self.state.active:
            return False
        if not self.seed.eligible(soil_mass_fraction):
            return False
        if self.state.eligibility_time is None:
            self.state.eligibility_time = time
        j = self._flux_soil_top(soil_mass_fraction)
        pc.activate(self.state, self.crust, self.seed, j, time, self.law)
        return True

    # -- time integration --------------------------------------------------

    def _rates(self, soil_mass_fraction: float) -> tuple[pc.InterfaceFluxes, dict]:
        chain = self.diagnose(soil_mass_fraction)
        fluxes = pc.interface_fluxes(
            self.crust, self.state.height, chain.evaporation_flux, self.law
        )
        rates = {
            "height": fluxes.growth_velocity,
            "solid_salt": fluxes.top_skeleton_rate + fluxes.internal_reaction,
            "dissolved_salt": fluxes.dissolved_salt_rate,
            "water": fluxes.water_storage_rate,
            "cumulative_salt_from_soil": fluxes.salt_flux,
            "cumulative_water_from_soil": fluxes.water_flux,
            "cumulative_evaporation": fluxes.evaporation_flux,
            "cumulative_vapour_sweep": -fluxes.vapour_sweep_rate,
        }
        return fluxes, rates

    def step(
        self,
        soil_mass_fraction: float,
        dt: float,
        scheme: str = "heun",
    ) -> dict:
        """Advance the crust state by ``dt``; returns a record for the history.

        ``scheme`` is ``"euler"`` (explicit Euler, first order) or ``"heun"``
        (explicit trapezoid, second order). Both are explicit so that the
        reported fluxes are exactly the fluxes used to advance the inventories.
        """
        if dt <= 0.0 or not math.isfinite(dt):
            raise ValueError("time step must be positive and finite")
        if not self.state.active:
            chain = self.diagnose(soil_mass_fraction)
            self.state.time += dt
            self.state.cumulative_evaporation += chain.evaporation_flux * dt
            # No crust: the soil water flux goes straight to the atmosphere and
            # the soil salt flux is zero, exactly the P1-R2 top condition.
            self.state.cumulative_water_from_soil += chain.evaporation_flux * dt
            record = _inactive_record(self.state, chain, dt)
            self.state.history.append(record)
            return record

        start = self.state.copy()
        fluxes_a, rates_a = self._rates(soil_mass_fraction)

        if scheme == "euler":
            rates = rates_a
            fluxes = fluxes_a
        elif scheme == "heun":
            for key, rate in rates_a.items():
                setattr(self.state, key, getattr(start, key) + rate * dt)
            self.state.time = start.time + dt
            try:
                fluxes_b, rates_b = self._rates(soil_mass_fraction)
            finally:
                for key in rates_a:
                    setattr(self.state, key, getattr(start, key))
                self.state.time = start.time
            rates = {k: 0.5 * (rates_a[k] + rates_b[k]) for k in rates_a}
            fluxes = fluxes_a
        else:
            raise ValueError("scheme must be 'euler' or 'heun'")

        for key, rate in rates.items():
            setattr(self.state, key, getattr(start, key) + rate * dt)
        self.state.time = start.time + dt

        chain = self.diagnose(soil_mass_fraction)
        qs_dissolved = pc.dissolved_salt_inventory(
            self.crust, self.state.height, chain.evaporation_flux, self.law
        )
        record = {
            "time_s": self.state.time,
            "dt_s": dt,
            "active": True,
            "soil_mass_fraction": soil_mass_fraction,
            "soil_water_activity": chain.water_activity,
            "interface_mass_fraction": chain.interface_mass_fraction,
            "interface_water_activity": chain.interface_activity,
            "evaporation_flux_kg_m2_s": chain.evaporation_flux,
            "growth_velocity_m_s": rates["height"],
            "crust_height_m": self.state.height,
            "crust_porosity": self.crust.crust_porosity,
            "damkohler": chain.damkohler,
            "damkohler_top": self.crust.damkohler_top(self.state.height),
            "damkohler_star": chain.damkohler_star,
            "peclet": chain.peclet,
            "efficiency": chain.efficiency,
            "top_supersaturation": chain.top_supersaturation,
            "internal_reaction_kg_m2_s": fluxes.internal_reaction,
            "salt_flux_soil_kg_m2_s": rates["cumulative_salt_from_soil"],
            "water_flux_soil_kg_m2_s": rates["cumulative_water_from_soil"],
            "salt_flux_base_route_kg_m2_s": fluxes.salt_flux_base_route,
            "route_difference_kg_m2_s": fluxes.route_difference,
            "route_difference_identity_kg_m2_s": fluxes.route_difference_identity,
            "vapour_sweep_kg_m2_s": fluxes.vapour_sweep_rate,
            "solid_salt_kg_m2": self.state.solid_salt,
            "dissolved_salt_kg_m2": self.state.dissolved_salt,
            "dissolved_salt_quasi_steady_kg_m2": qs_dissolved,
            "dissolved_storage_gap_kg_m2": self.state.dissolved_salt - qs_dissolved,
            "water_kg_m2": self.state.water,
            "cumulative_salt_from_soil_kg_m2": self.state.cumulative_salt_from_soil,
            "cumulative_water_from_soil_kg_m2": self.state.cumulative_water_from_soil,
            "cumulative_evaporation_kg_m2": self.state.cumulative_evaporation,
        }
        record.update(
            {
                f"validity_{k}": v
                for k, v in pc.validity_flags(
                    self.crust, self.state.height, chain.evaporation_flux
                ).items()
            }
        )
        self.state.history.append(record)
        return record

    # -- budgets -----------------------------------------------------------

    def budget(self, initial: pc.CrustState) -> dict:
        """Closure residuals of the component budgets since ``initial``."""
        s = self.state
        d_solid = s.solid_salt - initial.solid_salt
        d_dissolved = s.dissolved_salt - initial.dissolved_salt
        d_water = s.water - initial.water
        salt_transfer = s.cumulative_salt_from_soil - initial.cumulative_salt_from_soil
        water_transfer = (
            s.cumulative_water_from_soil - initial.cumulative_water_from_soil
        )
        evaporated = s.cumulative_evaporation - initial.cumulative_evaporation
        swept = s.cumulative_vapour_sweep - initial.cumulative_vapour_sweep
        salt_residual = salt_transfer - (d_solid + d_dissolved)
        water_residual = water_transfer - (evaporated - swept + d_water)
        salt_scale = max(abs(salt_transfer), abs(d_solid + d_dissolved), 1.0e-300)
        water_scale = max(abs(water_transfer), abs(evaporated), 1.0e-300)
        return {
            "delta_solid_salt_kg_m2": d_solid,
            "delta_dissolved_salt_kg_m2": d_dissolved,
            "delta_water_kg_m2": d_water,
            "integrated_salt_from_soil_kg_m2": salt_transfer,
            "integrated_water_from_soil_kg_m2": water_transfer,
            "integrated_evaporation_kg_m2": evaporated,
            "vapour_swept_kg_m2": swept,
            "seed_salt_kg_m2": s.seed_salt - initial.seed_salt,
            "seed_water_kg_m2": s.seed_water - initial.seed_water,
            "salt_residual_kg_m2": salt_residual,
            "water_residual_kg_m2": water_residual,
            "salt_residual_relative": abs(salt_residual) / salt_scale,
            "water_residual_relative": abs(water_residual) / water_scale,
        }


def _inactive_record(state: pc.CrustState, chain: ChainDiagnostic, dt: float) -> dict:
    return {
        "time_s": state.time,
        "dt_s": dt,
        "active": False,
        "soil_mass_fraction": chain.soil_mass_fraction,
        "soil_water_activity": chain.water_activity,
        "interface_mass_fraction": chain.interface_mass_fraction,
        "interface_water_activity": chain.interface_activity,
        "evaporation_flux_kg_m2_s": chain.evaporation_flux,
        "growth_velocity_m_s": 0.0,
        "crust_height_m": 0.0,
        "salt_flux_soil_kg_m2_s": 0.0,
        "water_flux_soil_kg_m2_s": chain.evaporation_flux,
        "solid_salt_kg_m2": 0.0,
        "dissolved_salt_kg_m2": 0.0,
        "water_kg_m2": 0.0,
        "cumulative_evaporation_kg_m2": state.cumulative_evaporation,
    }


# --- interface comparison (P2-A item 9) --------------------------------------


@dataclass(frozen=True)
class InterfaceComparison:
    """Quantified difference between evaporation-interface options A and B."""

    height: float
    soil_mass_fraction: float
    flux_a: float
    flux_b: float
    top_mass_fraction_b: float
    top_supersaturation_b: float
    relative_flux_difference: float
    growth_velocity_a: float
    growth_velocity_b: float
    relative_growth_difference: float
    peclet_a: float
    acceptable: bool
    tolerance: float


def compare_interfaces(
    film: GasFilm,
    crust: pc.CrustParameters,
    height: float,
    soil_mass_fraction: float,
    law: pc.GrowthLaw = pc.GrowthLaw.FINITE_DA,
    tolerance: float = 0.01,
) -> InterfaceComparison:
    """Compare option A (soil-top C) against option B (crust-top C).

    ``acceptable`` is True when the relative evaporation-flux difference is at
    or below ``tolerance``. The tolerance is a declared acceptance level, not a
    physical result.
    """
    seed = pc.SeedEvent(
        trigger_mass_fraction=crust.equilibrium_mass_fraction,
        seed_height=height,
        seed_porosity=crust.crust_porosity,
    )
    patch_a = SurfacePatch(
        film=film,
        crust=crust,
        seed=seed,
        law=law,
        interface=EvaporationInterface.SOIL_TOP,
    )
    patch_a.state.active = True
    patch_a.state.height = height
    patch_b = SurfacePatch(
        film=film,
        crust=crust,
        seed=seed,
        law=law,
        interface=EvaporationInterface.CRUST_TOP,
    )
    patch_b.state.active = True
    patch_b.state.height = height

    chain_a = patch_a.diagnose(soil_mass_fraction)
    chain_b = patch_b.diagnose(soil_mass_fraction)
    rel_flux = abs(chain_b.evaporation_flux - chain_a.evaporation_flux) / abs(
        chain_a.evaporation_flux
    )
    rel_growth = abs(chain_b.growth_velocity - chain_a.growth_velocity) / abs(
        chain_a.growth_velocity
    )
    return InterfaceComparison(
        height=height,
        soil_mass_fraction=soil_mass_fraction,
        flux_a=chain_a.evaporation_flux,
        flux_b=chain_b.evaporation_flux,
        top_mass_fraction_b=chain_b.interface_mass_fraction,
        top_supersaturation_b=chain_b.top_supersaturation,
        relative_flux_difference=rel_flux,
        growth_velocity_a=chain_a.growth_velocity,
        growth_velocity_b=chain_b.growth_velocity,
        relative_growth_difference=rel_growth,
        peclet_a=chain_a.peclet,
        acceptable=rel_flux <= tolerance,
        tolerance=tolerance,
    )
