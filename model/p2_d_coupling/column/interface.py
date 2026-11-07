"""P2-D: the soil/crust interface problem and its compatible flux solve.

Task items 6 to 9.

The flux conditions (item 7)
----------------------------
At the interface `z = L` (soil top = crust base) the two sides exchange one
solution-mass flux and one NaCl-mass flux. Writing `w_i` for the NaCl mass
fraction **at the interface**, the crust side gives, from the quasi-steady
reduced transport problem with `C(0) = w_i` and the moving-interface balances
(`graded_crust.solve_profile`),

    F_tot,crust(w_i) = rho U(0)
                     = [rho eps_0 S + rho_cr (1-eps_0) - rho_v] v_s + j
                       + R_int (1 - rho S / rho_cr)
    F_s,crust(w_i)   = w_i rho U(0) - rho eps_0 S D_eff dC/dz|_0

in which `j` is itself a function of `w_i`, because the gas film is driven by
the water activity of the **crust top** concentration, which the same solve
produces. That inner dependence is closed by bisection, exactly as P2-A's
`SurfacePatch` does for its option-B interface.

The soil side gives, from the same face reconstruction the column uses
everywhere else, with the interface half a cell above the top cell centre,

    F_s,soil(w_i) = F_tot w_up - rho_face D_face (w_i - w_N) / (dz/2) .

**One unknown, one equation.** The task forbids prescribing both the
interface concentration and the interface salt flux. Neither is prescribed
here: `w_i` is solved for from

    g(w_i) = F_s,soil(w_i) - F_s,crust(w_i) = 0 ,

and the salt flux is whatever that `w_i` implies. The total mass flux is not
an independent condition: it is handed to the soil as the crust's water
demand, so the two sides share one solution flux by construction.

What is independent, and therefore an actual check (item 8)
-----------------------------------------------------------
Three routes to the interface salt flux are evaluated. Two of them are
equal by construction and the third is not:

* route 1, **crust transport**: `w_i rho U(0) + D_base`, from the crust's
  concentration profile;
* route 3, **soil component flux**: the soil's own face reconstruction. Equal
  to route 1 at the solved `w_i` to the bisection tolerance, by definition of
  `g`;
* route 2, **crust inventory**: `d(M_solid + M_dissolved)/dt` measured from
  the crust state before and after the step. This is an independent
  measurement, and its agreement with route 1 is the real closure test.

Declared validity stops (item 9)
--------------------------------
The model stops, rather than extrapolating, when

* the crust loses hydraulic connectivity or clogs (`graded_crust` thresholds);
* the solved interface is more than a declared relative amount below `C_eq`,
  so the Prat reduction's halite-equilibrated base is not defensible;
* the top closure has no positive supersaturation root (growth arrested);
* the surface molality leaves the Pitzer activity correlation's range;
* `Pe << Da` fails, so advection inside the crust is not negligible;
* the soil Newton solve fails.
"""

from __future__ import annotations

import math
import sys
from dataclasses import dataclass, field
from pathlib import Path

_P2A = Path(__file__).resolve().parent.parent.parent / "p2_a_prat_surface"
if str(_P2A) not in sys.path:
    sys.path.insert(0, str(_P2A))

from surface import graded_crust as gc  # noqa: E402
from surface.prat_crust import CrustModelDomainExit  # noqa: E402
from surface.water_activity import MOLALITY_LIMIT  # noqa: E402

from . import soil  # noqa: E402

#: NaCl molar mass used by the activity correlation's molality [kg/mol].
MOLAR_MASS_NACL = 0.05844
#: Largest NaCl mass fraction the Pitzer correlation is declared valid at.
MAXIMUM_MASS_FRACTION = (
    MOLALITY_LIMIT * MOLAR_MASS_NACL / (1.0 + MOLALITY_LIMIT * MOLAR_MASS_NACL)
)


def mass_fraction_molality(mass_fraction: float) -> float:
    """Molality of a NaCl solution of the given mass fraction [mol/kg]."""
    if mass_fraction >= 1.0:
        return float("inf")
    return mass_fraction / ((1.0 - mass_fraction) * MOLAR_MASS_NACL)


class InterfaceDomainExit(RuntimeError):
    """The coupled state left the declared validity envelope."""


@dataclass
class CoupledState:
    column: soil.ColumnState
    crust: gc.GradedCrustState
    time: float = 0.0
    interface_mass_fraction: float = 0.0
    stopped: bool = False
    stop_reason: str | None = None
    seed_salt_kg_m2: float = 0.0
    seed_water_kg_m2: float = 0.0
    cumulative_evaporation_kg_m2: float = 0.0
    cumulative_salt_to_crust_kg_m2: float = 0.0
    cumulative_water_to_crust_kg_m2: float = 0.0
    history: list = field(default_factory=list)


def seed_from_soil(
    params: soil.SoilParameters,
    column: soil.ColumnState,
    crust_params: gc.GradedCrustParameters,
    film,
    seed_height: float,
    cells_per_delta: float = 8.0,
    window_lengths: float = 40.0,
    ratio: float = 1.25,
    maximum_cell_lengths: float = 10.0,
) -> tuple:
    """Create the initial crust seed out of the soil's own top-cell salt.

    Item 6 asks for a *physically consistent* initial seed. A seed conjured
    from nothing would make the coupled salt budget meaningless, and charging
    it to the soil as a transfer (what P2-A and P2-B do) would make the
    initial state depend on the activation event. Here the seed is formed at
    `t = 0` from salt and water that are **removed from the soil's top cells**,
    so the combined initial inventory is exactly the soil's original
    inventory and the coupled budget has a single, unambiguous origin.

    Raises if the soil top cannot supply the seed.
    """
    crust = gc.GradedCrustState(
        top_cell_size=min(
            crust_params.reaction_length() / cells_per_delta, seed_height / 8.0
        )
    )
    delta = crust_params.reaction_length()
    # Activate against a halite-equilibrated base; the first coupled solve
    # immediately replaces that with the solved interface value.
    flux_guess = film.flux_from_mass_fraction(
        crust_params.equilibrium_mass_fraction
    )
    gc.activate(
        crust,
        crust_params,
        seed_height,
        flux_guess,
        0.0,
        window=window_lengths * delta,
        ratio=ratio,
        maximum_cell=maximum_cell_lengths * delta,
    )
    seed_salt = crust.seed_salt
    seed_water = crust.seed_water
    # Undo P2-A/P2-B's "charge it to the soil" bookkeeping: here the seed is
    # an initial condition of the combined system, not a transfer.
    crust.cumulative_salt_from_soil = 0.0
    crust.cumulative_water_from_soil = 0.0

    dz = column.dz
    remaining_salt, remaining_water = seed_salt, seed_water
    for index in range(column.cells - 1, -1, -1):
        total, salt = soil.storage(
            params, column.pressure[index], column.mole_fraction[index]
        )
        cell_salt = salt * dz
        cell_water = (total - salt) * dz
        take_salt = min(remaining_salt, 0.9 * cell_salt)
        take_water = min(remaining_water, 0.9 * cell_water)
        if take_salt <= 0.0 and take_water <= 0.0:
            continue
        new_salt = cell_salt - take_salt
        new_water = cell_water - take_water
        new_w = new_salt / (new_salt + new_water)
        new_x = params.mole_fraction(new_w)
        new_total = (new_salt + new_water) / dz
        # Invert theta*rho = new_total for the saturation, then for pressure.
        sw = new_total / (params.porosity * params.rho(new_x))
        if not (params.residual_saturation < sw <= 1.0):
            raise InterfaceDomainExit(
                "the soil top cannot supply the declared seed without leaving "
                f"the saturation range (would need S_w = {sw!r})"
            )
        column.mole_fraction[index] = new_x
        column.pressure[index] = params.pw(sw)
        remaining_salt -= take_salt
        remaining_water -= take_water
        if remaining_salt <= 1.0e-18 and remaining_water <= 1.0e-18:
            break
    if remaining_salt > 1.0e-15 or remaining_water > 1.0e-15:
        raise InterfaceDomainExit(
            "the soil column does not hold enough salt or water for the "
            f"declared seed (short by {remaining_salt!r} kg/m2 salt, "
            f"{remaining_water!r} kg/m2 water)"
        )
    return crust, seed_salt, seed_water


def crust_response(
    crust: gc.GradedCrustState,
    crust_params: gc.GradedCrustParameters,
    film,
    interface_mass_fraction: float,
) -> dict:
    """Crust-side fluxes for a given interface concentration.

    The gas film is evaluated at the **crust top** mass fraction, which the
    crust solve itself produces, so `j` and `C_top` are closed by bisection
    rather than prescribed.
    """
    def top_for(flux: float) -> float:
        solution = gc.solve_profile(
            crust, crust_params, flux, interface_mass_fraction
        )
        return solution.top_mass_fraction

    def mismatch(flux: float) -> float:
        """``film(C_top(j)) - j``, decreasing in ``j``.

        A crust-top concentration beyond the Pitzer correlation's range is
        reported as "this flux is too large" rather than extrapolated, which
        keeps the bisection monotone without ever evaluating the activity
        outside its declared validity.
        """
        w_top = top_for(flux)
        if mass_fraction_molality(w_top) > MOLALITY_LIMIT:
            return -flux if flux > 0.0 else -1.0
        return film.flux_from_mass_fraction(w_top) - flux

    # j = 0 gives no supersaturation, so film(C_base) > 0 = j and the
    # mismatch is positive; j = film(0) is the free-water maximum, which no
    # saline surface can sustain, so the mismatch is negative there. The
    # bracket is therefore guaranteed and needs no widening.
    low = 0.0
    high = film.flux_from_mass_fraction(0.0)
    for _ in range(200):
        mid = 0.5 * (low + high)
        if mismatch(mid) > 0.0:
            low = mid
        else:
            high = mid
        if high - low <= 1.0e-16 * max(high, 1.0e-30):
            break
    flux = 0.5 * (low + high)
    solution = gc.solve_profile(crust, crust_params, flux, interface_mass_fraction)
    if mass_fraction_molality(solution.top_mass_fraction) > MOLALITY_LIMIT:
        raise InterfaceDomainExit(
            "the crust-top molality at the gas-film fixed point exceeds the "
            f"declared Pitzer limit {MOLALITY_LIMIT!r}; the water activity is "
            "not extrapolated"
        )
    return {
        "evaporation_flux": flux,
        "solution": solution,
        "total_mass_flux": solution.solution_mass_flux_base,
        "salt_mass_flux": solution.salt_flux_base_route,
    }


def soil_interface_salt_flux(
    params: soil.SoilParameters,
    column: soil.ColumnState,
    total_mass_flux: float,
    interface_mass_fraction: float,
) -> float:
    """Soil-side salt flux through the top face, upward positive.

    Uses the same upwinding and the same mass-averaged Fickian closure as
    every interior face, with the interface half a cell above the top cell
    centre. The face diffusivity is one-sided (the top cell's saturation),
    because the soil has no saturation state at the interface; that choice is
    declared and is the only place this module departs from the interior face
    formula.
    """
    x_top = column.mole_fraction[-1]
    w_top = params.mass_fraction(x_top)
    rho_top = params.rho(x_top)
    rho_interface = params.rho(params.mole_fraction(interface_mass_fraction))
    rho_face = 0.5 * (rho_top + rho_interface)
    d_face = params.effective_diffusion_coefficient(params.sw(column.pressure[-1]))
    w_up = w_top if total_mass_flux > 0.0 else interface_mass_fraction
    advective = total_mass_flux * w_up
    diffusive = (
        -rho_face * d_face * (interface_mass_fraction - w_top) / (0.5 * column.dz)
    )
    return advective + diffusive


def step(
    params: soil.SoilParameters,
    state: CoupledState,
    crust_params: gc.GradedCrustParameters,
    film,
    dt: float,
    undersaturation_limit: float = 0.02,
    tolerance: float = 1.0e-9,
    maximum_iterations: int = 40,
) -> dict:
    """One coupled step: solve for the compatible interface, then commit.

    ``tolerance`` is the relative interface-flux residual the Illinois
    iteration must reach. It is two orders tighter than the registered
    acceptance gate and one order above the noise floor set by the soil
    Newton solve, below which further iterations buy nothing.

    Which flux the soil is given (item 8)
    -------------------------------------
    The crust offers two expressions for the NaCl flux through its base:

    * the **transport route** ``w_i rho U(0) + D_base``, reconstructed from
      the concentration profile;
    * the **inventory route** ``d(M_solid + M_dissolved)/dt``, measured from
      the crust state before and after the step.

    They differ by the residual of the audited linearised sink
    (``P2_B_SALT_GROW_PARAMETER_COMPARISON.md`` section 7.2), which is of
    order ``s L_top / F_s``. The soil is handed the **inventory route**,
    because that is the quantity the crust actually accumulates, so the
    coupled pair conserves NaCl exactly rather than to the linearisation
    error. The transport route is then a genuinely independent audit of the
    same flux, and its departure is reported, never absorbed.

    The same argument applies to water: the soil is handed
    ``F_w = j - rho_v v_s + dM_w,c/dt``, so the total solution flux it loses
    is exactly what the crust gains plus what evaporates.
    """
    if state.stopped:
        raise InterfaceDomainExit(f"coupled model stopped: {state.stop_reason}")
    c_eq = crust_params.equilibrium_mass_fraction
    column_before = state.column.copy()
    cache = {}

    def evaluate(w_i: float) -> dict:
        response = crust_response(state.crust, crust_params, film, w_i)
        trial_crust = state.crust.copy()
        crust_record = gc.step(
            trial_crust,
            crust_params,
            response["evaporation_flux"],
            dt,
            base_mass_fraction=w_i,
        )
        salt_inventory = crust_record["salt_flux_soil_kg_m2_s"]
        water_inventory = crust_record["water_flux_soil_kg_m2_s"]
        total = salt_inventory + water_inventory
        trial_column = column_before.copy()
        boundary = soil.BoundaryCondition(
            top_flux=lambda _p, _s, t=total, q=salt_inventory: (t, q),
            bottom="no_flux",
        )
        info = soil.step(params, trial_column, boundary, dt)
        soil_salt = soil_interface_salt_flux(params, trial_column, total, w_i)
        record = {
            "response": response,
            "crust_record": crust_record,
            "trial_crust": trial_crust,
            "trial_column": trial_column,
            "info": info,
            "salt_inventory": salt_inventory,
            "water_inventory": water_inventory,
            "total_mass_flux": total,
            "soil_salt_flux": soil_salt,
            "residual": soil_salt - salt_inventory,
        }
        cache[len(cache)] = w_i
        return record

    # The crust is a very stiff concentration clamp at its base: a 2 %
    # departure from C_eq changes its solution demand by two orders of
    # magnitude and reverses its sign, because the reduced sink acts over a
    # half-micrometre reaction length. The search therefore starts tight and
    # widens only as needed, never past the DECLARED VALIDITY ENVELOPE.
    def try_bracket(offset: float):
        lo = c_eq * (1.0 - offset)
        hi = min(c_eq * (1.0 + offset), 0.99 * MAXIMUM_MASS_FRACTION)
        try:
            r_lo, r_hi = evaluate(lo), evaluate(hi)
        except (soil.ColumnSolverError, soil.ColumnDomainExit,
                CrustModelDomainExit, InterfaceDomainExit):
            return None
        if r_lo["residual"] * r_hi["residual"] > 0.0:
            return None
        return lo, hi, r_lo, r_hi

    offsets = sorted(
        {min(o, undersaturation_limit)
         for o in (1.0e-5, 1.0e-4, 1.0e-3, 1.0e-2, undersaturation_limit)}
    )
    found = None
    offset = offsets[-1]
    for offset in offsets:
        found = try_bracket(offset)
        if found is not None:
            break
    if found is None:
        state.stopped = True
        state.stop_reason = (
            "no compatible interface concentration could be bracketed inside "
            f"the declared envelope C_eq*(1 +/- {undersaturation_limit!r}); "
            "the reduced Prat closure is not defensible here"
        )
        raise InterfaceDomainExit(state.stop_reason)

    # Illinois (safeguarded regula falsi): keeps the bracket, converges
    # superlinearly, so the interface residual reaches round-off in a dozen
    # evaluations instead of the forty bisection would need.
    low, high, record_low, record_high = found
    g_low, g_high = record_low["residual"], record_high["residual"]
    record = record_low if abs(g_low) < abs(g_high) else record_high
    scale_guess = max(abs(record["salt_inventory"]), 1.0e-300)
    for _ in range(maximum_iterations):
        if g_high == g_low:
            break
        mid = high - g_high * (high - low) / (g_high - g_low)
        if not (min(low, high) < mid < max(low, high)):
            mid = 0.5 * (low + high)
        record = evaluate(mid)
        g_mid = record["residual"]
        if abs(g_mid) <= tolerance * scale_guess:
            break
        if g_mid * g_low < 0.0:
            high, g_high = mid, g_mid
            g_low *= 0.5
        else:
            low, g_low = mid, g_mid
            g_high *= 0.5
        if abs(high - low) <= 1.0e-14 * max(abs(high), 1.0e-30):
            break
    w_i = mid if "mid" in dir() else 0.5 * (low + high)
    response = record["response"]
    solution = response["solution"]

    # --- declared validity gates, before anything is committed -------------
    if solution.growth_arrested:
        state.stopped = True
        state.stop_reason = "growth_arrested_no_positive_supersaturation"
        raise InterfaceDomainExit(state.stop_reason)
    if w_i < c_eq * (1.0 - 0.999 * undersaturation_limit):
        state.stopped = True
        state.stop_reason = (
            f"interface undersaturated: w_i = {w_i!r} has reached the lower "
            f"edge of the declared envelope below C_eq = {c_eq!r}, so the "
            "Prat reduction's halite-equilibrated base is not defensible"
        )
        raise InterfaceDomainExit(state.stop_reason)
    peclet = (
        response["evaporation_flux"]
        / crust_params.rho
        * state.crust.height
        / (crust_params.deposition_porosity
           * crust_params.diffusivity(crust_params.deposition_porosity))
    )
    damkohler = (state.crust.height / solution.reaction_length) ** 2
    if peclet >= 0.01 * damkohler:
        state.stopped = True
        state.stop_reason = (
            f"Pe = {peclet!r} is not small against Da = {damkohler!r}; the "
            "audited Pe << Da reduction is invalid"
        )
        raise InterfaceDomainExit(state.stop_reason)

    # --- commit the states the converged fluxes actually produced ----------
    state.column = record["trial_column"]
    state.crust = record["trial_crust"]
    crust_record = record["crust_record"]
    state.interface_mass_fraction = w_i
    state.time += dt
    state.cumulative_evaporation_kg_m2 += response["evaporation_flux"] * dt
    state.cumulative_salt_to_crust_kg_m2 += record["salt_inventory"] * dt
    state.cumulative_water_to_crust_kg_m2 += record["water_inventory"] * dt
    if state.crust.stopped:
        state.stopped = True
        state.stop_reason = state.crust.stop_reason

    scale = max(abs(record["salt_inventory"]), 1.0e-300)
    transport_route = response["salt_mass_flux"]
    out = {
        "time_s": state.time,
        "dt_s": dt,
        "interface_mass_fraction": w_i,
        "interface_over_c_eq": w_i / c_eq,
        "interface_residual_kg_m2_s": record["residual"],
        "interface_residual_relative": abs(record["residual"]) / scale,
        "interface_iterations": len(cache),
        "interface_bracket_offset": offset,
        "evaporation_flux_kg_m2_s": response["evaporation_flux"],
        "crust_top_mass_fraction": solution.top_mass_fraction,
        "crust_top_supersaturation": solution.top_supersaturation,
        "growth_velocity_m_s": solution.growth_velocity,
        "crust_height_m": state.crust.height,
        "crust_minimum_porosity": crust_record["minimum_porosity"],
        "crust_cells": crust_record["cells"],
        "internal_reaction_kg_m2_s": solution.internal_reaction,
        "route1_crust_transport_kg_m2_s": transport_route,
        "route2_crust_inventory_kg_m2_s": record["salt_inventory"],
        "route3_soil_component_kg_m2_s": record["soil_salt_flux"],
        "route_2_minus_3_relative": abs(
            record["salt_inventory"] - record["soil_salt_flux"]
        ) / scale,
        "route_1_minus_2_relative": abs(transport_route - record["salt_inventory"])
        / scale,
        "route_1_minus_2_predicted_relative": abs(
            solution.top_supersaturation
            * solution.solution_mass_flux_top_relative
        ) / scale,
        "total_mass_flux_kg_m2_s": record["total_mass_flux"],
        "water_flux_kg_m2_s": record["water_inventory"],
        "soil_top_mass_fraction": params.mass_fraction(
            state.column.mole_fraction[-1]
        ),
        "soil_top_saturation": params.sw(state.column.pressure[-1]),
        "soil_top_pressure_pa": state.column.pressure[-1],
        "peclet": peclet,
        "damkohler": damkohler,
        "peclet_over_damkohler": peclet / damkohler,
        "newton_iterations": record["info"]["iterations"],
        "stopped": state.stopped,
        "stop_reason": state.stop_reason,
    }
    state.history.append(out)
    return out


def audit(
    params: soil.SoilParameters,
    state: CoupledState,
    initial_column: soil.ColumnState,
    seed_salt: float,
    seed_water: float,
) -> dict:
    """Combined soil + crust component budget.

    The seed is an initial condition formed from the soil, so the combined
    initial inventory is the soil's pre-seed inventory. Nothing enters or
    leaves the pair except water vapour at the crust top.
    """
    soil_now = soil.inventory(params, state.column)
    soil_before = soil.inventory(params, initial_column)
    crust_salt = state.crust.solid_salt_mass + state.crust.dissolved_salt
    crust_water = state.crust.water

    salt_now = soil_now["salt_kg_m2"] + crust_salt
    salt_before = soil_before["salt_kg_m2"]
    water_now = soil_now["water_kg_m2"] + crust_water
    water_before = soil_before["water_kg_m2"]
    evaporated = state.cumulative_evaporation_kg_m2

    salt_residual = salt_now - salt_before
    water_residual = water_now - water_before + evaporated
    return {
        "soil_salt_kg_m2": soil_now["salt_kg_m2"],
        "crust_salt_kg_m2": crust_salt,
        "total_salt_kg_m2": salt_now,
        "initial_total_salt_kg_m2": salt_before,
        "salt_residual_kg_m2": salt_residual,
        "salt_residual_relative": abs(salt_residual) / max(abs(salt_before), 1e-300),
        "soil_water_kg_m2": soil_now["water_kg_m2"],
        "crust_water_kg_m2": crust_water,
        "total_water_kg_m2": water_now,
        "initial_total_water_kg_m2": water_before,
        "evaporated_kg_m2": evaporated,
        "water_residual_kg_m2": water_residual,
        "water_residual_relative": abs(water_residual)
        / max(abs(water_before), abs(evaporated), 1e-300),
        "seed_salt_kg_m2": seed_salt,
        "seed_water_kg_m2": seed_water,
        "transferred_salt_to_crust_kg_m2": state.cumulative_salt_to_crust_kg_m2,
        "transferred_water_to_crust_kg_m2": state.cumulative_water_to_crust_kg_m2,
    }
