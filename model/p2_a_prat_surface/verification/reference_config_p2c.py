"""The declared P2-C configuration: salt_grow parameters in the P2-B model.

Every number here comes either from `config/acceptance_criteria.json` (the
P2-A declared reference), from `config/acceptance_criteria_p2c.json` (the
parameters transferred from salt_grow, registered before any run), or from the
salt_grow baseline config whose path and values are recorded there. Nothing is
fitted and nothing is calibrated against an observed growth rate.
"""

from __future__ import annotations

import json
import math

from reference_config import HERE, PARAMS, reference_film  # noqa: F401
from surface import graded_crust as gc

CRITERIA_P2C_PATH = HERE / "config" / "acceptance_criteria_p2c.json"

with CRITERIA_P2C_PATH.open("r", encoding="utf-8") as _handle:
    CRITERIA_P2C = json.load(_handle)

TRANSFER = {
    k: v["value"] for k, v in CRITERIA_P2C["transferred_parameters"].items()
}
TOL_P2C = CRITERIA_P2C["tolerances"]
REFERENCE = CRITERIA_P2C["reference_project"]

C_EQ = PARAMS["equilibrium_mass_fraction"]
H_SEED = PARAMS["seed_height"]

#: The declared evaporation flux of the P2-A gas film at halite equilibrium.
#: Identical to the one P2-B used, so Case 0 is the committed P2-B setting.
FILM_FLUX = reference_film().flux_from_mass_fraction(C_EQ)

#: salt_grow baseline, read from configs/baseline.json.
SALT_GROW = {
    "C_eq": 0.264,
    "rho_solution": 1200.0,
    "rho_crystal": 2170.0,
    "epsilon_ref": 0.43,
    "epsilon_birth": 0.43,
    "D_molecular": 1.3e-9,
    "k_reaction": 2.3e-3,
    "beta": 2.0,
    "d_sub": 25e-6,
    "h0": 1e-3,
    "U_evap": 1e-2 / 86400.0,
}
SALT_GROW["a_v_ref"] = 36.0 * (1.0 - SALT_GROW["epsilon_ref"]) / SALT_GROW["d_sub"]
SALT_GROW["D_eff_ref"] = math.sqrt(SALT_GROW["epsilon_ref"]) * SALT_GROW["D_molecular"]
SALT_GROW["evaporation_mass_flux"] = (
    SALT_GROW["rho_solution"] * SALT_GROW["U_evap"]
)


# --- the one-factor case ladder ----------------------------------------------

#: Case 0 is the committed P2-B configuration. Cases 1-3 change exactly one
#: quantity each, cumulatively, and nothing else.
CASES = {
    0: {"label": "P2-B original", "overrides": {}},
    1: {"label": "+ eps_0 = 0.43", "overrides": {"deposition_porosity": 0.43}},
    2: {
        "label": "+ salt_grow k_r",
        "overrides": {
            "deposition_porosity": 0.43,
            "reaction_coefficient": TRANSFER["reaction_coefficient"],
        },
    },
    3: {
        "label": "+ salt_grow a_v",
        "overrides": {
            "deposition_porosity": 0.43,
            "reaction_coefficient": TRANSFER["reaction_coefficient"],
            "specific_area": TRANSFER["specific_area"],
        },
    },
}


def case_parameters(case: int, **extra) -> gc.GradedCrustParameters:
    """Case 0-3 of the controlled one-factor comparison."""
    kwargs = dict(
        equilibrium_mass_fraction=C_EQ,
        deposition_porosity=PARAMS["crust_porosity"],
        crystal_density=PARAMS["crystal_density"],
        specific_area=PARAMS["specific_area"],
        reaction_coefficient=PARAMS["reaction_coefficient"],
        effective_salt_diffusivity=PARAMS["effective_salt_diffusivity"],
        beta=PARAMS["beta"],
        vapour_density=PARAMS["vapour_density"],
        clogging_porosity=0.01,
        connectivity_porosity=0.02,
    )
    kwargs.update(CASES[case]["overrides"])
    kwargs.update(extra)
    return gc.GradedCrustParameters(**kwargs)


# --- the fully matched salt_grow benchmark ------------------------------------


def matched_parameters(closure: str = "B2", **extra) -> gc.GradedCrustParameters:
    """Every physical input matched to salt_grow, including its closures.

    ``B1`` freezes ``a_v`` and ``D_eff`` at their reference-porosity values;
    ``B2`` uses ``a_v = 36(1-eps)/d_sub`` and ``D_eff = sqrt(eps) D_mol``.
    The vapour-displacement term of P2-A has no counterpart in salt_grow and
    is set to zero here so the comparison is like for like.
    """
    if closure not in ("B1", "B2"):
        raise ValueError(f"unknown closure {closure!r}")
    kwargs = dict(
        equilibrium_mass_fraction=SALT_GROW["C_eq"],
        deposition_porosity=SALT_GROW["epsilon_birth"],
        crystal_density=SALT_GROW["rho_crystal"],
        specific_area=SALT_GROW["a_v_ref"],
        reaction_coefficient=SALT_GROW["k_reaction"],
        effective_salt_diffusivity=SALT_GROW["D_eff_ref"],
        beta=SALT_GROW["beta"],
        solution_density=SALT_GROW["rho_solution"],
        vapour_density=0.0,
        clogging_porosity=0.01,
        connectivity_porosity=0.02,
    )
    if closure == "B2":
        kwargs["area_closure"] = "salt_grow_b2"
        kwargs["substrate_size"] = SALT_GROW["d_sub"]
        kwargs["diffusivity_closure"] = "sqrt_eps"
        kwargs["molecular_diffusivity"] = SALT_GROW["D_molecular"]
    kwargs.update(extra)
    return gc.GradedCrustParameters(**kwargs)


# --- analytic diagnostics -----------------------------------------------------


def dimensionless_groups(
    params: gc.GradedCrustParameters,
    height: float,
    evaporation_flux: float,
) -> dict:
    """Pe, Da, Da_top, Da*, delta_r and the closed-form Prat quantities.

    Evaluated at the deposition porosity, i.e. at the state of the crust the
    instant it is created, which is the only state at which the closed-form
    law is defined.
    """
    eps = params.deposition_porosity
    one = 1.0 - params.equilibrium_mass_fraction
    rho = params.rho
    d = params.diffusivity(eps)
    a_v = params.area(eps)
    k_r = params.reaction_coefficient
    beta = params.beta
    pe = (evaporation_flux / rho) * height / (eps * d)
    da = height * height * one * a_v * k_r / (eps * d)
    da_top = beta * (1.0 - eps) * one * height * k_r / (eps * d)
    da_star = beta**2 * (1.0 - eps) ** 2 * one * k_r / (eps * d * a_v)
    delta = params.reaction_length(eps)
    root_da = math.sqrt(da)
    tanh = math.tanh(root_da) if root_da < 350.0 else 1.0
    factor = math.sqrt(da_star) * tanh
    c_star_top = (
        params.equilibrium_mass_fraction * pe * tanh / (root_da * (1.0 + factor))
    )
    v_s = beta * k_r * rho * c_star_top / params.crystal_density
    efficiency = factor / (1.0 + factor)
    return {
        "Pe": pe,
        "Da": da,
        "Da_top": da_top,
        "Da_star": da_star,
        "sqrt_Da_star": math.sqrt(da_star),
        "delta_r_m": delta,
        "h_over_delta_r": height / delta,
        "c_star_top": c_star_top,
        "top_mass_fraction": params.equilibrium_mass_fraction + c_star_top,
        "growth_velocity_m_s": v_s,
        "growth_velocity_mm_day": v_s * 86400.0 * 1e3,
        "efficiency_f": efficiency,
        "internal_over_top": (1.0 - efficiency) / efficiency,
        "tau_porosity_s": (
            eps * params.crystal_density / (a_v * rho * k_r * c_star_top)
        ),
        "tau_burial_s": delta / v_s,
        "tau_growth_s": height / v_s,
        "predicted_wake_drop": params.predicted_wake_porosity_drop(),
        "predicted_wake_porosity": eps - params.predicted_wake_porosity_drop(),
        "predicted_sealing_growth_m": params.predicted_sealing_growth(),
        "seals_predicted": math.sqrt(da_star)
        <= (1.0 - eps) / (eps - params.connectivity_porosity),
    }


# --- run control --------------------------------------------------------------


def seeded_state(
    params: gc.GradedCrustParameters,
    evaporation_flux: float,
    height: float,
    cells_per_delta: float = 8.0,
    minimum_cells: int = 32,
    window_lengths: float = 40.0,
    ratio: float = 1.25,
    maximum_cell_lengths: float = 10.0,
) -> gc.GradedCrustState:
    """A freshly activated crust on a top-refined grid.

    The top cell is ``delta_r / cells_per_delta``, but never coarser than
    ``height / minimum_cells``, so a crust whose reaction length exceeds its
    own height (the P2-A regime) still gets a usable grid.
    """
    delta = params.reaction_length()
    top_cell = min(delta / cells_per_delta, height / minimum_cells)
    state = gc.GradedCrustState(top_cell_size=top_cell)
    gc.activate(
        state,
        params,
        height,
        evaporation_flux,
        0.0,
        window=window_lengths * delta,
        ratio=ratio,
        maximum_cell=maximum_cell_lengths * delta,
    )
    return state


def run_case(
    params: gc.GradedCrustParameters,
    evaporation_flux: float,
    height: float,
    cells_per_delta: float = 8.0,
    steps_per_tau: float = 1500.0,
    maximum_steps: int = 40000,
    plateau_windows: float = 5.0,
    velocity_drift: float | None = None,
    porosity_drift: float = 1.0e-5,
    horizon: float | None = None,
    snapshot_times: list | None = None,
    on_snapshot=None,
    freeze_lengths: float = 0.0,
    maximum_cell_lengths: float = 10.0,
    time_step: float | None = None,
    **seed_kwargs,
):
    """Integrate until the crust stops, plateaus, or hits the declared cap.

    The time step is ``tau_eps / steps_per_tau`` with
    ``tau_eps = eps_0 rho_cr / (a_v rho k_r C*_top)``, i.e. the time for the
    top porosity to be consumed at its initial rate. That estimate only sets a
    step size; it is never used as a result. A plateau is declared when the
    growth velocity and the minimum porosity both stop drifting over
    ``plateau_windows`` characteristic times; the declaration is recorded so
    it can be checked against the trajectory.
    """
    drift = TOL_P2C["velocity_plateau_drift_relative"] if velocity_drift is None else velocity_drift
    state = seeded_state(
        params,
        evaporation_flux,
        height,
        cells_per_delta=cells_per_delta,
        maximum_cell_lengths=maximum_cell_lengths,
        **seed_kwargs,
    )
    initial = state.copy()
    first = gc.solve_profile(state, params, evaporation_flux)
    delta = params.reaction_length()
    tau_eps = (
        params.deposition_porosity
        * params.crystal_density
        / (
            params.area(params.deposition_porosity)
            * params.rho
            * params.reaction_coefficient
            * max(first.top_supersaturation, 1.0e-300)
        )
    )
    tau_burial = delta / max(first.growth_velocity, 1.0e-300)
    scale = max(tau_eps, tau_burial)
    dt = tau_eps / steps_per_tau if time_step is None else time_step
    cap = horizon if horizon is not None else 60.0 * scale
    freeze = freeze_lengths * delta
    maximum_cell = maximum_cell_lengths * delta if freeze > 0.0 else 0.0
    window_steps = max(1, int(plateau_windows * scale / dt))
    pending = sorted(snapshot_times) if snapshot_times else []
    reason = "step_cap"
    for index in range(maximum_steps):
        record = gc.step(
            state,
            params,
            evaporation_flux,
            dt,
            freeze_depth=freeze,
            maximum_cell=maximum_cell,
        )
        while pending and record["time_s"] >= pending[0]:
            if on_snapshot is not None:
                on_snapshot(pending[0], state)
            pending.pop(0)
        if state.stopped:
            reason = state.stop_reason
            break
        if record["time_s"] >= cap:
            reason = "horizon"
            break
        if index >= 2 * window_steps:
            past = state.history[-1 - window_steps]
            now = record
            dv = abs(now["growth_velocity_m_s"] - past["growth_velocity_m_s"]) / max(
                abs(now["growth_velocity_m_s"]), 1.0e-300
            )
            de = abs(now["minimum_porosity"] - past["minimum_porosity"])
            if dv < drift and de < porosity_drift:
                reason = "plateau"
                break
    return {
        "state": state,
        "initial": initial,
        "first": first,
        "dt_s": dt,
        "tau_porosity_s": tau_eps,
        "tau_burial_s": tau_burial,
        "stop_reason": reason,
        "steps": len(state.history),
        "delta_r_m": delta,
        "top_cell_m": state.top_cell_size,
        "cells_per_delta_actual": delta / state.top_cell_size,
    }
