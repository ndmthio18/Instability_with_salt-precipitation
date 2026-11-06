"""The declared P2-B reference configuration (evolving crust porosity).

Reuses every P2-A declared value and adds only what the depth-resolved model
needs. Nothing is fitted.
"""

from __future__ import annotations

import json

from reference_config import HERE, PARAMS, reference_film  # noqa: F401
from surface import evolving_crust as ec

CRITERIA_P2B_PATH = HERE / "config" / "acceptance_criteria_p2b.json"

with CRITERIA_P2B_PATH.open("r", encoding="utf-8") as _handle:
    CRITERIA_P2B = json.load(_handle)

P2B = {k: v["value"] for k, v in CRITERIA_P2B["declared_reference_parameters"].items()}
TOL_P2B = CRITERIA_P2B["tolerances"]


def evolving_parameters(**overrides) -> ec.EvolvingCrustParameters:
    kwargs = dict(
        equilibrium_mass_fraction=PARAMS["equilibrium_mass_fraction"],
        deposition_porosity=P2B["deposition_porosity"],
        crystal_density=PARAMS["crystal_density"],
        specific_area=PARAMS["specific_area"],
        reaction_coefficient=PARAMS["reaction_coefficient"],
        effective_salt_diffusivity=PARAMS["effective_salt_diffusivity"],
        beta=PARAMS["beta"],
        vapour_density=PARAMS["vapour_density"],
        clogging_porosity=P2B["clogging_porosity"],
        connectivity_porosity=P2B["connectivity_porosity"],
        diffusivity_exponent=P2B["diffusivity_exponent"],
        specific_area_exponent=P2B["specific_area_exponent"],
    )
    kwargs.update(overrides)
    return ec.EvolvingCrustParameters(**kwargs)


def seeded_state(
    params: ec.EvolvingCrustParameters,
    evaporation_flux: float,
    cells: int | None = None,
    seed_height: float | None = None,
) -> ec.EvolvingCrustState:
    """A freshly activated crust at the declared seed height."""
    h0 = PARAMS["seed_height"] if seed_height is None else seed_height
    n = int(P2B["cells_per_seed_height"]) if cells is None else cells
    state = ec.EvolvingCrustState(cell_size=h0 / n)
    ec.activate(state, params, h0, evaporation_flux, 0.0)
    return state


def snapshot(state: ec.EvolvingCrustState) -> ec.EvolvingCrustState:
    """A deep-enough copy to take a budget against."""
    clone = ec.EvolvingCrustState(cell_size=state.cell_size)
    for key, value in state.__dict__.items():
        setattr(clone, key, list(value) if isinstance(value, list) else value)
    return clone


def run_to_sealing(
    params: ec.EvolvingCrustParameters,
    evaporation_flux: float,
    dt: float | None = None,
    cells: int | None = None,
    seed_height: float | None = None,
    maximum_steps: int = 60000,
    steps_target: int = 1500,
) -> tuple[ec.EvolvingCrustState, ec.EvolvingCrustState, float]:
    """Integrate until the crust seals, or until `maximum_steps`.

    When `dt` is omitted it is set from the analytic sealing estimate
    ``t_seal ~ eps_0 rho_cr / (a_v rho k_r C*_top)`` so that every parameter
    combination is resolved with a comparable number of steps. That estimate
    is only used to choose a step size; it is never used as a result.
    """
    state = seeded_state(params, evaporation_flux, cells=cells, seed_height=seed_height)
    initial = snapshot(state)
    if dt is None:
        first = ec.solve_profile(state, params, evaporation_flux)
        estimate = (
            params.deposition_porosity
            * params.crystal_density
            / (
                params.specific_area
                * params.rho
                * params.reaction_coefficient
                * max(first.top_supersaturation, 1.0e-30)
            )
        )
        dt = estimate / steps_target
    for _ in range(maximum_steps):
        ec.step(state, params, evaporation_flux, dt)
        if state.sealed:
            break
    return state, initial, dt
