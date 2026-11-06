#!/usr/bin/env python3
"""Generate the P2-B evolving-porosity evidence.

Writes CSV tables and a JSON summary into
``results/p2b_evolving_porosity/``. Every parameter comes from
``config/acceptance_criteria_p2b.json`` and the P2-A reference; nothing is
fitted.

Usage::

    python run_p2b.py [output_directory]
"""

from __future__ import annotations

import csv
import json
import math
import platform
import sys
from datetime import date
from pathlib import Path

from reference_config import PARAMS, reference_crust, reference_film
from reference_config_p2b import (
    CRITERIA_P2B,
    P2B,
    evolving_parameters,
    run_to_sealing,
    seeded_state,
    snapshot,
)
from surface import evolving_crust as ec
from surface import prat_crust as pc

C_EQ = PARAMS["equilibrium_mass_fraction"]
H_SEED = PARAMS["seed_height"]
MODULE_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_STRIDE = 20


def write_csv(path: Path, rows: list) -> None:
    if not rows:
        raise ValueError(f"refusing to write an empty table: {path}")
    fieldnames: list = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, restval="")
        writer.writeheader()
        writer.writerows(rows)


def flux(relative_humidity: float | None = None) -> float:
    film = (
        reference_film()
        if relative_humidity is None
        else reference_film(relative_humidity=relative_humidity)
    )
    return film.flux_from_mass_fraction(C_EQ)


# --- 1. reduction to the audited closed form ---------------------------------


def table_eq28_reduction() -> list:
    j = flux()
    frozen = evolving_parameters(vapour_density=0.0).frozen()
    reference = pc.growth_velocity(frozen, H_SEED, j)
    rows = []
    for linearised in (True, False):
        previous = None
        for cells in (8, 16, 32, 64, 128, 256, 512):
            params = evolving_parameters(
                vapour_density=0.0, linearised_top_closure=linearised
            )
            state = seeded_state(params, j, cells=cells)
            solution = ec.solve_profile(state, params, j)
            error = abs(solution.growth_velocity - reference) / reference
            order = (
                math.log2(previous / error)
                if previous is not None and error > 0.0
                else float("nan")
            )
            rows.append(
                {
                    "linearised_top_closure": linearised,
                    "cells_per_seed_height": cells,
                    "growth_velocity_m_s": solution.growth_velocity,
                    "equation_28_velocity_m_s": reference,
                    "relative_difference": error,
                    "observed_order": order,
                    "top_supersaturation": solution.top_supersaturation,
                    "predicted_offset_c_star_over_c_eq": (
                        pc.top_supersaturation(frozen, H_SEED, j) / C_EQ
                    ),
                }
            )
            previous = error
    return rows


# --- 2. reference trajectory and profiles ------------------------------------


def table_reference_trajectory() -> tuple[list, list, dict]:
    params = evolving_parameters()
    j = flux()
    state = seeded_state(params, j)
    initial = snapshot(state)
    dt = 100.0
    snapshots: list = []
    snapshot_times = set()
    rows: list = []
    index = 0
    while not state.sealed and index < 20000:
        record = ec.step(state, params, j, dt)
        if index % OUTPUT_STRIDE == 0 or state.sealed:
            row = {"output_stride": OUTPUT_STRIDE}
            row.update(record)
            rows.append(row)
        index += 1
    final_time = state.time
    for fraction in (0.0, 0.25, 0.5, 0.75, 1.0):
        snapshot_times.add(round(fraction * final_time, 6))

    # Re-run to capture depth profiles at the requested fractions.
    replay = seeded_state(params, j)
    targets = sorted(snapshot_times)
    cursor = 0
    solution = ec.solve_profile(replay, params, j)
    if targets and targets[0] <= 0.0:
        snapshots.extend(_profile_rows(replay, solution, 0.0))
        cursor = 1
    for _ in range(20000):
        if replay.sealed or cursor >= len(targets):
            break
        ec.step(replay, params, j, dt)
        if replay.time >= targets[cursor] - 1.0e-9:
            snapshots.extend(
                _profile_rows(
                    replay, ec.solve_profile(replay, params, j), replay.time
                )
            )
            cursor += 1
    if not replay.sealed:
        snapshots.extend(
            _profile_rows(replay, ec.solve_profile(replay, params, j), replay.time)
        )

    budget = ec.budget(state, initial)
    last = state.history[-1]
    summary = {
        "dt_s": dt,
        "cells_per_seed_height": P2B["cells_per_seed_height"],
        "sealed": state.sealed,
        "sealing_time_s": state.sealed_time,
        "seed_height_m": H_SEED,
        "final_height_m": state.height,
        "growth_before_sealing_m": state.height - H_SEED,
        "predicted_growth_m": params.beta * params.deposition_porosity
        / params.specific_area,
        "growth_ratio_to_prediction": (state.height - H_SEED)
        / (params.beta * params.deposition_porosity / params.specific_area),
        "relative_height_increase": (state.height - H_SEED) / H_SEED,
        "initial_growth_velocity_m_s": state.history[0]["growth_velocity_m_s"],
        "final_growth_velocity_m_s": last["growth_velocity_m_s"],
        "growth_velocity_amplification": last["growth_velocity_m_s"]
        / state.history[0]["growth_velocity_m_s"],
        "final_minimum_porosity": last["minimum_porosity"],
        "final_base_porosity": last["base_porosity"],
        "final_top_porosity": last["top_porosity"],
        "final_mean_porosity": last["mean_porosity"],
        "cumulative_internal_over_top": (
            budget["internal_reaction_kg_m2"] / budget["top_skeleton_kg_m2"]
        ),
        "salt_residual_relative": budget["salt_residual_relative"],
        "water_residual_relative": budget["water_residual_relative"],
        "solid_split_residual_relative": budget["solid_split_residual_relative"],
        "seed_salt_kg_m2": state.seed_salt,
        "seed_water_kg_m2": state.seed_water,
        "seed_note": (
            "The budget is taken from the post-activation state, so the "
            "seed appears as an absolute inventory here rather than as a "
            "delta. It is charged to the soil at activation."
        ),
    }
    return rows, snapshots, summary


def _profile_rows(state, solution, time: float) -> list:
    rows = []
    for position, centre in enumerate(solution.centres):
        rows.append(
            {
                "time_s": time,
                "height_m": state.height,
                "z_m": centre,
                "z_over_seed_height": centre / H_SEED,
                "porosity": solution.porosity[position],
                "supersaturation": solution.profile[position],
                "mass_fraction": C_EQ + solution.profile[position],
                "thickness_m": solution.thicknesses[position],
            }
        )
    return rows


# --- 3. grid and time-step convergence ---------------------------------------


def table_grid_convergence() -> list:
    params = evolving_parameters()
    j = flux()
    rows = []
    for cells in (16, 32, 64, 128):
        for dt in (400.0, 200.0, 100.0):
            state, _, _ = run_to_sealing(params, j, dt=dt, cells=cells)
            profile = state.porosity_profile()
            index = min(range(len(profile)), key=lambda i: profile[i])
            rows.append(
                {
                    "cells_per_seed_height": cells,
                    "dt_s": dt,
                    "sealing_time_s": state.sealed_time,
                    "final_height_m": state.height,
                    "growth_m": state.height - H_SEED,
                    "clogging_depth_m": state.centre(index),
                    "clogging_depth_over_seed_height": state.centre(index) / H_SEED,
                }
            )
    return rows


# --- 4. the growth-before-sealing law ----------------------------------------


def table_sealing_law() -> list:
    rows = []
    for a_v in (1.0e4, 1.0e5, 1.0e6, 1.0e7):
        for k_r in (1.0e-8, 1.0e-7, 1.0e-6):
            for humidity in (0.2, 0.4, 0.6):
                for porosity in (0.3, 0.4, 0.6):
                    params = evolving_parameters(
                        specific_area=a_v,
                        reaction_coefficient=k_r,
                        deposition_porosity=porosity,
                    )
                    j = flux(humidity)
                    state, _, dt = run_to_sealing(params, j)
                    growth = state.height - H_SEED
                    predicted = params.beta * porosity / a_v
                    # The reaction boundary layer, not the growth increment,
                    # is what the grid has to resolve: C* decays over
                    # h/sqrt(Da) below the top.
                    damkohler = params.frozen().damkohler(H_SEED)
                    layer = H_SEED / math.sqrt(damkohler)
                    cell = H_SEED / P2B["cells_per_seed_height"]
                    rows.append(
                        {
                            "specific_area_1_m": a_v,
                            "reaction_coefficient_m_s": k_r,
                            "relative_humidity": humidity,
                            "deposition_porosity": porosity,
                            "evaporation_flux_kg_m2_s": j,
                            "dt_s": dt,
                            "sealed": state.sealed,
                            "sealing_time_s": state.sealed_time,
                            "growth_m": growth,
                            "predicted_growth_m": predicted,
                            "ratio_to_prediction": growth / predicted,
                            "predicted_over_seed_height": predicted / H_SEED,
                            "damkohler_at_seed": damkohler,
                            "boundary_layer_m": layer,
                            "boundary_layer_cells": layer / cell,
                            "resolved": layer / cell >= 4.0,
                            "thin_crust": predicted / H_SEED <= 0.2,
                            "crystal_size_cubic_m": 6.0
                            * (1.0 - porosity)
                            / a_v,
                            "growth_over_crystal_size": growth
                            / (6.0 * (1.0 - porosity) / a_v),
                        }
                    )
    return rows


# --- 5. declared closure sensitivities ---------------------------------------


def predicted_growth(params, seal_porosity: float) -> float:
    """Closed form for the a_v(eps) power-law family.

    With ``a_v(eps) = a_v0 (eps/eps_0)**m`` the two rates

        d eps/dt = -(a_v0/eps_0**m) eps**m rho k_r C* / rho_cr
        d h  /dt =  beta k_r rho C* / rho_cr

    divide to give ``dh/d eps = -beta eps_0**m / (a_v0 eps**m)``, in which
    ``k_r``, ``rho``, ``rho_cr`` and ``C*`` have all cancelled. Integrating
    from ``eps_0`` down to the sealing porosity,

        dh = (beta eps_0**m / a_v0) * integral_{eps_seal}^{eps_0} eps**-m d eps

        m = 0:  (beta/a_v0)(eps_0 - eps_seal)
        m = 1:  (beta eps_0/a_v0) ln(eps_0/eps_seal)
        m = 2:  (beta eps_0/a_v0)(eps_0/eps_seal - 1)

    The growth before sealing is therefore independent of the evaporation
    flux, the reaction coefficient and the supersaturation for **any** power
    law, not only the constant case. It is exact while the clogging point
    stays at the top (``dh << h``), and a lower bound otherwise, because
    deeper material sees a smaller ``C*`` and takes longer to seal.
    """
    m = params.specific_area_exponent
    eps0 = params.deposition_porosity
    scale = params.beta * eps0**m / params.specific_area
    if m == 1.0:
        integral = math.log(eps0 / seal_porosity)
    else:
        integral = (eps0 ** (1.0 - m) - seal_porosity ** (1.0 - m)) / (1.0 - m)
    return scale * integral


def table_closure_sensitivity() -> list:
    j = flux()
    rows = []
    for area_exponent in (0.0, 1.0, 2.0):
        for diffusivity_exponent in (0.0, 1.5):
            params = evolving_parameters(
                specific_area_exponent=area_exponent,
                diffusivity_exponent=diffusivity_exponent,
            )
            state, initial, dt = run_to_sealing(
                params, j, dt=200.0, maximum_steps=20000
            )
            budget = ec.budget(state, initial)
            last = state.history[-1]
            rows.append(
                {
                    "specific_area_exponent": area_exponent,
                    "diffusivity_exponent": diffusivity_exponent,
                    "dt_s": dt,
                    "steps": len(state.history),
                    "sealed": state.sealed,
                    "sealing_time_s": state.sealed_time,
                    "elapsed_time_s": state.time,
                    "growth_m": state.height - H_SEED,
                    "predicted_growth_m": predicted_growth(
                        params, params.clogging_porosity
                    ),
                    "ratio_to_prediction": (state.height - H_SEED)
                    / predicted_growth(params, params.clogging_porosity),
                    "growth_relative_to_constant_closure": float("nan"),
                    "minimum_porosity": last["minimum_porosity"],
                    "growth_velocity_m_s": last["growth_velocity_m_s"],
                    "top_supersaturation": last["top_supersaturation"],
                    "salt_residual_relative": budget["salt_residual_relative"],
                    "water_residual_relative": budget["water_residual_relative"],
                }
            )
    baseline = rows[0]["growth_m"]
    for row in rows:
        row["growth_relative_to_constant_closure"] = row["growth_m"] / baseline
    return rows


# --- 6. frozen versus evolving porosity --------------------------------------


def table_frozen_versus_evolving() -> tuple[list, dict]:
    params = evolving_parameters()
    frozen = reference_crust()
    j = flux()
    state = seeded_state(params, j)
    dt = 100.0
    frozen_height = H_SEED
    rows = []
    index = 0
    worst = 0.0
    while not state.sealed and index < 20000:
        record = ec.step(state, params, j, dt)
        frozen_velocity = pc.growth_velocity(frozen, frozen_height, j)
        frozen_height += frozen_velocity * dt
        difference = abs(state.height - frozen_height) / (state.height - H_SEED)
        worst = max(worst, difference)
        if index % OUTPUT_STRIDE == 0 or state.sealed:
            rows.append(
                {
                    "output_stride": OUTPUT_STRIDE,
                    "time_s": record["time_s"],
                    "height_evolving_m": state.height,
                    "height_frozen_m": frozen_height,
                    "growth_evolving_m": state.height - H_SEED,
                    "growth_frozen_m": frozen_height - H_SEED,
                    "growth_ratio_evolving_over_frozen": (
                        (state.height - H_SEED) / (frozen_height - H_SEED)
                    ),
                    "velocity_evolving_m_s": record["growth_velocity_m_s"],
                    "velocity_frozen_m_s": frozen_velocity,
                    "velocity_ratio": record["growth_velocity_m_s"] / frozen_velocity,
                    "minimum_porosity": record["minimum_porosity"],
                    "mean_porosity": record["mean_porosity"],
                }
            )
        index += 1
    last = rows[-1]
    summary = {
        "window_s": last["time_s"],
        "growth_evolving_m": last["growth_evolving_m"],
        "growth_frozen_m": last["growth_frozen_m"],
        "growth_ratio": last["growth_ratio_evolving_over_frozen"],
        "velocity_ratio_at_sealing": last["velocity_ratio"],
        "maximum_relative_height_difference": worst,
        "note": (
            "The frozen-porosity model keeps growing past the point at which "
            "the resolved crust has sealed, so beyond the sealing time the "
            "comparison is no longer between two admissible models."
        ),
    }
    return rows, summary


# --- summary ------------------------------------------------------------------


def build_summary(tables: dict, extra: dict) -> dict:
    reduction = tables["eq28_reduction"]
    linearised = [r for r in reduction if r["linearised_top_closure"]]
    exact = [r for r in reduction if not r["linearised_top_closure"]]
    orders = [
        r["observed_order"]
        for r in linearised
        if not math.isnan(r["observed_order"])
    ]
    law = tables["sealing_law"]
    thin = [r for r in law if r["thin_crust"] and r["sealed"] and r["resolved"]]
    ratios = [r["ratio_to_prediction"] for r in thin]
    unresolved = [r for r in law if not r["resolved"]]

    def _spread(key_fields: tuple, varied: str) -> float:
        groups: dict = {}
        for row in law:
            if not row["resolved"]:
                continue
            key = tuple(row[field] for field in key_fields)
            groups.setdefault(key, []).append(row["growth_m"])
        worst = 0.0
        for values in groups.values():
            if len(values) > 1 and min(values) > 0.0:
                worst = max(worst, (max(values) - min(values)) / min(values))
        return worst
    return {
        "generated": date.today().isoformat(),
        "python": platform.python_version(),
        "criteria_file": "verification/config/acceptance_criteria_p2b.json",
        "declared_reference": P2B,
        "reduction_to_equation_28": {
            "linearised_observed_orders": orders,
            "linearised_relative_at_512": linearised[-1]["relative_difference"],
            "exact_relative_at_512": exact[-1]["relative_difference"],
            "predicted_offset_c_star_over_c_eq": exact[-1][
                "predicted_offset_c_star_over_c_eq"
            ],
            "interpretation": (
                "With the same C ~ C_eq top linearisation the variable-porosity "
                "solver reproduces Prat equation (28) at second order. Without "
                "it the answer differs by C*_top/C_eq, which is the "
                "linearisation error carried by equation (28) itself."
            ),
        },
        "reference_trajectory": extra["reference_trajectory"],
        "sealing_law": {
            "rows": len(law),
            "all_sealed": all(r["sealed"] for r in law),
            "under_resolved_rows_excluded": len(unresolved),
            "resolution_criterion": (
                "at least 4 cells across the reaction boundary layer "
                "h/sqrt(Da) at the seed height"
            ),
            "thin_crust_resolved_rows": len(thin),
            "ratio_min": min(ratios) if ratios else None,
            "ratio_max": max(ratios) if ratios else None,
            "ratio_mean": sum(ratios) / len(ratios) if ratios else None,
            "growth_spread_across_relative_humidity": _spread(
                (
                    "specific_area_1_m",
                    "reaction_coefficient_m_s",
                    "deposition_porosity",
                ),
                "relative_humidity",
            ),
            "growth_spread_across_reaction_coefficient": _spread(
                (
                    "specific_area_1_m",
                    "relative_humidity",
                    "deposition_porosity",
                ),
                "reaction_coefficient_m_s",
            ),
            "growth_over_crystal_size_range": [
                min(r["growth_over_crystal_size"] for r in law),
                max(r["growth_over_crystal_size"] for r in law),
            ],
            "law": "dh = beta * eps_0 / a_v",
            "derivation": (
                "A material point seals after t = eps_0 rho_cr / "
                "(a_v rho k_r C*_top); the top advances at "
                "v_s = beta k_r rho C*_top / rho_cr; the product removes "
                "k_r, rho and C*_top, hence also the evaporation flux."
            ),
            "general_power_law": (
                "With a_v(eps) = a_v0 (eps/eps_0)**m the same cancellation "
                "holds for any m: dh = (beta eps_0**m / a_v0) * "
                "integral_{eps_seal}^{eps_0} eps**-m d eps."
            ),
            "cubic_corollary": (
                "With a_v = 6 (1-eps_0)/d for the ideal cubic construction "
                "that gives beta = 2 in the audit, dh/d = beta eps_0 / "
                "(6 (1-eps_0)) = 0.222 at eps_0 = 0.4."
            ),
        },
        "frozen_versus_evolving": extra["frozen_versus_evolving"],
        "closure_sensitivity": tables["closure_sensitivity"],
        "grid_convergence": {
            "sealing_times_s": {
                f"n{row['cells_per_seed_height']}_dt{int(row['dt_s'])}": row[
                    "sealing_time_s"
                ]
                for row in tables["grid_convergence"]
            },
            "clogging_depth_over_seed_height": [
                row["clogging_depth_over_seed_height"]
                for row in tables["grid_convergence"]
            ],
        },
        "limitations": CRITERIA_P2B["limitations"],
    }


def main(argv: list) -> int:
    out = (
        Path(argv[1])
        if len(argv) > 1
        else MODULE_ROOT / "results" / "p2b_evolving_porosity"
    )
    out.mkdir(parents=True, exist_ok=True)

    trajectory, snapshots, trajectory_summary = table_reference_trajectory()
    comparison, comparison_summary = table_frozen_versus_evolving()
    tables = {
        "eq28_reduction": table_eq28_reduction(),
        "reference_trajectory": trajectory,
        "porosity_snapshots": snapshots,
        "grid_convergence": table_grid_convergence(),
        "sealing_law": table_sealing_law(),
        "closure_sensitivity": table_closure_sensitivity(),
        "frozen_versus_evolving": comparison,
    }
    for name, rows in tables.items():
        write_csv(out / f"{name}.csv", rows)

    summary = build_summary(
        tables,
        {
            "reference_trajectory": trajectory_summary,
            "frozen_versus_evolving": comparison_summary,
        },
    )
    with (out / "summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2)
        handle.write("\n")

    print(f"wrote {len(tables)} tables and summary.json to {out}")
    for name, rows in tables.items():
        print(f"  {name}.csv: {len(rows)} rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
