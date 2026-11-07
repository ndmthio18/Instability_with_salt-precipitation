#!/usr/bin/env python3
"""Generate the P2-D coupling evidence.

Writes CSV tables and a JSON summary into `results/p2d_coupling/`. Every
parameter comes from `config/acceptance_criteria_p2d.json` or from the P1-R2
Case 36877 inputs. Nothing is fitted.

Usage::

    python run_p2d.py [output_directory]
"""

from __future__ import annotations

import csv
import json
import platform
import sys
import time
from datetime import date
from pathlib import Path

import reference_config_p2d as rc
from column import interface, soil

MODULE_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_ROWS = 200


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


def thin(history: list) -> list:
    if len(history) <= OUTPUT_ROWS:
        return list(history)
    stride = max(1, len(history) // OUTPUT_ROWS)
    kept = [history[i] for i in range(0, len(history), stride)]
    if kept[-1] is not history[-1]:
        kept.append(history[-1])
    return kept


# --- stage 1: the activity-controlled film --------------------------------


def stage_one(cells: int, dt: float, verbose: bool, horizon: float | None = None,
              maximum_steps: int = 4000) -> dict:
    state = rc.hydrostatic_state(cells=cells)
    boundary = soil.BoundaryCondition(
        top_flux=rc.film_flux(rc.film()), bottom="no_flux"
    )
    steps = maximum_steps if horizon is None else int(round(horizon / dt))
    start = time.time()
    report = soil.run(
        rc.SOIL, state, boundary, dt, steps, record_stride=5,
        stop=rc.solubility_stop,
    )
    report["audit"] = soil.audit(
        rc.SOIL, report["state"], report["initial"], report["transferred"]
    )
    report["wall_s"] = time.time() - start
    report["cells"] = cells
    report["dt_s"] = dt
    report["converged"] = bool(report["history"]) and "Error" not in report["stop_reason"]
    if verbose:
        print(
            f"   cells={cells:3d} dt={dt:6.0f} -> {report['stop_reason'][:46]} "
            f"at t={report['state'].time:.6g} s [{report['wall_s']:.0f} s]"
        )
    return report


def film_row(report: dict, label: str) -> dict:
    audit = report["audit"]
    ok = report["converged"]
    first = report["history"][0] if report["history"] else {}
    last = report["history"][-1] if report["history"] else {}
    row = {
        "study": label,
        "cells": report["cells"],
        "dt_s": report["dt_s"],
        "stop_reason": report["stop_reason"],
        "converged": ok,
        "reached_solubility": report["stop_reason"] == "soil_top_reached_solubility",
        "final_time_s": report["state"].time,
        "steps": len(report["history"]),
        "wall_s": report["wall_s"],
        "initial_flux_kg_m2_s": first.get("top_total_kg_m2_s"),
        "final_flux_kg_m2_s": last.get("top_total_kg_m2_s"),
        "final_top_mass_fraction": last.get("top_mass_fraction"),
        "final_top_over_c_eq": (
            last["top_mass_fraction"] / rc.C_EQ if last else None
        ),
        "final_top_saturation": last.get("top_saturation"),
        "evaporated_kg_m2": audit["transferred_water_kg_m2"] if ok else None,
        "salt_transferred_kg_m2": audit["transferred_salt_kg_m2"] if ok else None,
        "initial_salt_kg_m2": audit["initial_salt_kg_m2"],
        "final_salt_kg_m2": audit["final_salt_kg_m2"] if ok else None,
        "liquid_residual_relative": audit["liquid_residual_relative"] if ok else None,
        "salt_residual_relative": audit["salt_residual_relative"] if ok else None,
        "water_residual_relative": audit["water_residual_relative"] if ok else None,
        "worst_newton_residual_relative": (
            max(r["newton_residual_relative"] for r in report["history"])
            if report["history"] else None
        ),
        "note": "" if ok else (
            "the solve failed, so the inventories are not a converged state "
            "and no budget is reported"
        ),
    }
    if first and last:
        row["flux_decline_percent"] = (
            100.0 * (first["top_total_kg_m2_s"] - last["top_total_kg_m2_s"])
            / first["top_total_kg_m2_s"]
        )
    return row


def table_film(verbose: bool) -> tuple:
    """Grid and time-step studies at a common horizon, plus the step limit.

    Grid and time step are varied separately. Comparing 20/40/80 cells each at
    its own largest stable step would confound the two, so the grid study
    fixes dt at 250 s, which every resolution tested can take.
    """
    rows = []
    horizon = 50000.0

    print("   grid study (dt = 250 s, horizon 50000 s)")
    previous = None
    for cells in (20, 40, 80):
        report = stage_one(cells, 250.0, verbose, horizon=horizon)
        row = film_row(report, "grid")
        if previous is not None and row["final_flux_kg_m2_s"] and previous["final_flux_kg_m2_s"]:
            for key in ("final_flux_kg_m2_s", "final_top_mass_fraction"):
                row[f"successive_{key}"] = abs(row[key] - previous[key]) / abs(row[key])
        previous = row
        rows.append(row)

    print("   time-step study (40 cells, horizon 50000 s)")
    previous = None
    for dt in (1000.0, 500.0, 250.0, 125.0):
        report = stage_one(40, dt, verbose, horizon=horizon)
        row = film_row(report, "timestep")
        if previous is not None and row["final_flux_kg_m2_s"] and previous["final_flux_kg_m2_s"]:
            for key in ("final_flux_kg_m2_s", "final_top_mass_fraction"):
                row[f"successive_{key}"] = abs(row[key] - previous[key]) / abs(row[key])
        previous = row
        rows.append(row)

    print("   measured step limit")
    for cells, dt in ((20, 500.0), (80, 500.0)):
        rows.append(film_row(stage_one(cells, dt, verbose), "step_limit"))

    print("   production run to solubility (40 cells, dt = 500 s)")
    production = stage_one(40, 500.0, verbose)
    if not production["converged"]:
        raise RuntimeError("the production film run did not reach solubility")
    rows.append(film_row(production, "production"))
    return rows, production


# --- stages 2 and 3: the seed and the coupled interface --------------------


def table_coupled(production: dict, steps: int, dt: float, verbose: bool) -> tuple:
    column = production["state"].copy()
    pre_seed = column.copy()
    params = rc.crust_parameters()
    film = rc.film()
    salt_before = soil.inventory(rc.SOIL, column)["salt_kg_m2"]
    crust, seed_salt, seed_water = interface.seed_from_soil(
        rc.SOIL, column, params, film, rc.DECLARED["seed_height_m"]
    )
    salt_after = soil.inventory(rc.SOIL, column)["salt_kg_m2"]
    seed = {
        "seed_height_m": crust.height,
        "seed_cells": len(crust.thicknesses),
        "seed_top_cell_m": crust.top_cell_size,
        "reaction_length_m": params.reaction_length(),
        "seed_over_reaction_length": crust.height / params.reaction_length(),
        "seed_salt_kg_m2": seed_salt,
        "seed_water_kg_m2": seed_water,
        "soil_salt_before_kg_m2": salt_before,
        "soil_salt_after_kg_m2": salt_after,
        "seed_fraction_of_column_salt": seed_salt / salt_before,
        "seed_neutrality_residual_relative": abs(
            (salt_after + seed_salt) - salt_before
        )
        / salt_before,
        "crust_minimum_porosity": crust.minimum_porosity,
        "connectivity_threshold": params.connectivity_porosity,
    }
    if verbose:
        print(
            f"   seed: h0={crust.height:.3g} m, {seed['seed_cells']} cells, "
            f"{100 * seed['seed_fraction_of_column_salt']:.2f} % of the column "
            f"salt, neutrality {seed['seed_neutrality_residual_relative']:.2e}"
        )

    state = interface.CoupledState(
        column=column, crust=crust,
        seed_salt_kg_m2=seed_salt, seed_water_kg_m2=seed_water,
    )
    start = time.time()
    reason = "steps"
    for index in range(steps):
        try:
            record = interface.step(
                rc.SOIL, state, params, film, dt,
                undersaturation_limit=rc.DECLARED["interface_undersaturation_limit"],
            )
        except interface.InterfaceDomainExit as error:
            reason = f"InterfaceDomainExit: {error}"
            break
        except (soil.ColumnSolverError, soil.ColumnDomainExit) as error:
            reason = f"{type(error).__name__}: {error}"
            break
        if verbose and (index % 10 == 0 or index == steps - 1):
            print(
                f"      t={record['time_s']:7.1f} w_i/Ceq={record['interface_over_c_eq']:.9f} "
                f"v_s={record['growth_velocity_m_s']:.4e} h={record['crust_height_m']:.8g} "
                f"eps_min={record['crust_minimum_porosity']:.5f} "
                f"r23={record['route_2_minus_3_relative']:.2e}"
            )
    budget = interface.audit(rc.SOIL, state, pre_seed, seed_salt, seed_water)
    summary = {
        "steps_requested": steps,
        "steps_taken": len(state.history),
        "dt_s": dt,
        "stop_reason": reason,
        "wall_s": time.time() - start,
        **seed,
        **budget,
    }
    if state.history:
        last = state.history[-1]
        summary.update(
            {
                "final_time_s": last["time_s"],
                "final_interface_over_c_eq": last["interface_over_c_eq"],
                "final_growth_velocity_m_s": last["growth_velocity_m_s"],
                "final_crust_height_m": last["crust_height_m"],
                "crust_growth_m": last["crust_height_m"] - seed["seed_height_m"],
                "final_minimum_porosity": last["crust_minimum_porosity"],
                "final_evaporation_kg_m2_s": last["evaporation_flux_kg_m2_s"],
                "worst_interface_residual_relative": max(
                    row["interface_residual_relative"] for row in state.history
                ),
                "worst_route_2_minus_3_relative": max(
                    row["route_2_minus_3_relative"] for row in state.history
                ),
                "worst_route_1_minus_2_relative": max(
                    row["route_1_minus_2_relative"] for row in state.history
                ),
                "mean_route_1_minus_2_predicted_relative": sum(
                    row["route_1_minus_2_predicted_relative"] for row in state.history
                )
                / len(state.history),
                "worst_peclet_over_damkohler": max(
                    row["peclet_over_damkohler"] for row in state.history
                ),
            }
        )
    return summary, state


def table_timestep(production: dict, verbose: bool) -> list:
    rows = []
    previous = None
    for dt in (20.0, 10.0, 5.0, 2.5):
        steps = int(100.0 / dt)
        summary, _ = table_coupled(production, steps, dt, False)
        row = {
            "dt_s": dt,
            "steps": summary["steps_taken"],
            "final_time_s": summary.get("final_time_s"),
            "final_crust_height_m": summary.get("final_crust_height_m"),
            "final_growth_velocity_m_s": summary.get("final_growth_velocity_m_s"),
            "final_minimum_porosity": summary.get("final_minimum_porosity"),
            "salt_residual_relative": summary["salt_residual_relative"],
            "water_residual_relative": summary["water_residual_relative"],
            "worst_route_2_minus_3_relative": summary.get(
                "worst_route_2_minus_3_relative"
            ),
            "stop_reason": summary["stop_reason"],
        }
        if previous is not None:
            for key in ("final_crust_height_m", "final_growth_velocity_m_s",
                        "final_minimum_porosity"):
                if row[key] and previous[key]:
                    row[f"successive_{key}"] = abs(row[key] - previous[key]) / abs(
                        row[key]
                    )
        previous = row
        rows.append(row)
        if verbose:
            print(
                f"   dt={dt:5.1f} -> salt residual {row['salt_residual_relative']:.3e}, "
                f"h={row['final_crust_height_m']:.9g}"
            )
    return rows


def main(argv: list) -> int:
    output = (
        Path(argv[1]) if len(argv) > 1
        else MODULE_ROOT / "results" / "p2d_coupling"
    )
    output.mkdir(parents=True, exist_ok=True)
    verbose = True
    started = time.time()

    print("1/4 activity-controlled film evaporation (items 4 and 5)")
    film_rows, production = table_film(verbose)
    write_csv(output / "film_evaporation.csv", film_rows)
    write_csv(output / "film_trajectory.csv", thin(production["history"]))

    print("2/4 seed and coupled interface (items 6, 7, 8, 9)")
    summary, state = table_coupled(production, 60, 5.0, verbose)
    write_csv(output / "coupled_summary.csv", [summary])
    if state.history:
        write_csv(output / "coupled_trajectory.csv", thin(state.history))

    print("3/4 coupled time-step refinement")
    timestep_rows = table_timestep(production, verbose)
    write_csv(output / "coupled_timestep.csv", timestep_rows)

    print("4/4 environment and provenance")
    import shutil

    environment = {
        "cpp_tools": {
            name: shutil.which(name)
            for name in ("g++", "clang++", "cl", "cmake", "make", "ninja")
        },
        "python": platform.python_version(),
        "platform": platform.platform(),
    }

    headline = {
        "solubility_reached_s": production["state"].time,
        "film_production_flux_decline_percent": film_rows[-1]["flux_decline_percent"],
        "film_production_salt_residual_relative": film_rows[-1][
            "salt_residual_relative"
        ],
        "film_production_liquid_residual_relative": film_rows[-1][
            "liquid_residual_relative"
        ],
        "film_step_limit_note": (
            "dt = 500 s converges on 40 cells and fails on 80; dt = 250 s "
            "converges on 80. The limit scales with the surface cell "
            "thickness."
        ),
        "seed_fraction_of_column_salt": summary["seed_fraction_of_column_salt"],
        "seed_neutrality_residual_relative": summary[
            "seed_neutrality_residual_relative"
        ],
        "coupled_steps": summary["steps_taken"],
        "coupled_stop_reason": summary["stop_reason"],
        "interface_over_c_eq": summary.get("final_interface_over_c_eq"),
        "coupled_salt_residual_relative": summary["salt_residual_relative"],
        "coupled_water_residual_relative": summary["water_residual_relative"],
        "worst_route_2_minus_3_relative": summary.get(
            "worst_route_2_minus_3_relative"
        ),
        "worst_route_1_minus_2_relative": summary.get(
            "worst_route_1_minus_2_relative"
        ),
        "predicted_route_1_minus_2_relative": summary.get(
            "mean_route_1_minus_2_predicted_relative"
        ),
        "cpp_compiled": False,
    }
    report = {
        "generated": date.today().isoformat(),
        "wall_seconds": time.time() - started,
        "criteria": "config/acceptance_criteria_p2d.json",
        "environment": environment,
        "headline": headline,
        "tables": sorted(p.name for p in output.glob("*.csv")),
    }
    with (output / "summary.json").open("w", encoding="utf-8") as handle:
        json.dump(report, handle, indent=2, sort_keys=True)
    print(json.dumps(headline, indent=2, sort_keys=True))
    print(f"\nwrote {len(report['tables'])} tables to {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
