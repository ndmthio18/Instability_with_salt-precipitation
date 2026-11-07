#!/usr/bin/env python3
"""Generate the P2-C salt_grow parameter-transfer evidence.

Writes CSV tables and a JSON summary into
``results/p2c_salt_grow_parameters/``. Every parameter comes from
``config/acceptance_criteria_p2c.json`` (registered before any run) or from
the P2-A declared reference. Nothing is fitted.

Usage::

    python run_p2c.py [output_directory]
"""

from __future__ import annotations

import csv
import json
import math
import platform
import sys
import time
from datetime import date
from pathlib import Path

import reference_config_p2c as rc
from reference_config import PARAMS
from surface import graded_crust as gc

MODULE_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_ROWS = 200
REF_A = rc.REFERENCE["published_velocities_m_s"]["prat_analytical"]
REF_C = rc.REFERENCE["published_velocities_m_s"]["m3_b2_corrected_6h"]
REF_B2 = rc.REFERENCE["published_b2_6h"]


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
    """Committed trajectories are thinned; the integration never is.

    The stride is chosen so every trajectory keeps about ``OUTPUT_ROWS``
    samples whatever its length, and the final record - which carries the stop
    event - is always kept.
    """
    if len(history) <= OUTPUT_ROWS:
        return list(history)
    stride = max(1, len(history) // OUTPUT_ROWS)
    kept = [history[index] for index in range(0, len(history), stride)]
    if kept[-1] is not history[-1]:
        kept.append(history[-1])
    return kept


def summarise(run: dict, params: gc.GradedCrustParameters, height0: float) -> dict:
    state = run["state"]
    record = state.history[-1]
    report = gc.budget(state, run["initial"])
    return {
        "stop_reason": run["stop_reason"],
        "steps": run["steps"],
        "time_s": record["time_s"],
        "dt_s": run["dt_s"],
        "tau_porosity_s": run["tau_porosity_s"],
        "tau_burial_s": run["tau_burial_s"],
        "top_cell_m": run["top_cell_m"],
        "cells_per_delta": run["cells_per_delta_actual"],
        "cells_final": record["cells"],
        "merges": state.merges,
        "crust_height_m": state.height,
        "growth_m": state.height - height0,
        "growth_relative_to_seed": (state.height - height0) / height0,
        "growth_velocity_m_s": record["growth_velocity_m_s"],
        "growth_velocity_mm_day": record["growth_velocity_m_s"] * 86400.0 * 1e3,
        "growth_velocity_ratio_to_initial": (
            record["growth_velocity_m_s"] / state.history[0]["growth_velocity_m_s"]
        ),
        "top_supersaturation": record["top_supersaturation"],
        "top_mass_fraction": record["top_mass_fraction"],
        "minimum_porosity": record["minimum_porosity"],
        "top_porosity": record["top_porosity"],
        "mean_porosity": record["mean_porosity"],
        "reaction_length_m": record["reaction_length_m"],
        "cells_across_reaction_length": record["cells_across_reaction_length"],
        "reaction_depth_95_m": record["reaction_depth_95_m"],
        "internal_reaction_kg_m2_s": record["internal_reaction_kg_m2_s"],
        "internal_reaction_clamped": record["internal_reaction_clamped"],
        "top_skeleton_kg_m2_s": record["top_skeleton_kg_m2_s"],
        "internal_over_top": record["internal_over_top"],
        "internal_fraction_of_new_solid": report["internal_fraction_of_new_solid"],
        "salt_flux_soil_kg_m2_s": record["salt_flux_soil_kg_m2_s"],
        "water_flux_soil_kg_m2_s": record["water_flux_soil_kg_m2_s"],
        "solid_salt_kg_m2": record["solid_salt_kg_m2"],
        "dissolved_salt_kg_m2": record["dissolved_salt_kg_m2"],
        "total_salt_kg_m2": record["solid_salt_kg_m2"] + record["dissolved_salt_kg_m2"],
        "water_kg_m2": record["water_kg_m2"],
        "route_difference_relative": record["route_difference_relative"],
        "route_difference_legacy_relative": record["route_difference_legacy_relative"],
        "route_difference_identity_kg_m2_s": record["route_difference_identity_kg_m2_s"],
        "route_identity_error_relative": record["route_identity_error_relative"],
        "salt_residual_relative": report["salt_residual_relative"],
        "water_residual_relative": report["water_residual_relative"],
        "solid_split_residual_relative": report["solid_split_residual_relative"],
        "disconnected": state.disconnected,
        "disconnected_time_s": state.disconnected_time,
        "sealed": state.sealed,
        "sealed_time_s": state.sealed_time,
    }


def profile_rows(state: gc.GradedCrustState, params, flux, label: str) -> list:
    solution = gc.solve_profile(state, params, flux)
    rows = []
    for index in range(len(solution.centres)):
        rows.append(
            {
                "label": label,
                "z_m": solution.centres[index],
                "depth_below_front_m": state.height - solution.centres[index],
                "thickness_m": solution.thicknesses[index],
                "porosity": solution.porosity[index],
                "supersaturation": solution.profile[index],
                "mass_fraction": params.equilibrium_mass_fraction
                + solution.profile[index],
                "reaction_kg_m3_s": params.area(solution.porosity[index])
                * params.rho
                * params.reaction_coefficient
                * solution.profile[index],
            }
        )
    # Thin deep cells only; the reaction layer is always kept in full.
    keep = [r for r in rows if r["depth_below_front_m"] <= 20.0 * solution.reaction_length]
    deep = [r for r in rows if r["depth_below_front_m"] > 20.0 * solution.reaction_length]
    return keep + deep[:: max(1, len(deep) // 60)]


# --- table 1: parameters -------------------------------------------------------


def table_parameters() -> list:
    rows = [
        {
            "quantity": "deposition / crust porosity eps_0",
            "symbol": "eps_0",
            "p2b_value": PARAMS["crust_porosity"],
            "salt_grow_value": rc.SALT_GROW["epsilon_ref"],
            "unit": "-",
            "transferred": True,
            "p2b_provenance": "declared in acceptance_criteria.json",
            "salt_grow_provenance": "configs/baseline.json physical.epsilon_ref",
        },
        {
            "quantity": "reaction coefficient",
            "symbol": "k_r",
            "p2b_value": PARAMS["reaction_coefficient"],
            "salt_grow_value": rc.SALT_GROW["k_reaction"],
            "unit": "m/s",
            "transferred": True,
            "p2b_provenance": "declared placeholder; the audit supplies no value",
            "salt_grow_provenance": "configs/baseline.json physical.k_reaction",
        },
        {
            "quantity": "specific crystal area",
            "symbol": "a_v",
            "p2b_value": PARAMS["specific_area"],
            "salt_grow_value": rc.SALT_GROW["a_v_ref"],
            "unit": "1/m",
            "transferred": True,
            "p2b_provenance": "declared placeholder; the audit supplies no value",
            "salt_grow_provenance": "36*(1-eps_0)/d_sub with d_sub = 25 um",
        },
        {
            "quantity": "substrate / crystal size",
            "symbol": "d_sub",
            "p2b_value": None,
            "salt_grow_value": rc.SALT_GROW["d_sub"],
            "unit": "m",
            "transferred": True,
            "p2b_provenance": "not used; P2-B's cubic corollary used a_v = 6(1-eps_0)/d",
            "salt_grow_provenance": "configs/baseline.json physical.d_sub",
        },
        {
            "quantity": "equilibrium mass fraction",
            "symbol": "C_eq",
            "p2b_value": rc.C_EQ,
            "salt_grow_value": rc.SALT_GROW["C_eq"],
            "unit": "-",
            "transferred": False,
            "p2b_provenance": "frozen DuMuX SolubilityMoleFraction = 0.1",
            "salt_grow_provenance": "configs/baseline.json physical.C_eq",
        },
        {
            "quantity": "effective salt diffusivity",
            "symbol": "D_s,eff",
            "p2b_value": PARAMS["effective_salt_diffusivity"],
            "salt_grow_value": rc.SALT_GROW["D_eff_ref"],
            "unit": "m2/s",
            "transferred": False,
            "p2b_provenance": "declared in acceptance_criteria.json",
            "salt_grow_provenance": "sqrt(eps)*D_mol with D_mol = 1.3e-9",
        },
        {
            "quantity": "solution density",
            "symbol": "rho",
            "p2b_value": rc.case_parameters(0).rho,
            "salt_grow_value": rc.SALT_GROW["rho_solution"],
            "unit": "kg/m3",
            "transferred": False,
            "p2b_provenance": "BrineDensityLaw at C_eq",
            "salt_grow_provenance": "configs/baseline.json physical.rho_solution",
        },
        {
            "quantity": "crystal density",
            "symbol": "rho_cr",
            "p2b_value": PARAMS["crystal_density"],
            "salt_grow_value": rc.SALT_GROW["rho_crystal"],
            "unit": "kg/m3",
            "transferred": False,
            "p2b_provenance": "DuMuX nacl.hh",
            "salt_grow_provenance": "configs/baseline.json physical.rho_crystal",
        },
        {
            "quantity": "upscaling factor",
            "symbol": "beta",
            "p2b_value": PARAMS["beta"],
            "salt_grow_value": rc.SALT_GROW["beta"],
            "unit": "-",
            "transferred": False,
            "p2b_provenance": "Prat appendix cubic construction",
            "salt_grow_provenance": "configs/baseline.json physical.beta",
        },
        {
            "quantity": "initial crust height",
            "symbol": "h0",
            "p2b_value": rc.H_SEED,
            "salt_grow_value": rc.SALT_GROW["h0"],
            "unit": "m",
            "transferred": False,
            "p2b_provenance": "declared seed height",
            "salt_grow_provenance": "configs/baseline.json physical.h0",
        },
        {
            "quantity": "evaporation mass flux",
            "symbol": "j",
            "p2b_value": rc.FILM_FLUX,
            "salt_grow_value": rc.SALT_GROW["evaporation_mass_flux"],
            "unit": "kg/m2/s",
            "transferred": False,
            "p2b_provenance": "P2-A gas film at C_eq, physical inputs only",
            "salt_grow_provenance": "rho * U_evap with U_evap = 1 cm/day",
        },
    ]
    for row in rows:
        if row["p2b_value"] and row["salt_grow_value"]:
            row["ratio_salt_grow_over_p2b"] = (
                row["salt_grow_value"] / row["p2b_value"]
            )
    return rows


# --- table 2: dimensionless groups ---------------------------------------------


def table_groups() -> list:
    rows = []
    for case in sorted(rc.CASES):
        params = rc.case_parameters(case)
        groups = rc.dimensionless_groups(params, rc.H_SEED, rc.FILM_FLUX)
        rows.append(
            {
                "configuration": f"Case {case}: {rc.CASES[case]['label']}",
                "height_m": rc.H_SEED,
                "evaporation_flux_kg_m2_s": rc.FILM_FLUX,
                **groups,
            }
        )
    for closure in ("B1", "B2"):
        params = rc.matched_parameters(closure)
        groups = rc.dimensionless_groups(
            params, rc.SALT_GROW["h0"], rc.SALT_GROW["evaporation_mass_flux"]
        )
        rows.append(
            {
                "configuration": f"matched salt_grow {closure}",
                "height_m": rc.SALT_GROW["h0"],
                "evaporation_flux_kg_m2_s": rc.SALT_GROW["evaporation_mass_flux"],
                **groups,
            }
        )
    base = rows[0]
    for row in rows:
        for key in ("Da", "Da_star", "delta_r_m", "growth_velocity_m_s"):
            row[f"{key}_over_case0"] = row[key] / base[key]
    return rows


# --- table 3: grid convergence -------------------------------------------------


def table_grid_convergence(verbose: bool) -> list:
    rows = []
    settings = [
        ("Case 3", rc.case_parameters(3), rc.FILM_FLUX, rc.H_SEED, 0.0),
        (
            "matched B2",
            rc.matched_parameters("B2"),
            rc.SALT_GROW["evaporation_mass_flux"],
            rc.SALT_GROW["h0"],
            40.0,
        ),
    ]
    for label, params, flux, height, freeze in settings:
        previous = None
        for per in (1.0, 2.0, 4.0, 8.0, 16.0):
            start = time.time()
            run = rc.run_case(
                params,
                flux,
                height,
                cells_per_delta=per,
                freeze_lengths=freeze,
            )
            row = {
                "configuration": label,
                "cells_per_delta_requested": per,
                "wall_s": time.time() - start,
                **summarise(run, params, height),
            }
            row["resolved"] = row["cells_per_delta"] >= rc.TOL_P2C[
                "boundary_layer_cells_minimum"
            ]
            if previous is not None:
                for key in (
                    "growth_velocity_m_s",
                    "top_supersaturation",
                    "internal_reaction_kg_m2_s",
                    "minimum_porosity",
                    "reaction_depth_95_m",
                ):
                    row[f"successive_{key}"] = abs(
                        row[key] - previous[key]
                    ) / max(abs(row[key]), 1.0e-300)
            previous = row
            rows.append(row)
            if verbose:
                print(
                    f"   grid {label} per={per:4.1f} -> v={row['growth_velocity_m_s']:.8e} "
                    f"eps_min={row['minimum_porosity']:.6f} "
                    f"({row['wall_s']:.1f} s)"
                )
    # The committed P2-B grid, shown explicitly as under-resolved.
    params = rc.case_parameters(3)
    run = rc.run_case(params, rc.FILM_FLUX, rc.H_SEED, cells_per_delta=1.0e-9)
    row = {
        "configuration": "Case 3 on the committed 32-cell P2-B seed grid",
        "cells_per_delta_requested": None,
        "wall_s": 0.0,
        **summarise(run, params, rc.H_SEED),
    }
    row["resolved"] = False
    rows.append(row)
    return rows


# --- table 4: time-step convergence -------------------------------------------


def table_timestep_convergence(verbose: bool) -> list:
    rows = []
    settings = [
        ("Case 0", rc.case_parameters(0), rc.FILM_FLUX, rc.H_SEED, 0.0),
        ("Case 3", rc.case_parameters(3), rc.FILM_FLUX, rc.H_SEED, 0.0),
        (
            "matched B2",
            rc.matched_parameters("B2"),
            rc.SALT_GROW["evaporation_mass_flux"],
            rc.SALT_GROW["h0"],
            40.0,
        ),
    ]
    for label, params, flux, height, freeze in settings:
        previous = None
        for steps in (375.0, 750.0, 1500.0, 3000.0, 6000.0):
            run = rc.run_case(
                params,
                flux,
                height,
                steps_per_tau=steps,
                maximum_steps=200000,
                freeze_lengths=freeze,
            )
            row = {
                "configuration": label,
                "steps_per_tau": steps,
                **summarise(run, params, height),
            }
            if previous is not None:
                for key in (
                    "growth_velocity_m_s",
                    "minimum_porosity",
                    "route_difference_relative",
                    "route_difference_legacy_relative",
                ):
                    row[f"successive_{key}"] = abs(
                        row[key] - previous[key]
                    ) / max(abs(row[key]), 1.0e-300)
            previous = row
            rows.append(row)
            if verbose:
                print(
                    f"   dt  {label} steps/tau={steps:6.0f} -> dt={row['dt_s']:.5g} s "
                    f"route_corr={row['route_difference_relative']:.4e} "
                    f"legacy={row['route_difference_legacy_relative']:.4e}"
                )
    return rows


# --- table 5: the one-factor case ladder ---------------------------------------


def table_cases(output: Path, verbose: bool) -> tuple:
    rows = []
    profiles = []
    for case in sorted(rc.CASES):
        params = rc.case_parameters(case)
        groups = rc.dimensionless_groups(params, rc.H_SEED, rc.FILM_FLUX)
        run = rc.run_case(params, rc.FILM_FLUX, rc.H_SEED)
        row = {
            "case": case,
            "label": rc.CASES[case]["label"],
            "eps_0": params.deposition_porosity,
            "k_r": params.reaction_coefficient,
            "a_v": params.specific_area,
            "Da": groups["Da"],
            "Da_star": groups["Da_star"],
            "sqrt_Da_star": groups["sqrt_Da_star"],
            "delta_r_m": groups["delta_r_m"],
            "analytic_velocity_m_s": groups["growth_velocity_m_s"],
            "predicted_seals": groups["seals_predicted"],
            "predicted_wake_porosity": groups["predicted_wake_porosity"],
            "predicted_sealing_growth_m": groups["predicted_sealing_growth_m"],
            **summarise(run, params, rc.H_SEED),
        }
        row["wake_prediction_error_relative"] = (
            abs(row["minimum_porosity"] - row["predicted_wake_porosity"])
            / max(abs(row["minimum_porosity"]), 1.0e-300)
        )
        row["sealing_growth_ratio"] = (
            row["growth_m"] / row["predicted_sealing_growth_m"]
            if row["predicted_seals"]
            else None
        )
        row["sealing_prediction_correct"] = bool(
            row["predicted_seals"] == bool(row["disconnected"])
        )
        rows.append(row)
        write_csv(output / f"case_{case}_trajectory.csv", thin(run["state"].history))
        profiles.extend(
            profile_rows(run["state"], params, rc.FILM_FLUX, f"Case {case}")
        )
        if verbose:
            print(
                f"   case {case}: stop={row['stop_reason']} "
                f"v={row['growth_velocity_m_s']:.6e} eps_min={row['minimum_porosity']:.6f} "
                f"predicted seal={row['predicted_seals']} observed={bool(row['disconnected'])}"
            )
    return rows, profiles


# --- table 6: the fully matched benchmark --------------------------------------


def table_matched(output: Path, verbose: bool) -> tuple:
    rows = []
    profiles = []
    height = rc.SALT_GROW["h0"]
    flux = rc.SALT_GROW["evaporation_mass_flux"]
    analytic = rc.dimensionless_groups(rc.matched_parameters("B2"), height, flux)
    rows.append(
        {
            "model": "A  Prat analytical growth law (this repository)",
            "closure": "uniform frozen eps_0",
            "growth_velocity_m_s": analytic["growth_velocity_m_s"],
            "growth_velocity_mm_day": analytic["growth_velocity_mm_day"],
            "top_supersaturation": analytic["c_star_top"],
            "minimum_porosity": rc.SALT_GROW["epsilon_ref"],
            "internal_fraction_of_new_solid": analytic["internal_over_top"]
            / (1.0 + analytic["internal_over_top"]),
            "source": "closed form, equations (28)/(29)",
        }
    )
    rows.append(
        {
            "model": "A' Prat analytical growth law (salt_grow published)",
            "closure": "uniform frozen eps_0",
            "growth_velocity_m_s": REF_A,
            "growth_velocity_mm_day": REF_A * 86400e3,
            "source": "salt_grow RESULTS.md",
        }
    )
    for closure in ("B1", "B2"):
        params = rc.matched_parameters(closure)
        run = rc.run_case(params, flux, height, cells_per_delta=16.0, freeze_lengths=40.0)
        row = {
            "model": f"B  P2-B quasi-steady evolving crust ({closure})",
            "closure": (
                "a_v, D_eff frozen at eps_0"
                if closure == "B1"
                else "a_v = 36(1-eps)/d_sub, D_eff = sqrt(eps) D_mol"
            ),
            "source": "this repository, surface/graded_crust.py",
            **summarise(run, params, height),
        }
        rows.append(row)
        profiles.extend(profile_rows(run["state"], params, flux, f"matched {closure}"))
        write_csv(output / f"matched_{closure}_trajectory.csv", thin(run["state"].history))
        if verbose:
            print(
                f"   matched {closure}: v={row['growth_velocity_m_s']:.8e} "
                f"({100*(row['growth_velocity_m_s']-REF_C)/REF_C:+.4f}% vs salt_grow M3)"
            )
    rows.append(
        {
            "model": "C  salt_grow M3 transient ALE (B2, 640 cells, 6 h)",
            "closure": "a_v = 36(1-eps)/d_sub, D_eff = sqrt(eps) D_mol",
            "growth_velocity_m_s": REF_C,
            "growth_velocity_mm_day": REF_C * 86400e3,
            "top_supersaturation": REF_B2["top_excess_mass_fraction"],
            "minimum_porosity": REF_B2["minimum_porosity"],
            "crust_height_m": REF_B2["height_m"],
            "reaction_depth_95_m": REF_B2["L95_m"],
            "internal_fraction_of_new_solid": REF_B2["internal_fraction_of_new_solid"],
            "source": "salt_grow RESULTS.md, M3 corrected baseline",
        }
    )
    for row in rows:
        velocity = row.get("growth_velocity_m_s")
        if velocity:
            row["difference_from_A_percent"] = 100.0 * (velocity - REF_A) / REF_A
            row["difference_from_C_percent"] = 100.0 * (velocity - REF_C) / REF_C
    return rows, profiles


# --- table 7: conservation and connectivity ------------------------------------


def table_conservation(verbose: bool) -> list:
    rows = []
    settings = [
        ("Case 0", rc.case_parameters(0), rc.FILM_FLUX, rc.H_SEED, 0.0),
        ("Case 2", rc.case_parameters(2), rc.FILM_FLUX, rc.H_SEED, 0.0),
        ("Case 3", rc.case_parameters(3), rc.FILM_FLUX, rc.H_SEED, 0.0),
        (
            "matched B2",
            rc.matched_parameters("B2"),
            rc.SALT_GROW["evaporation_mass_flux"],
            rc.SALT_GROW["h0"],
            40.0,
        ),
    ]
    for label, params, flux, height, freeze in settings:
        run = rc.run_case(params, flux, height, freeze_lengths=freeze)
        state = run["state"]
        report = gc.budget(state, run["initial"])
        corrected = [r["route_difference_relative"] for r in state.history]
        legacy = [r["route_difference_legacy_relative"] for r in state.history]
        last = state.history[-1]
        predicted_legacy = (
            params.equilibrium_mass_fraction
            * params.rho
            * params.crust_saturation
            / params.crystal_density
            * last["internal_reaction_kg_m2_s"]
            / max(abs(last["salt_flux_soil_kg_m2_s"]), 1.0e-300)
        )
        rows.append(
            {
                "configuration": label,
                "steps": run["steps"],
                "dt_s": run["dt_s"],
                "stop_reason": run["stop_reason"],
                "route_difference_corrected_final": corrected[-1],
                "route_difference_corrected_max": max(corrected),
                "route_difference_corrected_median": sorted(corrected)[
                    len(corrected) // 2
                ],
                "route_difference_legacy_final": legacy[-1],
                "route_difference_legacy_max": max(legacy),
                "route_difference_legacy_predicted": predicted_legacy,
                "legacy_prediction_ratio": legacy[-1] / max(predicted_legacy, 1.0e-300),
                "route_difference_kg_m2_s": last["route_difference_kg_m2_s"],
                "route_difference_identity_kg_m2_s": last[
                    "route_difference_identity_kg_m2_s"
                ],
                "route_identity_error_relative": last["route_identity_error_relative"],
                "solution_mass_flux_top_relative_kg_m2_s": last[
                    "solution_mass_flux_top_relative_kg_m2_s"
                ],
                "salt_flux_base_transport_kg_m2_s": last[
                    "salt_flux_base_route_kg_m2_s"
                ],
                "salt_flux_base_inventory_kg_m2_s": last["salt_flux_soil_kg_m2_s"],
                "water_flux_base_kg_m2_s": last["water_flux_soil_kg_m2_s"],
                "solution_mass_flux_base_kg_m2_s": last[
                    "solution_mass_flux_base_kg_m2_s"
                ],
                "salt_residual_relative": report["salt_residual_relative"],
                "water_residual_relative": report["water_residual_relative"],
                "solid_split_residual_relative": report["solid_split_residual_relative"],
                "internal_reaction_clamped_ever": any(
                    r["internal_reaction_clamped"] for r in state.history
                ),
                "connectivity_loss_time_s": state.disconnected_time,
                "clogging_time_s": state.sealed_time,
                "connectivity_before_clogging": bool(
                    state.disconnected and not state.sealed
                ),
                "minimum_porosity": last["minimum_porosity"],
                "connectivity_threshold": params.connectivity_porosity,
                "clogging_threshold": params.clogging_porosity,
            }
        )
        if verbose:
            print(
                f"   conservation {label}: corrected={corrected[-1]:.3e} "
                f"legacy={legacy[-1]:.3e} predicted={predicted_legacy:.3e}"
            )
    return rows


# --- table 8: the long matched run ---------------------------------------------


def table_long_run(output: Path, verbose: bool) -> list:
    """One matched B2 run carried to the salt_grow reporting time of 6 h.

    Uses a declared larger time step and wake coarsening, both justified by
    the convergence tables, so the comparison is made at the same time as the
    published reference rather than at a plateau declared earlier.
    """
    params = rc.matched_parameters("B2")
    height = rc.SALT_GROW["h0"]
    flux = rc.SALT_GROW["evaporation_mass_flux"]
    rows = []
    for dt in (0.25, 0.5):
        start = time.time()
        run = rc.run_case(
            params,
            flux,
            height,
            cells_per_delta=8.0,
            time_step=dt,
            maximum_steps=200000,
            horizon=21600.0,
            plateau_windows=1.0e9,
            freeze_lengths=40.0,
        )
        state = run["state"]
        row = {
            "configuration": "matched B2 to 6 h",
            "time_step_s": dt,
            "wall_s": time.time() - start,
            **summarise(run, params, height),
        }
        row["velocity_difference_from_C_percent"] = (
            100.0 * (row["growth_velocity_m_s"] - REF_C) / REF_C
        )
        row["height_difference_from_C_percent"] = (
            100.0 * (row["crust_height_m"] - REF_B2["height_m"]) / REF_B2["height_m"]
        )
        row["minimum_porosity_difference_from_C_percent"] = (
            100.0
            * (row["minimum_porosity"] - REF_B2["minimum_porosity"])
            / REF_B2["minimum_porosity"]
        )
        rows.append(row)
        if dt == 0.5:
            write_csv(output / "matched_B2_6h_trajectory.csv", thin(state.history))
        if verbose:
            print(
                f"   6 h run dt={dt}: t={row['time_s']:.6g} s v={row['growth_velocity_m_s']:.8e} "
                f"({row['velocity_difference_from_C_percent']:+.4f}% vs C) "
                f"h={row['crust_height_m']:.8g} ({row['wall_s']:.0f} s)"
            )
    return rows


# --- table 9: porosity and supersaturation profiles in time --------------------


def table_profile_snapshots(verbose: bool) -> list:
    """eps(z, t) and C*(z, t) at declared fractions of each case's run.

    The declared times are fractions of the time each case reached in the
    production run of table 5, so every case is sampled over its own history
    rather than at an arbitrary wall-clock time.
    """
    rows = []
    settings = [
        (f"Case {case}", rc.case_parameters(case), rc.FILM_FLUX, rc.H_SEED, 0.0)
        for case in sorted(rc.CASES)
    ]
    settings.append(
        (
            "matched B2",
            rc.matched_parameters("B2"),
            rc.SALT_GROW["evaporation_mass_flux"],
            rc.SALT_GROW["h0"],
            40.0,
        )
    )
    for label, params, flux, height, freeze in settings:
        probe = rc.run_case(params, flux, height, cells_per_delta=4.0,
                            freeze_lengths=freeze)
        total = probe["state"].history[-1]["time_s"]
        times = [fraction * total for fraction in (0.05, 0.25, 0.5, 1.0)]
        captured = []

        def capture(when, state, _label=label, _params=params, _flux=flux,
                    _captured=captured):
            _captured.append(
                (when, profile_rows(state, _params, _flux, f"{_label} t={when:.6g}s"))
            )

        rc.run_case(params, flux, height, cells_per_delta=4.0,
                    freeze_lengths=freeze, snapshot_times=times,
                    on_snapshot=capture)
        for when, block in captured:
            for row in block:
                row["time_s"] = when
                row["configuration"] = label
            rows.extend(block)
        if verbose:
            print(f"   profiles {label}: {len(captured)} snapshots to t={total:.6g} s")
    return rows


# --- main ----------------------------------------------------------------------


def main(argv: list) -> int:
    output = (
        Path(argv[1])
        if len(argv) > 1
        else MODULE_ROOT / "results" / "p2c_salt_grow_parameters"
    )
    output.mkdir(parents=True, exist_ok=True)
    verbose = True
    started = time.time()

    print("1/9 parameter comparison")
    parameters = table_parameters()
    write_csv(output / "parameter_comparison.csv", parameters)

    print("2/9 dimensionless groups")
    groups = table_groups()
    write_csv(output / "dimensionless_groups.csv", groups)

    print("3/9 grid convergence")
    grid = table_grid_convergence(verbose)
    write_csv(output / "grid_convergence.csv", grid)

    print("4/9 time-step convergence")
    timestep = table_timestep_convergence(verbose)
    write_csv(output / "timestep_convergence.csv", timestep)

    print("5/9 one-factor cases")
    cases, case_profiles = table_cases(output, verbose)
    write_csv(output / "case_summary.csv", cases)

    print("6/9 fully matched benchmark")
    matched, matched_profiles = table_matched(output, verbose)
    write_csv(output / "matched_benchmark.csv", matched)
    write_csv(output / "profiles.csv", case_profiles + matched_profiles)

    print("7/9 conservation and connectivity")
    conservation = table_conservation(verbose)
    write_csv(output / "conservation.csv", conservation)

    print("8/9 porosity and supersaturation profiles in time")
    snapshots = table_profile_snapshots(verbose)
    write_csv(output / "profile_snapshots.csv", snapshots)

    print("9/9 matched run carried to 6 h")
    long_run = table_long_run(output, verbose)
    write_csv(output / "matched_6h.csv", long_run)

    matched_b2 = next(
        r for r in matched if r["model"].startswith("B  P2-B") and "(B2)" in r["model"]
    )
    six_hour = min(long_run, key=lambda row: row["time_step_s"])
    summary = {
        "generated": date.today().isoformat(),
        "python": platform.python_version(),
        "wall_seconds": time.time() - started,
        "criteria": "config/acceptance_criteria_p2c.json",
        "reference_project": rc.REFERENCE,
        "headline": {
            "matched_B2_plateau_velocity_m_s": matched_b2["growth_velocity_m_s"],
            "matched_B2_plateau_difference_from_salt_grow_percent": matched_b2[
                "difference_from_C_percent"
            ],
            "matched_B2_6h_velocity_m_s": six_hour["growth_velocity_m_s"],
            "matched_B2_6h_difference_from_salt_grow_percent": six_hour[
                "velocity_difference_from_C_percent"
            ],
            "matched_B2_6h_minimum_porosity": six_hour["minimum_porosity"],
            "salt_grow_6h_minimum_porosity": REF_B2["minimum_porosity"],
            "sealing_predictions_correct": all(
                row["sealing_prediction_correct"] for row in cases
            ),
            "cases_that_lose_connectivity": [
                row["case"] for row in cases if row["disconnected"]
            ],
            "worst_salt_residual_relative": max(
                row["salt_residual_relative"] for row in cases + [matched_b2]
            ),
            "worst_water_residual_relative": max(
                row["water_residual_relative"] for row in cases + [matched_b2]
            ),
            "route_difference_corrected_worst": max(
                row["route_difference_corrected_final"] for row in conservation
            ),
            "route_difference_legacy_worst": max(
                row["route_difference_legacy_final"] for row in conservation
            ),
        },
        "tables": sorted(p.name for p in output.glob("*.csv")),
    }
    with (output / "summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True)
    print(json.dumps(summary["headline"], indent=2, sort_keys=True))
    print(f"\nwrote {len(summary['tables'])} tables to {output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
