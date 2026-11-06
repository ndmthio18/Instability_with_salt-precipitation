#!/usr/bin/env python3
"""Evaluate P1 V1--V7 against the frozen, predeclared tolerances."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Any


COMPONENTS = ("water", "air", "salt")


def rows(path: Path) -> list[dict[str, float]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return [{key: float(value) for key, value in row.items() if key}
                for row in csv.DictReader(handle)]


def independent_inventory(path: Path, projected_area: float) -> dict[str, float]:
    result = {component: 0.0 for component in COMPONENTS}
    with path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            volume = float(row["bulk_volume_m3"])
            porosity = float(row["porosity"])
            for phase in (0, 1):
                saturation = float(row[f"s{phase}"])
                molar_density = float(row[f"molar_density{phase}_mol_m3"])
                for component in COMPONENTS:
                    fraction = float(row[f"x{phase}_{component}"])
                    result[component] += volume * porosity * saturation * molar_density * fraction
    return {key: value / projected_area for key, value in result.items()}


def max_budget_residual(summary: dict[str, Any]) -> float:
    return max(summary[f"max_abs_{component}_residual_mol_m2"] for component in COMPONENTS)


def transition_time(case_dir: Path) -> float | None:
    data = rows(case_dir / "component_budget.csv")
    initial = tuple(int(data[0][key]) for key in
                    ("state_first_only", "state_both", "state_second_only"))
    for row in data[1:]:
        current = tuple(int(row[key]) for key in
                        ("state_first_only", "state_both", "state_second_only"))
        if current != initial:
            return row["time_s"]
    return None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("run_dir", type=Path)
    parser.add_argument("--cases", type=Path, default=Path(__file__).with_name("verification_cases.json"))
    parser.add_argument("--criteria", type=Path, default=Path(__file__).with_name("acceptance_criteria.json"))
    args = parser.parse_args()

    criteria = json.loads(args.criteria.read_text(encoding="utf-8"))
    manifest = json.loads(args.cases.read_text(encoding="utf-8"))
    raw = json.loads((args.run_dir / "verification_summary.json").read_text(encoding="utf-8"))
    summaries = {case["id"]: case for case in raw["cases"]}
    definitions = {case["id"]: case for case in manifest["cases"]}
    tolerance = criteria["absolute_component_budget_residual_mol_m2"]
    gates: dict[str, Any] = {}

    v1 = summaries["V1_closed_2d"]
    gates["V1"] = {"pass": v1["completed"] and max_budget_residual(v1) <= tolerance,
                   "max_abs_budget_residual_mol_m2": max_budget_residual(v1)}

    v2_cases = [summaries[f"V2_{component}_flux"] for component in COMPONENTS]
    gates["V2"] = {"pass": all(case["completed"] and max_budget_residual(case) <= tolerance
                                for case in v2_cases),
                   "max_abs_budget_residual_mol_m2": max(max_budget_residual(case) for case in v2_cases)}

    v3 = summaries["V3_hydrostatic_2d"]
    velocity_limit = (criteria["V3_max_spurious_velocity_fraction_of_evaporation"]
                      * manifest["reference_evaporation_velocity_m_s"])
    max_velocity = max(v3["max_liquid_darcy_velocity_m_s"], v3["max_gas_darcy_velocity_m_s"])
    gates["V3"] = {
        "pass": (v3["completed"]
                 and v3["max_pressure_endpoint_drift_pa"] <= criteria["V3_max_pressure_drift_pa"]
                 and max_velocity <= velocity_limit),
        "max_pressure_endpoint_drift_pa": v3["max_pressure_endpoint_drift_pa"],
        "max_phase_darcy_velocity_m_s": max_velocity,
        "velocity_limit_m_s": velocity_limit,
    }

    diffusion: dict[str, Any] = {}
    diffusivity = 1.5e-9 * 0.4057687912195745 ** (1.0 / 3.0)
    for grid in (40, 80):
        case_id = f"V4_diffusion_1d_n{grid}"
        data = rows(args.run_dir / case_id / "component_budget.csv")
        initial_amplitude = 0.5 * (data[0]["max_salt_x"] - data[0]["min_salt_x"])
        final_amplitude = 0.5 * (data[-1]["max_salt_x"] - data[-1]["min_salt_x"])
        expected = initial_amplitude * math.exp(-diffusivity * math.pi**2 * data[-1]["time_s"] / 0.20**2)
        diffusion[str(grid)] = {"initial_amplitude": initial_amplitude,
                                "final_amplitude": final_amplitude,
                                "analytical_amplitude": expected,
                                "relative_error": abs(final_amplitude - expected) / expected}
    coarse_error = diffusion["40"]["relative_error"]
    fine_error = diffusion["80"]["relative_error"]
    refinement_ratio = coarse_error / max(fine_error, 1e-30)
    gates["V4"] = {"pass": (summaries["V4_diffusion_1d_n40"]["completed"]
                             and summaries["V4_diffusion_1d_n80"]["completed"]
                             and coarse_error <= criteria["V4_relative_amplitude_error_n40"]
                             and refinement_ratio >= criteria["V4_refinement_ratio_required"]),
                   "effective_diffusivity_m2_s": diffusivity,
                   "results": diffusion,
                   "coarse_to_fine_error_ratio": refinement_ratio}

    one_d = rows(args.run_dir / "V5_uniform_1d" / "component_budget.csv")[-1]
    two_d = rows(args.run_dir / "V5_uniform_2d" / "component_budget.csv")[-1]
    compare_keys = [f"{component}_inventory_mol_m2" for component in COMPONENTS]
    compare_keys += ["min_salt_x", "max_salt_x", "min_sw", "max_sw",
                     "min_liquid_pressure_pa", "max_liquid_pressure_pa"]
    differences = {key: abs(one_d[key] - two_d[key]) / max(abs(one_d[key]), abs(two_d[key]), 1e-30)
                   for key in compare_keys}
    gates["V5"] = {"pass": (summaries["V5_uniform_1d"]["completed"]
                             and summaries["V5_uniform_2d"]["completed"]
                             and max(differences.values()) <= criteria["V5_max_1d_2d_relative_difference"]),
                   "max_relative_difference": max(differences.values()),
                   "relative_differences": differences}

    transitions = {}
    for case_id in ("V6_appearance_dt10", "V6_appearance_dt5",
                    "V6_disappearance_dt10", "V6_disappearance_dt5"):
        transitions[case_id] = transition_time(args.run_dir / case_id)
    v6_complete = all(summaries[key]["completed"] for key in transitions)
    v6_changed = all(value is not None for value in transitions.values())
    if v6_changed:
        appearance_difference = abs(transitions["V6_appearance_dt10"] - transitions["V6_appearance_dt5"])
        appearance_difference /= max(transitions["V6_appearance_dt5"], 1e-30)
        disappearance_difference = abs(transitions["V6_disappearance_dt10"] - transitions["V6_disappearance_dt5"])
        disappearance_difference /= max(transitions["V6_disappearance_dt5"], 1e-30)
    else:
        appearance_difference = math.inf
        disappearance_difference = math.inf
    finite_states = all(summaries[key]["min_sw"] >= criteria["saturation_lower_bound"]
                        and summaries[key]["max_sw"] <= criteria["saturation_upper_bound"]
                        and summaries[key]["min_salt_x"] >= criteria["mole_fraction_lower_bound"]
                        and summaries[key]["max_salt_x"] <= criteria["mole_fraction_upper_bound"]
                        for key in transitions)
    gates["V6"] = {"pass": (v6_complete and v6_changed and finite_states
                             and appearance_difference <= criteria["V6_dt_refinement_max_transition_time_difference_fraction"]
                             and disappearance_difference <= criteria["V6_dt_refinement_max_transition_time_difference_fraction"]),
                   "transition_times_s": transitions,
                   "appearance_dt_relative_difference": appearance_difference,
                   "disappearance_dt_relative_difference": disappearance_difference,
                   "physical_bounds_pass": finite_states}

    independent = {}
    for case_id in ("V1_closed_2d", "V2_water_flux", "V2_air_flux", "V2_salt_flux"):
        case_dir = args.run_dir / case_id
        metadata = json.loads((case_dir / "model_metadata.json").read_text(encoding="utf-8"))
        area = metadata["projected_area_m2_per_extrusion"]
        initial = independent_inventory(case_dir / "cell_state_initial.csv", area)
        final = independent_inventory(case_dir / "cell_state_final.csv", area)
        definition = definitions[case_id]["overrides"]
        expected_out = {component: 0.0 for component in COMPONENTS}
        if definition["Problem.BoundaryMode"] == "Flux":
            component = COMPONENTS[int(definition["Problem.TestComponent"])]
            expected_out[component] = float(definition["Problem.TestMolarFlux"]) * float(definition["TimeLoop.TEnd"])
        residual = {component: final[component] - initial[component] + expected_out[component]
                    for component in COMPONENTS}
        independent[case_id] = {"initial_mol_m2": initial, "final_mol_m2": final,
                                "analytic_boundary_out_mol_m2": expected_out,
                                "residual_mol_m2": residual,
                                "max_abs_residual_mol_m2": max(abs(value) for value in residual.values())}
    max_independent = max(value["max_abs_residual_mol_m2"] for value in independent.values())
    gates["V7"] = {"pass": max_independent <= tolerance,
                   "max_abs_independent_residual_mol_m2": max_independent,
                   "cases": independent}

    overall = all(gate["pass"] for gate in gates.values())
    result = {"overall_pass": overall, "gates": gates,
              "dependent_benchmarks_authorized": overall}
    output_path = args.run_dir / "gate_assessment.json"
    output_path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0 if overall else 1


if __name__ == "__main__":
    raise SystemExit(main())
