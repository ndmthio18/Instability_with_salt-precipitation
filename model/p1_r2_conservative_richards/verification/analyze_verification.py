#!/usr/bin/env python3
"""Analyze P1-R2 V1-V3 solver budgets against preregistered gates."""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import defaultdict
from pathlib import Path


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def max_abs(rows: list[dict[str, str]], key: str) -> float:
    return max(abs(float(row[key])) for row in rows)


def analyze_run(
    name: str,
    run: Path,
    criteria: dict[str, object],
) -> dict[str, object]:
    budget = read_csv(run / "component_budget.csv")
    faces = read_csv(run / "independent_face_fluxes.csv")
    cuts = read_csv(run / "cut_fluxes.csv")
    steps = read_csv(run / "accepted_steps.csv")
    if not budget or not steps:
        raise ValueError(f"{name}: missing accepted budget records")

    initial_total = abs(float(budget[0]["total_mass_kg"]))
    initial_salt = abs(float(budget[0]["salt_mass_kg"]))
    initial_water = abs(float(budget[0]["water_mass_kg"]))
    absolute = float(criteria["controlled_test_budget_absolute_kg"])
    relative = float(criteria["controlled_test_budget_relative"])
    tolerances = {
        "total": max(absolute, relative * initial_total),
        "salt": max(absolute, relative * initial_salt),
        "water": max(absolute, relative * initial_water),
    }
    residuals = {
        "total": max_abs(budget, "total_residual_kg"),
        "salt": max_abs(budget, "salt_residual_kg"),
        "water": max_abs(budget, "water_residual_kg"),
    }
    budget_pass = all(
        residuals[key] <= tolerances[key] for key in residuals
    )

    boundary_sum: dict[int, list[float]] = defaultdict(
        lambda: [0.0, 0.0]
    )
    cut_sum: dict[tuple[int, str], list[float]] = defaultdict(
        lambda: [0.0, 0.0]
    )
    for row in faces:
        step = int(row["step"])
        if row["kind"] == "boundary":
            boundary_sum[step][0] += float(row["total_flux_kg_s"])
            boundary_sum[step][1] += float(row["salt_flux_kg_s"])
        elif row["kind"] == "internal_cut":
            key = (step, row["cut"])
            cut_sum[key][0] += float(row["total_flux_kg_s"])
            cut_sum[key][1] += float(row["salt_flux_kg_s"])

    face_error = 0.0
    for row in budget:
        step = int(row["step"])
        reconstructed = boundary_sum[step]
        face_error = max(
            face_error,
            abs(
                reconstructed[0]
                - float(row["external_total_out_kg_s"])
            ),
            abs(
                reconstructed[1]
                - float(row["external_salt_out_kg_s"])
            ),
        )
    for row in cuts:
        key = (int(row["step"]), row["cut"])
        reconstructed = cut_sum[key]
        face_error = max(
            face_error,
            abs(reconstructed[0] - float(row["total_down_kg_s"])),
            abs(reconstructed[1] - float(row["salt_down_kg_s"])),
        )
    flux_scale = max(
        1e-30,
        max(
            abs(float(row["external_total_out_kg_s"]))
            + abs(float(row["external_salt_out_kg_s"]))
            for row in budget
        ),
    )
    face_tolerance = max(
        1e-14,
        float(criteria["face_reconstruction_relative"]) * flux_scale,
    )
    face_pass = face_error <= face_tolerance

    result: dict[str, object] = {
        "run": name,
        "accepted_steps": len(steps),
        "final_time_s": float(budget[-1]["time_s"]),
        "max_budget_residual_kg": residuals,
        "budget_tolerance_kg": tolerances,
        "budget_pass": budget_pass,
        "max_face_reconstruction_error_kg_s": face_error,
        "face_reconstruction_tolerance_kg_s": face_tolerance,
        "face_reconstruction_pass": face_pass,
        "retry_count": sum(int(row["retried"]) for row in steps),
    }

    if name == "V2_component_flux":
        expected_total = 0.06 * 1e-5
        expected_salt = 0.06 * 1e-7
        observed_total = max_abs(budget[1:], "external_total_out_kg_s")
        observed_salt = max_abs(budget[1:], "external_salt_out_kg_s")
        flux_relative = float(criteria["prescribed_flux_relative"])
        flux_pass = (
            math.isclose(
                observed_total,
                expected_total,
                rel_tol=flux_relative,
                abs_tol=1e-14,
            )
            and math.isclose(
                observed_salt,
                expected_salt,
                rel_tol=flux_relative,
                abs_tol=1e-14,
            )
        )
        result["prescribed_flux"] = {
            "expected_total_kg_s": expected_total,
            "observed_total_kg_s": observed_total,
            "expected_salt_kg_s": expected_salt,
            "observed_salt_kg_s": observed_salt,
            "pass": flux_pass,
        }
    result["pass"] = bool(
        budget_pass
        and face_pass
        and (
            name != "V2_component_flux"
            or result["prescribed_flux"]["pass"]
        )
    )
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--runs", type=Path, required=True)
    parser.add_argument("--criteria", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    criteria = json.loads(args.criteria.read_text(encoding="utf-8"))
    results = [
        analyze_run(name, args.runs / name, criteria)
        for name in (
            "V1_closed_box",
            "V2_component_flux",
            "V3_variable_density",
        )
    ]
    payload = {
        "schema_version": 1,
        "gates": results,
        "all_pass": all(item["pass"] for item in results),
    }
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps(payload, indent=2, sort_keys=True))
    if not payload["all_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
