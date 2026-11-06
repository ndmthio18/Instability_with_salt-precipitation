#!/usr/bin/env python3
"""Run the staged P1 verification cases without changing acceptance criteria."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
from typing import Any


COMPONENTS = ("water", "air", "salt")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def command_for(executable: Path, base_input: Path, overrides: dict[str, Any]) -> list[str]:
    command = [str(executable), str(base_input)]
    for key, value in overrides.items():
        if isinstance(value, bool):
            value = "true" if value else "false"
        command.extend([f"-{key}", str(value)])
    return command


def parse_budget(path: Path) -> list[dict[str, float]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return [{key: float(value) for key, value in row.items()}
                for row in csv.DictReader(handle)]


def summarize_case(case: dict[str, Any], case_dir: Path, returncode: int) -> dict[str, Any]:
    summary: dict[str, Any] = {"id": case["id"], "returncode": returncode, "completed": False}
    status_path = case_dir / "run_status.json"
    if status_path.exists():
        summary["run_status"] = json.loads(status_path.read_text(encoding="utf-8"))
        summary["completed"] = summary["run_status"].get("status") == "completed"
    budget_path = case_dir / "component_budget.csv"
    if not budget_path.exists():
        return summary
    rows = parse_budget(budget_path)
    summary["accepted_rows"] = len(rows) - 1
    summary["final_time_s"] = rows[-1]["time_s"]
    for component in COMPONENTS:
        residual_key = f"{component}_residual_mol_m2"
        inventory_key = f"{component}_inventory_mol_m2"
        out_key = f"{component}_cumulative_out_mol_m2"
        summary[f"max_abs_{residual_key}"] = max(abs(row[residual_key]) for row in rows)
        summary[f"inventory_change_{component}_mol_m2"] = rows[-1][inventory_key] - rows[0][inventory_key]
        summary[f"cumulative_out_{component}_mol_m2"] = rows[-1][out_key]
    summary["initial_pressure_range_pa"] = [rows[0]["min_liquid_pressure_pa"], rows[0]["max_liquid_pressure_pa"]]
    summary["final_pressure_range_pa"] = [rows[-1]["min_liquid_pressure_pa"], rows[-1]["max_liquid_pressure_pa"]]
    summary["max_pressure_endpoint_drift_pa"] = max(
        abs(rows[-1]["min_liquid_pressure_pa"] - rows[0]["min_liquid_pressure_pa"]),
        abs(rows[-1]["max_liquid_pressure_pa"] - rows[0]["max_liquid_pressure_pa"]),
    )
    summary["max_liquid_darcy_velocity_m_s"] = max(row.get("max_liquid_darcy_velocity_m_s", 0.0) for row in rows)
    summary["max_gas_darcy_velocity_m_s"] = max(row.get("max_gas_darcy_velocity_m_s", 0.0) for row in rows)
    summary["initial_state_counts"] = [int(rows[0][key]) for key in
                                       ("state_first_only", "state_both", "state_second_only")]
    summary["final_state_counts"] = [int(rows[-1][key]) for key in
                                     ("state_first_only", "state_both", "state_second_only")]
    summary["min_sw"] = min(row["min_sw"] for row in rows)
    summary["max_sw"] = max(row["max_sw"] for row in rows)
    summary["min_salt_x"] = min(row["min_salt_x"] for row in rows)
    summary["max_salt_x"] = max(row["max_salt_x"] for row in rows)
    summary["final_top_salt_mean"] = rows[-1]["top_salt_mean"]
    return summary


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--binary-1d", type=Path, required=True)
    parser.add_argument("--binary-2d", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--cases", type=Path, default=Path(__file__).with_name("verification_cases.json"))
    parser.add_argument("--criteria", type=Path, default=Path(__file__).with_name("acceptance_criteria.json"))
    args = parser.parse_args()

    cases_path = args.cases.resolve()
    criteria_path = args.criteria.resolve()
    manifest = json.loads(cases_path.read_text(encoding="utf-8"))
    criteria = json.loads(criteria_path.read_text(encoding="utf-8"))
    base_input = (cases_path.parent / manifest["base_input"]).resolve()
    binaries = {"1d": args.binary_1d.resolve(), "2d": args.binary_2d.resolve()}
    args.output.mkdir(parents=True, exist_ok=False)

    provenance = {
        "cases_sha256": sha256(cases_path),
        "criteria_sha256": sha256(criteria_path),
        "base_input_sha256": sha256(base_input),
        "binary_1d_sha256": sha256(binaries["1d"]),
        "binary_2d_sha256": sha256(binaries["2d"]),
        "environment": {key: os.environ.get(key) for key in
                        ("OMPI_MCA_pml", "OMPI_MCA_btl", "OMPI_MCA_mtl", "UCX_TLS")},
    }
    (args.output / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n", encoding="utf-8")

    summaries = []
    for case in manifest["cases"]:
        case_dir = args.output / case["id"]
        case_dir.mkdir()
        command = command_for(binaries[case["target"]], base_input, case["overrides"])
        (case_dir / "command.json").write_text(json.dumps(command, indent=2) + "\n", encoding="utf-8")
        with (case_dir / "run.log").open("w", encoding="utf-8") as log:
            completed = subprocess.run(command, cwd=case_dir, stdout=log, stderr=subprocess.STDOUT, check=False)
        summaries.append(summarize_case(case, case_dir, completed.returncode))

    result = {"criteria": criteria, "cases": summaries}
    (args.output / "verification_summary.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, indent=2))
    return 0 if all(case["completed"] for case in summaries) else 1


if __name__ == "__main__":
    sys.exit(main())
