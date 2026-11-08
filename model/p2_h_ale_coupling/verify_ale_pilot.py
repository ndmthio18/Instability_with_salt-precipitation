"""Conservative ALE grid/time pilot audit at a common accepted time."""
import csv
import json
import math
import pathlib

import numpy as np

root = pathlib.Path.cwd()
reference = root / "pilot/pilot_ref_ale64"
cases = {
    "reference_ale64_dt0p5": reference,
    "grid_ale128_dt0p5": root / "ale_pilot/pilot_grid_ale128",
    "ale64_dt0p25": root / "ale_pilot/pilot_dt_half_ale64",
}


def rows(path):
    with path.open(newline="") as handle:
        return [{key: float(value) for key, value in row.items()} for row in csv.DictReader(handle)]


histories = {name: rows(path / "coupled_history.csv") for name, path in cases.items()}
target = min(history[-1]["time_s"] for history in histories.values())


def at(history, time):
    return min(history, key=lambda row: abs(row["time_s"] - time))


metrics = ("h_m", "h_mass_m", "internal_solid_kg_m2", "top_solid_kg_m2", "epsilon_min", "evaporation_kg_m2_s")
states = {name: at(history, target) for name, history in histories.items()}
comparisons = {}
for name in ("grid_ale128_dt0p5", "ale64_dt0p25"):
    comparisons[name] = {
        metric: abs(states[name][metric] - states["reference_ale64_dt0p5"][metric])
        / max(abs(states["reference_ale64_dt0p5"][metric]), 1e-30)
        for metric in metrics
    }


def profiles(path):
    grouped = {}
    for row in rows(path / "crust_fields.csv"):
        grouped.setdefault(row["time_s"], []).append(row)
    return grouped


profile_sets = {name: profiles(path) for name, path in cases.items()}
common_times = set.intersection(*(set(grouped) for grouped in profile_sets.values()))
profile_time = max(time for time in common_times if time <= target)
eta = np.linspace(0, 1, 400)
profile_values = {}
for name, grouped in profile_sets.items():
    profile = grouped[profile_time]
    center = np.array([0.5 * (row["z_bottom_m"] + row["z_top_m"]) for row in profile])
    height = max(row["z_top_m"] for row in profile)
    epsilon = np.array([row["epsilon"] for row in profile])
    profile_values[name] = np.interp(eta, center / height, epsilon, left=epsilon[0], right=epsilon[-1])
profile_errors = {
    name: {
        "l1": float(np.mean(np.abs(values-profile_values["reference_ale64_dt0p5"]))),
        "linf": float(np.max(np.abs(values-profile_values["reference_ale64_dt0p5"]))),
    }
    for name, values in profile_values.items() if name != "reference_ale64_dt0p5"
}

checks = []
for name, path in cases.items():
    history = histories[name]
    outcome = json.loads((path / "outcome.json").read_text()) if (path / "outcome.json").exists() else {"exit_code": 0}
    initial_salt = abs(history[0]["soil_salt_kg"]); initial_water = abs(history[0]["soil_water_kg"])
    salt = max(abs(row["salt_residual_kg"]) for row in history) / initial_salt
    water = max(abs(row["water_residual_kg"]) for row in history) / initial_water
    checks.extend([
        {"name": f"{name}: completed common window", "passed": outcome["exit_code"] == 0 and history[-1]["time_s"] >= target-1e-8, "detail": history[-1]["time_s"]},
        {"name": f"{name}: activated", "passed": any(row["active"] == 1 for row in history), "detail": states[name]["activation_s"]},
        {"name": f"{name}: water balance", "passed": water <= 1e-8, "detail": water},
        {"name": f"{name}: salt balance", "passed": salt <= 1e-8, "detail": salt},
    ])
result = {
    "passed": all(item["passed"] for item in checks),
    "common_time_s": target,
    "profile_time_s": profile_time,
    "states": states,
    "relative_scalar_differences": comparisons,
    "porosity_profile_errors": profile_errors,
    "interpretation": "Diagnostic only. In particular, epsilon_min and local clogging time are not declared converged at ALE64.",
    "checks": checks,
}
print(json.dumps(result, indent=2))
raise SystemExit(0 if result["passed"] else 1)
