"""Compact, verify, and retain the reservoir-fed five-case matrix."""
import csv
import gzip
import hashlib
import json
import math
import pathlib
import shutil

root = pathlib.Path.cwd()
source = root / "matrix"
out = root / "analysis"
results = out / "results"
series = results / "time_series"
fields = results / "fields"
logs = out / "logs"
for path in (series, fields, logs):
    path.mkdir(parents=True, exist_ok=True)

names = ("kr0p1_av1", "kr1_av1_ref", "kr10_av1", "kr1_av0p5", "kr1_av2")


def floats(row):
    return {key: float(value) for key, value in row.items()}


def compact(path, destination, interval=100.0):
    first = last = activation = precipitation = None
    selected = []
    next_time = 0.0
    maxima = {"salt": 0.0, "water": 0.0, "solid": 0.0}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames
        for raw in reader:
            row = floats(raw)
            first = first or row
            last = row
            if activation is None and row["active"] == 1:
                activation = row
            if precipitation is None and row["top_solid_kg_m2"] > 0:
                precipitation = row
            maxima["salt"] = max(maxima["salt"], abs(row["salt_residual_kg"]))
            maxima["water"] = max(maxima["water"], abs(row["water_residual_kg"]))
            maxima["solid"] = max(maxima["solid"], abs(row["solid_ledger_residual_kg"]))
            if row["time_s"] + 1e-10 >= next_time or row is activation:
                selected.append(row)
                while next_time <= row["time_s"] + 1e-10:
                    next_time += interval
    if not last:
        raise RuntimeError(f"empty {path}")
    if selected[-1]["time_s"] != last["time_s"]:
        selected.append(last)
    with destination.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader(); writer.writerows(selected)
    return first, last, activation, precipitation, maxima, len(selected)


def compact_boundary(path, destination, interval=100.0):
    first = last = None
    selected = []
    next_time = 0.0
    max_local_error = 0.0
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames
        for raw in reader:
            row = floats(raw); first = first or row; last = row
            if row["dt_s"] > 0:
                liquid = (row["bottom_liquid_after_kg"] - row["bottom_liquid_before_kg"]) / row["dt_s"] - row["cut_liquid_down_kg_s"]
                salt = (row["bottom_salt_after_kg"] - row["bottom_salt_before_kg"]) / row["dt_s"] - row["cut_salt_down_kg_s"]
                max_local_error = max(max_local_error, abs(liquid-row["bottom_liquid_in_kg_s"]), abs(salt-row["bottom_salt_in_kg_s"]))
            if row["time_s"] + 1e-10 >= next_time:
                selected.append(row)
                while next_time <= row["time_s"] + 1e-10:
                    next_time += interval
    if selected[-1]["time_s"] != last["time_s"]:
        selected.append(last)
    with destination.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader(); writer.writerows(selected)
    return first, last, max_local_error


def compress(source_path, destination):
    with source_path.open("rb") as inp, gzip.open(destination, "wb", compresslevel=9) as output:
        shutil.copyfileobj(inp, output, 1024 * 1024)


summaries = []
for name in names:
    case = source / name
    spec = json.loads((case / "case_spec.json").read_text())
    command = json.loads((case / "command.json").read_text())
    first, last, active, precip, maxima, samples = compact(case / "coupled_history.csv", series / f"{name}.csv")
    _, bottom, local_error = compact_boundary(case / "open_boundary_budget.csv", series / f"{name}_boundary.csv")
    attempts = list(csv.DictReader((case / "donor_attempts.csv").open(newline="")))
    accepted = [row for row in attempts if int(row["sufficient"]) == 1]
    transfer = math.nan
    if accepted:
        row = accepted[-1]
        transfer = max(
            abs(float(row["salt_transfer_residual_kg"])) / max(abs(float(row["initial_salt_kg"])), 1e-30),
            abs(float(row["water_transfer_residual_kg"])) / max(abs(float(row["initial_water_kg"])), 1e-30),
        )
    config = next(csv.DictReader((case / "reservoir_config.csv").open(newline="")))
    completed = int((case / "exit_status.txt").read_text()) == 0 and last["time_s"] >= 179999.999
    passed = completed and active is not None and maxima["salt"] / abs(first["soil_salt_kg"]) <= 1e-8 and maxima["water"] / abs(first["soil_water_kg"]) <= 1e-8 and transfer <= 1e-12
    summaries.append({
        "case": name,
        "reaction_coefficient_m_s": spec["reaction_coefficient_m_s"],
        "av_multiplier": spec["av_multiplier"],
        "bottom_mode": spec["bottom_mode"],
        "water_level_m": 0.01,
        "initial_and_reservoir_salt_mole_fraction": float(config["configured_salt_mole_fraction"]),
        "activation_threshold_m": 1e-7,
        "soil_cells": spec["soil_cells"],
        "ale_cells": spec["ale_cells"],
        "dt_max_s": spec["dt_s"],
        "actual_end_s": last["time_s"],
        "completed_50h": int(completed),
        "ale_activated": int(active is not None),
        "status": "PASS" if passed else ("FAIL_NO_ALE_ACTIVATION" if completed and active is None else "FAIL"),
        "first_precipitation_s": precip["time_s"] if precip else math.nan,
        "activation_s": active["activation_s"] if active else math.nan,
        "final_height_m": last["h_m"],
        "final_mass_height_m": last["h_mass_m"],
        "final_epsilon_min": last["epsilon_min"],
        "final_top_solid_kg_m2": last["top_solid_kg_m2"],
        "final_internal_solid_kg_m2": last["internal_solid_kg_m2"],
        "final_internal_fraction": last["internal_fraction_net"],
        "final_evaporation_kg_m2_s": last["evaporation_kg_m2_s"],
        "cumulative_bottom_water_kg": bottom["cumulative_bottom_water_kg"],
        "cumulative_bottom_salt_kg": bottom["cumulative_bottom_salt_kg"],
        "salt_relative_residual": maxima["salt"] / abs(first["soil_salt_kg"]),
        "water_relative_residual": maxima["water"] / abs(first["soil_water_kg"]),
        "solid_ledger_abs_residual_kg": maxima["solid"],
        "activation_transfer_relative_residual": transfer,
        "bottom_local_reconstruction_abs_kg_s": local_error,
        "compact_rows": samples,
        "command": " ".join(command),
    })
    for filename in ("soil_fields.csv", "crust_fields.csv", "donor_profiles.csv"):
        compress(case / filename, fields / f"{name}_{filename}.gz")
    for filename in ("solver.out", "outcome.json", "exit_status.txt", "command.json", "case_spec.json", "params.input", "reservoir_config.csv", "activation_config.csv", "boundary_faces.csv", "coupling_events.csv", "donor_attempts.csv"):
        shutil.copy2(case / filename, logs / f"{name}_{filename}")

with (results / "summary.csv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(summaries[0])); writer.writeheader(); writer.writerows(summaries)

checks = []
def check(name, condition, detail):
    checks.append({"name": name, "passed": bool(condition), "detail": detail})
for row in summaries:
    check(f"{row['case']}: requested horizon", row["completed_50h"] == 1, row["actual_end_s"])
    check(f"{row['case']}: ALE activated", row["ale_activated"] == 1, row["activation_s"])
    check(f"{row['case']}: water conservation", row["water_relative_residual"] <= 1e-8, row["water_relative_residual"])
    check(f"{row['case']}: salt conservation", row["salt_relative_residual"] <= 1e-8, row["salt_relative_residual"])
    check(f"{row['case']}: activation transfer", row["activation_transfer_relative_residual"] <= 1e-12, row["activation_transfer_relative_residual"])
    check(f"{row['case']}: bottom local reconstruction", row["bottom_local_reconstruction_abs_kg_s"] <= 1e-16, row["bottom_local_reconstruction_abs_kg_s"])
verification = {"passed": all(item["passed"] for item in checks), "checks": checks}
(out / "verification.json").write_text(json.dumps(verification, indent=2) + "\n")

def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

provenance = {
    "branch": "research/p2-h-kr-av-reservoir",
    "base_commit": "1e713ae2a41e6a8ec8ba1a42f332cbb74714074c",
    "hpcc_root": str(root),
    "source_hashes": {
        path.name: sha256(path)
        for path in (
            root / "repo/model/p2_h_ale_coupling/soil_methods.inc",
            root / "repo/model/p2_h_ale_coupling/run_reservoir.py",
            root / "build/soil_solver",
        )
    },
}
(out / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
print(json.dumps({"cases": summaries, "verification": verification}, indent=2))
raise SystemExit(0 if verification["passed"] else 1)
