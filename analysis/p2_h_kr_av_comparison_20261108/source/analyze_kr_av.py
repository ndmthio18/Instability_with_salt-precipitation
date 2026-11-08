"""Stream and compact the 50 h P2-H k_r/a_v comparison evidence."""
import csv
import gzip
import hashlib
import json
import math
import pathlib
import shutil

root = pathlib.Path.cwd()
source = root / "production_corrected"
analysis = root / "analysis"
results = analysis / "results"
time_series = results / "time_series"
fields = results / "fields"
logs = analysis / "logs"
for path in (time_series, fields, logs):
    path.mkdir(parents=True, exist_ok=True)

expected = (
    "kr0p1_av1",
    "kr1_av1_ref",
    "kr10_av1",
    "kr1_av0p5",
    "kr1_av2",
)

def numeric(row):
    return {key: float(value) for key, value in row.items()}

def mole_fraction_from_mass_fraction(mass_fraction):
    molar_mass_salt = 58.44e-3
    molar_mass_water = 18.01528e-3
    salt_moles = mass_fraction / molar_mass_salt
    water_moles = (1.0 - mass_fraction) / molar_mass_water
    return salt_moles / (salt_moles + water_moles)

def compact_history(path, destination):
    first = None
    last = None
    activation = None
    first_precipitation = None
    sampled = []
    next_sample = 0.0
    max_salt = 0.0
    max_water = 0.0
    max_ledger = 0.0
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        fields_in = reader.fieldnames
        for raw in reader:
            row = numeric(raw)
            if first is None:
                first = row
            last = row
            if activation is None and row["active"] == 1:
                activation = row
            if first_precipitation is None and row["top_solid_kg_m2"] > 0:
                first_precipitation = row
            max_salt = max(max_salt, abs(row["salt_residual_kg"]))
            max_water = max(max_water, abs(row["water_residual_kg"]))
            max_ledger = max(max_ledger, abs(row["solid_ledger_residual_kg"]))
            if row["time_s"] + 1e-10 >= next_sample or (
                activation is row
            ):
                sampled.append(row)
                while next_sample <= row["time_s"] + 1e-10:
                    next_sample += 100.0
    if first is None:
        raise RuntimeError(f"empty history {path}")
    if sampled[-1]["time_s"] != last["time_s"]:
        sampled.append(last)
    with destination.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields_in)
        writer.writeheader()
        writer.writerows(sampled)
    return first, last, activation, first_precipitation, max_salt, max_water, max_ledger, len(sampled)

def event_reason(path):
    if not path.exists():
        return "missing coupling_events.csv"
    last = None
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            last = row
    return last["event"] if last else "no recorded event"

def compress_file(source_path, destination):
    with source_path.open("rb") as source_handle, gzip.open(destination, "wb", compresslevel=9) as target:
        shutil.copyfileobj(source_handle, target, length=1024 * 1024)

def field_times(path):
    times = []
    previous = None
    rows = 0
    with path.open(newline="") as handle:
        for raw in csv.DictReader(handle):
            time = float(raw["time_s"])
            rows += 1
            if previous is None or time != previous:
                times.append(time)
                previous = time
    return times, rows

summaries = []
histories = {}
for name in expected:
    case = source / name
    if not case.is_dir():
        raise FileNotFoundError(case)
    spec = json.loads((case / "case_spec.json").read_text())
    command = json.loads((case / "command.json").read_text())
    first, last, activation, first_precipitation, max_salt, max_water, max_ledger, samples = compact_history(
        case / "coupled_history.csv", time_series / f"{name}.csv"
    )
    with (time_series / f"{name}.csv").open(newline="") as handle:
        histories[name] = [numeric(row) for row in csv.DictReader(handle)]
    exit_code = int((case / "exit_status.txt").read_text())
    completed = exit_code == 0 and last["time_s"] >= 179999.999
    last_coupling_event = event_reason(case / "coupling_events.csv")
    stop_reason = "requested_horizon" if completed else last_coupling_event
    scientific_status = (
        "PASS" if completed and activation is not None else
        "FAIL_NO_ALE_ACTIVATION" if completed else
        "FAIL_EARLY_STOP"
    )
    soil_times, soil_rows = field_times(case / "soil_fields.csv")
    crust_times, crust_rows = field_times(case / "crust_fields.csv")
    compress_file(case / "soil_fields.csv", fields / f"{name}_soil_fields.csv.gz")
    compress_file(case / "crust_fields.csv", fields / f"{name}_crust_fields.csv.gz")

    transfer_relative = math.nan
    attempts = case / "donor_attempts.csv"
    if attempts.exists():
        successful = []
        with attempts.open(newline="") as handle:
            successful = [row for row in csv.DictReader(handle) if int(row["sufficient"]) == 1]
        if successful:
            row = successful[-1]
            transfer_relative = max(
                abs(float(row["salt_transfer_residual_kg"])) / float(row["initial_salt_kg"]),
                abs(float(row["water_transfer_residual_kg"])) / float(row["initial_water_kg"]),
            )

    summaries.append(
        {
            "case": name,
            "reaction_coefficient_m_s": spec["reaction_coefficient_m_s"],
            "av_multiplier": spec["av_multiplier"],
            "av_closure": "alpha_av*36*(1-epsilon)/d_sub; d_sub=25e-6 m",
            "initial_crystal_interpretation": "zero solid; conservative activation at h_eq=0.1 um",
            "initial_salt_mole_fraction": mole_fraction_from_mass_fraction(first["top_soil_mass_fraction"]),
            "activation_threshold_m": 1e-7,
            "soil_cells": spec["soil_cells"],
            "ale_cells": spec["ale_cells"],
            "maximum_soil_timestep_s": spec["dt_s"],
            "actual_end_s": last["time_s"],
            "completed_50h": int(completed),
            "ale_activated": int(activation is not None),
            "scientific_status": scientific_status,
            "exit_code": exit_code,
            "stop_reason": stop_reason,
            "last_coupling_event": last_coupling_event,
            "first_surface_precipitation_s": first_precipitation["time_s"] if first_precipitation else math.nan,
            "activation_s": activation["activation_s"] if activation else math.nan,
            "activation_initial_crust_concentration_mass_fraction": activation["crust_top_mass_fraction"] if activation else math.nan,
            "final_crust_height_m": last["h_m"],
            "final_minimum_porosity": last["epsilon_min"],
            "final_evaporation_kg_m2_s": last["evaporation_kg_m2_s"],
            "total_evaporated_water_kg_m2": last["vapour_kg"] / 0.001,
            "final_internal_solid_kg_m2": last["internal_solid_kg_m2"],
            "final_top_new_solid_kg_m2": last["top_new_solid_kg_m2"],
            "final_total_top_solid_kg_m2": last["top_solid_kg_m2"],
            "final_internal_fraction_net": last["internal_fraction_net"],
            "salt_relative_residual": max_salt / max(abs(first["soil_salt_kg"]), 1e-300),
            "water_relative_residual": max_water / max(abs(first["soil_water_kg"]), 1e-300),
            "solid_ledger_absolute_residual_kg": max_ledger,
            "activation_transfer_relative_residual": transfer_relative,
            "soil_field_snapshots": len(soil_times),
            "crust_field_snapshots": len(crust_times),
            "soil_field_max_gap_s": max((b - a for a, b in zip(soil_times, soil_times[1:])), default=math.nan),
            "crust_field_max_gap_s": max((b - a for a, b in zip(crust_times, crust_times[1:])), default=math.nan),
            "soil_field_rows": soil_rows,
            "crust_field_rows": crust_rows,
            "compact_history_rows": samples,
            "command": " ".join(command),
        }
    )

summary_path = results / "summary.csv"
with summary_path.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(summaries[0]))
    writer.writeheader()
    writer.writerows(summaries)

(results / "case_matrix.json").write_text(
    json.dumps(
        {
            "reference": {"reaction_coefficient_m_s": 1e-3, "av_multiplier": 1.0},
            "fixed_av_group": ["kr0p1_av1", "kr1_av1_ref", "kr10_av1"],
            "fixed_kr_group": ["kr1_av0p5", "kr1_av1_ref", "kr1_av2"],
            "unique_runs": list(expected),
            "reference_reused": True,
        },
        indent=2,
    )
    + "\n"
)

def interpolate(history, target):
    if target < history[0]["time_s"] or target > history[-1]["time_s"]:
        return None
    for left, right in zip(history, history[1:]):
        if left["time_s"] <= target <= right["time_s"]:
            if right["time_s"] == left["time_s"]:
                return left
            weight = (target - left["time_s"]) / (right["time_s"] - left["time_s"])
            return {key: left[key] + weight * (right[key] - left[key]) for key in left}
    return history[-1]

comparison_fields = (
    "h_m",
    "epsilon_min",
    "evaporation_kg_m2_s",
    "internal_solid_kg_m2",
    "top_new_solid_kg_m2",
    "internal_fraction_net",
    "top_soil_mass_fraction",
    "top_soil_sw",
)
comparisons = []
for group, names in (
    ("fixed_av", ("kr0p1_av1", "kr1_av1_ref", "kr10_av1")),
    ("fixed_kr", ("kr1_av0p5", "kr1_av1_ref", "kr1_av2")),
):
    common_end = min(histories[name][-1]["time_s"] for name in names)
    targets = {
        "initial": 0.0,
        "quarter_common_window": 0.25 * common_end,
        "half_common_window": 0.50 * common_end,
        "three_quarter_common_window": 0.75 * common_end,
        "common_accepted_end": common_end,
    }
    group_activation_times = [
        row["activation_s"]
        for name in names
        for row in histories[name]
        if row["active"] == 1
    ]
    if group_activation_times:
        first_activation = min(group_activation_times)
        targets["first_group_activation"] = first_activation
        if first_activation + 100.0 <= common_end:
            targets["100s_after_first_activation"] = first_activation + 100.0
    for label, target in sorted(targets.items(), key=lambda item: item[1]):
        for name in names:
            row = interpolate(histories[name], target)
            if row is not None:
                comparisons.append(
                    {"group": group, "comparison_time": label, "time_s": target, "case": name,
                     **{field: row[field] for field in comparison_fields}}
                )
with (results / "common_time_comparison.csv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(comparisons[0]))
    writer.writeheader()
    writer.writerows(comparisons)

spreads = []
for group in ("fixed_av", "fixed_kr"):
    labels = sorted({row["comparison_time"] for row in comparisons if row["group"] == group})
    for label in labels:
        selected = [row for row in comparisons if row["group"] == group and row["comparison_time"] == label]
        for field in comparison_fields:
            values = [row[field] for row in selected]
            mean = sum(values) / len(values)
            spread = max(values) - min(values)
            spreads.append({
                "group": group,
                "comparison_time": label,
                "time_s": selected[0]["time_s"],
                "quantity": field,
                "absolute_range": spread,
                "relative_range": spread / abs(mean) if mean != 0 else (0.0 if spread == 0 else math.inf),
            })
with (results / "comparison_spread.csv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(spreads[0]))
    writer.writeheader()
    writer.writerows(spreads)

def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

provenance = {
    "branch": "research/p2-h-kr-av-comparison",
    "base_commit": "63f267f9eae85bf434a2ad9ed0804a63723e9478",
    "hpcc_root": "/WORK/hulm_work/zhixin/p2_h_kr_av_comparison_20261108",
    "build_job": "16407548",
    "corrected_matrix_job": "16408105",
    "invalidated_wrong_namespace_matrix_job": "16407598",
    "solver_sha256": sha256(root / "build/soil_solver"),
    "case_jobs": sorted({json.loads((source / name / "outcome.json").read_text())["slurm_job_id"] for name in expected}),
    "interpretation": "zero initial crust and solid; conservative ALE activation threshold at 0.1 um equivalent real surface solid",
    "initial_salt_mole_fraction": 0.05,
    "requested_end_s": 180000.0,
    "source_hashes": {
        name: sha256(root / "repo/model/p2_h_ale_coupling" / name)
        for name in ("run_kr_av.py", "analyze_kr_av.py", "plot_kr_av.py", "verify_kr_av.py")
    },
}
(analysis / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
print(summary_path)
