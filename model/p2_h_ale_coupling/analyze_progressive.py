"""Postprocess progressive-activation evidence on a Slurm compute node."""
import csv
import hashlib
import json
import math
import os
import pathlib

if not os.environ.get("SLURM_JOB_ID"):
    raise RuntimeError("analyze_progressive.py requires a Slurm compute allocation")

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt


ROOT = pathlib.Path.cwd()
OUT = ROOT / "analysis_output"
FIG = OUT / "figures"
FIG.mkdir(parents=True, exist_ok=True)


def read_csv(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def f(row, key, default=math.nan):
    try:
        return float(row[key])
    except (KeyError, TypeError, ValueError):
        return default


def nearest(rows, target):
    return min(rows, key=lambda row: abs(f(row, "time_s") - target))


def initial_crust_c(case_dir, activation):
    rows = read_csv(case_dir / "crust_fields.csv")
    selected = [r for r in rows if abs(f(r, "time_s") - activation) < 1e-7]
    if not selected:
        return math.nan
    salt = sum(f(r, "concentration") * (f(r, "z_top_m") - f(r, "z_bottom_m")) * f(r, "epsilon") for r in selected)
    pore = sum((f(r, "z_top_m") - f(r, "z_bottom_m")) * f(r, "epsilon") for r in selected)
    return salt / pore


def final_epsilon_max(case_dir, final_time):
    path = case_dir / "ale_diagnostics.csv"
    if not path.exists():
        return math.nan
    latest = -math.inf
    maximum = math.nan
    with path.open(newline="") as stream:
        for row in csv.DictReader(stream):
            t = f(row, "time_s")
            if row["status"] != "accepted" or t > final_time + 1e-7:
                continue
            if t > latest + 1e-10:
                latest, maximum = t, f(row, "epsilon")
            elif abs(t-latest) <= 1e-10:
                maximum = max(maximum, f(row, "epsilon"))
    return maximum


records = []
histories = {}
for case_dir in sorted((ROOT / "runs").glob("*")):
    if not (case_dir / "case_spec.json").exists() or not (case_dir / "coupled_history.csv").exists():
        continue
    spec = json.loads((case_dir / "case_spec.json").read_text())
    outcome = json.loads((case_dir / "outcome.json").read_text()) if (case_dir / "outcome.json").exists() else {"exit_code": 999}
    rows = read_csv(case_dir / "coupled_history.csv")
    if not rows:
        continue
    histories[spec["name"]] = rows
    precip = next((r for r in rows if f(r, "precipitation_flux_kg_m2_s", 0) > 0), None)
    active = next((r for r in rows if int(float(r["active"])) == 1), None)
    pre = None
    if active:
        before = [r for r in rows if f(r, "time_s") < f(active, "time_s") and int(float(r["active"])) == 0]
        pre = before[-1] if before else None
    activation = f(active, "activation_s") if active else math.nan
    init_c = initial_crust_c(case_dir, activation) if active else math.nan
    post500 = nearest([r for r in rows if active and f(r, "time_s") >= activation], activation + 500) if active else None
    first = rows[0]
    initial_salt = abs(f(first, "soil_salt_kg"))
    initial_water = abs(f(first, "soil_water_kg"))
    max_salt_rel = max(abs(f(r, "salt_residual_kg")) for r in rows) / max(initial_salt, 1e-30)
    max_water_rel = max(abs(f(r, "water_residual_kg")) for r in rows) / max(initial_water, 1e-30)
    end = rows[-1]
    record = {
        "case": spec["name"], "mode": spec["donor_mode"],
        "soil_cells": spec["soil_cells"], "dt_s": spec["dt_s"], "ale_cells": spec["ale_cells"],
        "exit_code": outcome.get("exit_code", 999), "first_precipitation_s": f(precip, "time_s") if precip else math.nan,
        "activation_s": activation, "preactivation_top_c": f(pre, "top_soil_mass_fraction") if pre else math.nan,
        "preactivation_top_sw": f(pre, "top_soil_sw") if pre else math.nan,
        "initial_ale_c": init_c, "initial_concentration_jump": init_c - f(pre, "top_soil_mass_fraction") if pre else math.nan,
        "precursor_pore_salt_kg_m2": f(active, "precursor_pore_salt_kg_m2") if active else math.nan,
        "precursor_pore_water_kg_m2": f(active, "precursor_pore_water_kg_m2") if active else math.nan,
        "precursor_required_volume_m3_m2": f(active, "precursor_required_volume_m3_m2") if active else math.nan,
        "precursor_filled_volume_m3_m2": f(active, "precursor_filled_volume_m3_m2") if active else math.nan,
        "precursor_unfilled_volume_m3_m2": f(active, "precursor_unfilled_volume_m3_m2") if active else math.nan,
        "post500_time_s": f(post500, "time_s") if post500 else math.nan,
        "post500_height_m": f(post500, "h_m") if post500 else math.nan,
        "post500_epsilon_min": f(post500, "epsilon_min") if post500 else math.nan,
        "post500_internal_positive_kg_m2": f(post500, "internal_positive_kg_m2") if post500 else math.nan,
        "post500_internal_dissolution_kg_m2": f(post500, "internal_dissolution_kg_m2") if post500 else math.nan,
        "final_time_s": f(end, "time_s"), "final_height_m": f(end, "h_m"),
        "final_epsilon_min": f(end, "epsilon_min"), "final_epsilon_max": final_epsilon_max(case_dir, f(end, "time_s")),
        "final_internal_positive_kg_m2": f(end, "internal_positive_kg_m2"),
        "final_internal_dissolution_kg_m2": f(end, "internal_dissolution_kg_m2"),
        "max_salt_relative_residual": max_salt_rel, "max_water_relative_residual": max_water_rel,
        "completed_requested_horizon": int(outcome.get("exit_code", 999) == 0 and f(end, "time_s") >= spec["tend_s"] - 1e-6),
    }
    records.append(record)


def write_table(name, rows, fields=None):
    if fields is None:
        fields = list(rows[0]) if rows else []
    with (OUT / name).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader(); writer.writerows(rows)


activation = [r for r in records if r["soil_cells"] == 80 and r["dt_s"] == 0.5 and r["ale_cells"] == 64 and not r["case"].startswith("pilot_")]
write_table("ACTIVATION_COMPARISON.csv", activation)
convergence = [r for r in records if r["mode"] == "progressive_top_boundary_filling" and r["ale_cells"] == 64 and not r["case"].startswith("pilot_")]
write_table("GRID_CONVERGENCE.csv", convergence)
write_table("TIMESTEP_CONVERGENCE.csv", convergence)

budget_rows = []
for record in records:
    rows = histories[record["case"]]
    for index, row in enumerate(rows):
        if index == len(rows)-1 or index % max(1, round(100 / record["dt_s"])) == 0:
            budget_rows.append({
                "case": record["case"], "time_s": f(row, "time_s"),
                "salt_residual_kg": f(row, "salt_residual_kg"),
                "water_residual_kg": f(row, "water_residual_kg"),
                "soil_salt_kg": f(row, "soil_salt_kg"), "crust_salt_kg": f(row, "crust_salt_kg"),
                "soil_water_kg": f(row, "soil_water_kg"), "crust_water_kg": f(row, "crust_water_kg"),
                "vapour_kg": f(row, "vapour_kg"),
            })
write_table("WATER_SALT_BUDGET.csv", budget_rows)


plt.rcParams.update({"font.size": 8, "axes.linewidth": .8, "savefig.dpi": 300})
pilot_name = "pilot_progressive_n80_dt0p5_a64"
pilot = histories.get(pilot_name)
if pilot:
    t = [f(r, "time_s") / 1000 for r in pilot]
    fig, ax = plt.subplots(figsize=(3.5, 2.6))
    ax.plot(t, [f(r, "h_mass_m") * 1e6 for r in pilot], label="solid-equivalent height")
    ax.plot(t, [f(r, "precursor_required_volume_m3_m2") * 1e6 for r in pilot], label="required pore volume/area")
    ax.plot(t, [f(r, "precursor_filled_volume_m3_m2") * 1e6 for r in pilot], "--", label="filled pore volume/area")
    ax.set(xlabel="time (ks)", ylabel="equivalent length (µm)")
    ax.legend(frameon=False, fontsize=7); fig.tight_layout(); fig.savefig(FIG / "figure1_precursor_inventory.png"); plt.close(fig)

    fig, ax = plt.subplots(figsize=(3.5, 2.6)); ax2 = ax.twinx()
    ax.plot(t, [f(r, "top_soil_mass_fraction") for r in pilot], color="#0072B2", label="top NaCl mass fraction")
    ax2.plot(t, [f(r, "top_soil_sw") for r in pilot], color="#D55E00", label="top saturation")
    ax.set(xlabel="time (ks)", ylabel="NaCl mass fraction"); ax2.set_ylabel("water saturation")
    fig.tight_layout(); fig.savefig(FIG / "figure2_top_state.png"); plt.close(fig)

progressive = [r for r in convergence if math.isfinite(r["activation_s"])]
if progressive:
    fig, ax = plt.subplots(figsize=(3.5, 2.6))
    for dt in sorted({r["dt_s"] for r in progressive}, reverse=True):
        subset = sorted((r for r in progressive if r["dt_s"] == dt), key=lambda r: r["soil_cells"])
        ax.plot([r["soil_cells"] for r in subset], [r["activation_s"] / 1000 for r in subset], marker="o", label=f"dt={dt:g} s")
    ax.set(xlabel="soil vertical cells", ylabel="activation time (ks)"); ax.legend(frameon=False)
    fig.tight_layout(); fig.savefig(FIG / "figure3_activation_convergence.png"); plt.close(fig)

if activation:
    fig, ax = plt.subplots(figsize=(3.5, 2.6))
    labels = [r["mode"].replace("progressive_top_boundary_filling", "progressive").replace("adaptive_multilayer", "adaptive") for r in activation]
    ax.bar(labels, [r["initial_concentration_jump"] for r in activation], color=["#0072B2", "#E69F00", "#009E73"][:len(labels)])
    ax.axhline(0, color="black", lw=.8); ax.set_ylabel("ALE initial C − top-soil C")
    ax.tick_params(axis="x", rotation=20); fig.tight_layout(); fig.savefig(FIG / "figure4_activation_concentration_jump.png"); plt.close(fig)

pilot_dir = ROOT / "runs" / pilot_name
if (pilot_dir / "ale_diagnostics.csv").exists() and pilot:
    diag = [r for r in read_csv(pilot_dir / "ale_diagnostics.csv") if r["status"] == "accepted" and int(r["i"]) == 0]
    active_row = next((r for r in pilot if int(float(r["active"])) == 1), None)
    if diag and active_row:
        ta = f(active_row, "activation_s")
        fig, ax = plt.subplots(figsize=(3.5, 2.6)); ax2 = ax.twinx()
        ax.plot([(f(r, "time_s")-ta) for r in diag], [f(r, "epsilon") for r in diag], color="#0072B2", label="bottom porosity")
        post = [r for r in pilot if f(r, "time_s") >= ta]
        ax2.plot([(f(r, "time_s")-ta) for r in post], [f(r, "internal_dissolution_kg_m2")*1e3 for r in post], color="#D55E00", label="cumulative dissolution")
        ax.set(xlabel="time after activation (s)", ylabel="bottom-cell porosity"); ax2.set_ylabel("dissolution (g m⁻²)")
        fig.tight_layout(); fig.savefig(FIG / "figure5_postactivation_failure.png"); plt.close(fig)

if records:
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(3.5, 4.2), sharex=True)
    for record in records:
        if record["mode"] != "progressive_top_boundary_filling" or record["dt_s"] != 0.5 or record["ale_cells"] != 64:
            continue
        rows = histories[record["case"]]
        ax1.semilogy([f(r, "time_s")/1000 for r in rows], [max(abs(f(r, "salt_residual_kg")), 1e-20) for r in rows], label=f"N={record['soil_cells']}")
        ax2.semilogy([f(r, "time_s")/1000 for r in rows], [max(abs(f(r, "water_residual_kg")), 1e-20) for r in rows])
    ax1.set_ylabel("|salt residual| (kg)"); ax2.set(ylabel="|water residual| (kg)", xlabel="time (ks)")
    ax1.legend(frameon=False); fig.tight_layout(); fig.savefig(FIG / "figure6_conservation.png"); plt.close(fig)


def finite_max(key, subset):
    values = [r[key] for r in subset if math.isfinite(r[key])]
    return max(values) if values else math.inf


progressive_records = [r for r in records if r["mode"] == "progressive_top_boundary_filling"]
verification = {
    "slurm_job_id": os.environ["SLURM_JOB_ID"], "hostname": os.uname().nodename,
    "case_count": len(records), "progressive_case_count": len(progressive_records),
    "unit_test_passed": (ROOT / "unit" / "progressive_test.log").read_text().strip() == "progressive precursor tests passed",
    "all_progressive_activated": bool(progressive_records) and all(math.isfinite(r["activation_s"]) for r in progressive_records),
    "all_progressive_initial_volume_residual_le_1e_12_m": bool(progressive_records) and all(abs(r["precursor_unfilled_volume_m3_m2"]) <= 1e-12 for r in progressive_records if math.isfinite(r["activation_s"])),
    "max_progressive_salt_relative_residual": finite_max("max_salt_relative_residual", progressive_records),
    "max_progressive_water_relative_residual": finite_max("max_water_relative_residual", progressive_records),
    "global_mass_gate_le_1e_8": finite_max("max_salt_relative_residual", progressive_records) <= 1e-8 and finite_max("max_water_relative_residual", progressive_records) <= 1e-8,
    "all_requested_horizons_completed": bool(progressive_records) and all(r["completed_requested_horizon"] for r in progressive_records),
    "scientific_status": "INITIALIZATION_IMPROVED_BUT_INTERFACE_FAILURE_REMAINS" if any(math.isfinite(r["activation_s"]) and not r["completed_requested_horizon"] for r in progressive_records) else "PENDING_REVIEW",
}
(OUT / "verification.json").write_text(json.dumps(verification, indent=2) + "\n")

provenance = {"slurm_job_id": os.environ["SLURM_JOB_ID"], "hostname": os.uname().nodename, "files": {}}
for path in sorted([ROOT / "build" / "soil_solver", ROOT / "source.tar.gz", ROOT / "unit" / "source_binary.sha256"]):
    if path.exists():
        provenance["files"][str(path.relative_to(ROOT))] = hashlib.sha256(path.read_bytes()).hexdigest()
(OUT / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
print(json.dumps(verification, indent=2), flush=True)
