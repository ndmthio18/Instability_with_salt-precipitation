"""Compare censored x=0.08 reservoir pilots at 0.1 and 3.5 micrometres."""
import csv
import gzip
import hashlib
import json
import pathlib
import shutil

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image


root = pathlib.Path.cwd()
case = root / "pilot/pilot_ref_ale64"
diagnostic = root / "diagnostic_bound"
old = pathlib.Path(
    "/WORK/hulm_work/zhixin/p2_h_kr_av_reservoir_x008_20261108/analysis_x008_fail"
)
out = root / "analysis_h3p5"
results, figures, logs, source = (
    out / "results", out / "figures", out / "logs", out / "source"
)
for path in (results, figures, logs, source):
    path.mkdir(parents=True, exist_ok=True)


def read_rows(path):
    with path.open(newline="") as handle:
        return [{key: float(value) for key, value in row.items()} for row in csv.DictReader(handle)]


def compact(path, destination, interval=10.0):
    selected = []
    first = last = None
    next_time = 0.0
    maxima = {"salt": 0.0, "water": 0.0, "solid": 0.0}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames
        for raw in reader:
            row = {key: float(value) for key, value in raw.items()}
            first = first or row
            last = row
            maxima["salt"] = max(maxima["salt"], abs(row.get("salt_residual_kg", 0.0)))
            maxima["water"] = max(maxima["water"], abs(row.get("water_residual_kg", 0.0)))
            maxima["solid"] = max(maxima["solid"], abs(row.get("solid_ledger_residual_kg", 0.0)))
            if row["time_s"] + 1e-12 >= next_time:
                selected.append(row)
                while next_time <= row["time_s"] + 1e-12:
                    next_time += interval
    if selected[-1]["time_s"] != last["time_s"]:
        selected.append(last)
    with destination.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(selected)
    return first, last, selected, maxima


first, last, history, maxima = compact(
    case / "coupled_history.csv", results / "h3p5_history.csv"
)
_, boundary_last, boundary, _ = compact(
    case / "open_boundary_budget.csv", results / "h3p5_boundary.csv"
)
old_history = read_rows(old / "results/pilot_history.csv")
old_extrema = read_rows(old / "results/ale_extrema.csv")
old_summary = json.loads((old / "results/pilot_summary.json").read_text())

first_precipitation = activation = None
with (case / "coupled_history.csv").open(newline="") as handle:
    for raw in csv.DictReader(handle):
        row = {key: float(value) for key, value in raw.items()}
        if first_precipitation is None and row["top_solid_kg_m2"] > 0:
            first_precipitation = row["time_s"]
        if activation is None and row["active"] == 1:
            activation = row["activation_s"]
        if first_precipitation is not None and activation is not None:
            break

bound = json.loads((diagnostic / "bound_diagnostic.json").read_text())
last_accepted_time = bound["last_accepted"]["time_s"]
first_rejected_time = bound["first_rejected"]["time_s"]


def summarize_group(rows):
    eps = [float(row["epsilon"]) for row in rows]
    concentration = [float(row["c"]) for row in rows]
    min_index, max_index = eps.index(min(eps)), eps.index(max(eps))
    return {
        "time_s": float(rows[0]["time_s"]),
        "status": rows[0]["status"],
        "epsilon_min": eps[min_index],
        "epsilon_min_cell": int(float(rows[min_index]["i"])),
        "epsilon_min_z_m": 0.5 * (float(rows[min_index]["z_bottom_m"]) + float(rows[min_index]["z_top_m"])),
        "epsilon_max": eps[max_index],
        "epsilon_max_cell": int(float(rows[max_index]["i"])),
        "epsilon_max_z_m": 0.5 * (float(rows[max_index]["z_bottom_m"]) + float(rows[max_index]["z_top_m"])),
        "concentration_min": min(concentration),
        "concentration_max": max(concentration),
    }


extrema = []
accepted_profile = rejected_profile = None
with (diagnostic / "ale_diagnostics.csv").open(newline="") as handle:
    current_key, group = None, []
    for row in csv.DictReader(handle):
        key = (row["time_s"], row["status"])
        if group and key != current_key:
            summary = summarize_group(group)
            if summary["status"] == "accepted":
                extrema.append(summary)
            if summary["status"] == "accepted" and abs(summary["time_s"] - last_accepted_time) < 1e-9:
                accepted_profile = group
            if summary["status"] == "bound_rejected" and abs(summary["time_s"] - first_rejected_time) < 1e-9 and rejected_profile is None:
                rejected_profile = group
            group = []
        group.append(row)
        current_key = key
    if group:
        summary = summarize_group(group)
        if summary["status"] == "accepted":
            extrema.append(summary)
        if summary["status"] == "accepted" and abs(summary["time_s"] - last_accepted_time) < 1e-9:
            accepted_profile = group
        if summary["status"] == "bound_rejected" and abs(summary["time_s"] - first_rejected_time) < 1e-9 and rejected_profile is None:
            rejected_profile = group
if not accepted_profile or not rejected_profile:
    raise RuntimeError("missing h3p5 accepted/rejected profiles")

with (results / "h3p5_ale_extrema.csv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(extrema[0]))
    writer.writeheader()
    writer.writerows(extrema)
with (results / "h3p5_last_accepted_and_first_rejected.csv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(accepted_profile[0]))
    writer.writeheader()
    writer.writerows(accepted_profile)
    writer.writerows(rejected_profile)

attempts = list(csv.DictReader((case / "donor_attempts.csv").open(newline="")))
accepted_attempt = next(row for row in attempts if int(row["sufficient"]) == 1)
transfer_relative = max(
    abs(float(accepted_attempt["salt_transfer_residual_kg"])) / abs(float(accepted_attempt["initial_salt_kg"])),
    abs(float(accepted_attempt["water_transfer_residual_kg"])) / abs(float(accepted_attempt["initial_water_kg"])),
)
initial_salt, initial_water = abs(first["soil_salt_kg"]), abs(first["soil_water_kg"])
accepted_summary = summarize_group(accepted_profile)
rejected_summary = summarize_group(rejected_profile)
summary = {
    "status": "FAIL_ALE_UPPER_POROSITY_BOUND",
    "requested_end_s": 180000.0,
    "last_accepted_time_s": last["time_s"],
    "first_precipitation_time_s": first_precipitation,
    "activation_time_s": activation,
    "activation_height_m": float(accepted_attempt["h_m"]),
    "accepted_ale_window_s": last["time_s"] - activation,
    "donor_count": int(accepted_attempt["donor_count"]),
    "used_donor_depth_m": float(accepted_attempt["used_depth_m"]),
    "required_pore_liquid_kg": float(accepted_attempt["required_kg"]),
    "last_height_m": last["h_m"],
    "last_mass_height_m": last["h_mass_m"],
    "last_top_solid_kg_m2": last["top_solid_kg_m2"],
    "last_internal_solid_kg_m2": last["internal_solid_kg_m2"],
    "last_internal_positive_kg_m2": last["internal_positive_kg_m2"],
    "last_internal_dissolution_kg_m2": last["internal_dissolution_kg_m2"],
    "last_epsilon_min": accepted_summary["epsilon_min"],
    "last_epsilon_min_cell": accepted_summary["epsilon_min_cell"],
    "last_epsilon_max": accepted_summary["epsilon_max"],
    "last_epsilon_max_cell": accepted_summary["epsilon_max_cell"],
    "first_rejected_time_s": first_rejected_time,
    "first_rejected_epsilon_max": rejected_summary["epsilon_max"],
    "last_concentration_min": accepted_summary["concentration_min"],
    "last_concentration_max": accepted_summary["concentration_max"],
    "cumulative_evaporation_kg_m2": last["vapour_kg"] / 0.001,
    "cumulative_bottom_water_kg_m2": boundary_last["cumulative_bottom_water_kg"] / 0.001,
    "cumulative_bottom_salt_kg_m2": boundary_last["cumulative_bottom_salt_kg"] / 0.001,
    "salt_relative_residual": maxima["salt"] / initial_salt,
    "water_relative_residual": maxima["water"] / initial_water,
    "solid_ledger_abs_residual_kg": maxima["solid"],
    "activation_transfer_relative_residual": transfer_relative,
}
(results / "h3p5_summary.json").write_text(json.dumps(summary, indent=2) + "\n")
comparison = {
    "activation_delay_s": summary["activation_time_s"] - old_summary["activation_time_s"],
    "accepted_window_ratio_h3p5_over_h0p1": summary["accepted_ale_window_s"] / old_summary["accepted_ale_window_s"],
    "failure_height_ratio_h3p5_over_h0p1": summary["last_height_m"] / old_summary["last_height_m"],
    "failure_mass_height_ratio_h3p5_over_h0p1": summary["last_mass_height_m"] / old_summary["last_mass_height_m"],
    "same_failure_mechanism": "ALE soil-side cell reaches epsilon=0.95 after net internal dissolution",
}
(results / "comparison.json").write_text(json.dumps(comparison, indent=2) + "\n")

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8.5, "axes.labelsize": 9,
    "legend.fontsize": 8, "lines.linewidth": 1.5, "axes.linewidth": 0.8,
})
blue, orange, green, purple = "#0072B2", "#D55E00", "#009E73", "#CC79A7"


def save(fig, name):
    for extension in ("png", "pdf", "svg"):
        fig.savefig(figures / f"{name}.{extension}", dpi=300 if extension == "png" else None, bbox_inches="tight")
    plt.close(fig)


fig, axes = plt.subplots(2, 1, figsize=(7.2, 5.8), sharex=True, constrained_layout=True)
for data, label, color, style, act in (
    (old_history, "0.1 µm activation", blue, "-", old_summary["activation_time_s"]),
    (history, "3.5 µm activation", orange, "--", activation),
):
    time = [row["time_s"] / 3600 for row in data]
    axes[0].plot(time, [row["top_soil_mass_fraction"] for row in data], color=color, ls=style, label=label)
    axes[1].plot(time, [row["h_mass_m"] * 1e6 for row in data], color=color, ls=style, label=label)
    axes[0].axvline(act / 3600, color=color, ls=style, alpha=0.7)
axes[0].axhline(0.264940870223385, color="0.25", ls=":", label="equilibrium")
axes[0].set_ylabel("Surface NaCl mass fraction")
axes[1].set_ylabel("Solid-equivalent height (µm)")
axes[1].set_xlabel("Physical time (h)")
axes[0].legend(frameon=False, ncol=2)
save(fig, "figure1_activation_comparison")

fig, axes = plt.subplots(2, 2, figsize=(8.5, 6.2), sharex="col", constrained_layout=True)
for data, ext, label, color, style, act in (
    (old_history, old_extrema, "0.1 µm", blue, "-", old_summary["activation_time_s"]),
    (history, extrema, "3.5 µm", orange, "--", activation),
):
    active = [row for row in data if row["active"] == 1]
    elapsed = [(row["time_s"] - act) for row in active]
    axes[0, 0].plot(elapsed, [row["h_m"] * 1e6 for row in active], color=color, ls=style, label=label)
    axes[1, 0].plot(elapsed, [row["internal_solid_kg_m2"] * 1e3 for row in active], color=color, ls=style, label=label)
    ext_active = [row for row in ext if row["time_s"] >= act]
    axes[0, 1].plot([row["time_s"] - act for row in ext_active], [row["epsilon_max"] for row in ext_active], color=color, ls=style, label=label)
    axes[1, 1].plot([row["time_s"] - act for row in ext_active], [row["epsilon_min"] for row in ext_active], color=color, ls=style, label=label)
axes[0, 1].axhline(0.95, color="0.25", ls=":")
axes[1, 1].axhline(0.02, color="0.25", ls=":")
axes[0, 0].set_ylabel("Actual ALE height (µm)")
axes[1, 0].set_ylabel(r"Net internal solid (g m$^{-2}$)")
axes[0, 1].set_ylabel(r"$\epsilon_{max}$")
axes[1, 1].set_ylabel(r"$\epsilon_{min}$")
for axis in axes[1]:
    axis.set_xlabel("Time after activation (s)")
for axis in axes.flat:
    axis.legend(frameon=False)
save(fig, "figure2_post_activation_comparison")

old_profiles = list(csv.DictReader((old / "results/ale_last_accepted_and_first_rejected.csv").open(newline="")))
fig, axes = plt.subplots(2, 2, figsize=(8.5, 6.3), constrained_layout=True)
for column, (profiles, title) in enumerate(((old_profiles, "0.1 µm"), (accepted_profile + rejected_profile, "3.5 µm"))):
    statuses = []
    for status, label, color, style in (("accepted", "last accepted", blue, "-"), ("bound_rejected", "first rejected", orange, "--")):
        rows = [row for row in profiles if row["status"] == status]
        if not rows:
            continue
        height = max(float(row["z_top_m"]) for row in rows)
        z = [0.5 * (float(row["z_bottom_m"]) + float(row["z_top_m"])) / height for row in rows]
        axes[0, column].plot([float(row["epsilon"]) for row in rows], z, color=color, ls=style, label=label)
        axes[1, column].plot([(float(row["c"]) - 0.264940870223385) * 1e6 for row in rows], z, color=color, ls=style, label=label)
    axes[0, column].axvline(0.95, color="0.25", ls=":")
    axes[1, column].axvline(0, color="0.25", ls=":")
    axes[0, column].set_title(title)
    axes[0, column].set_xlabel("Porosity")
    axes[1, column].set_xlabel(r"$10^6(C-C_{eq})$")
    axes[0, column].legend(frameon=False)
for axis in axes[:, 0]:
    axis.set_ylabel("Relative height z/h")
save(fig, "figure3_failure_profiles")

fig, axes = plt.subplots(1, 2, figsize=(8.2, 3.7), constrained_layout=True)
labels = ["0.1 µm", "3.5 µm"]
axes[0].bar(labels, [old_summary["accepted_ale_window_s"], summary["accepted_ale_window_s"]], color=[blue, orange])
axes[0].set_ylabel("Accepted post-activation window (s)")
axes[1].bar(labels, [old_summary["last_height_m"] * 1e6, summary["last_height_m"] * 1e6], color=[blue, orange])
axes[1].set_ylabel("Actual height at failure (µm)")
for axis in axes:
    axis.set_ylim(bottom=0)
save(fig, "figure4_censored_endpoints")

for filename in ("case_spec.json", "command.json", "outcome.json", "params.input", "reservoir_config.csv", "coupling_events.csv", "donor_attempts.csv"):
    shutil.copy2(case / filename, logs / filename)
for filename in ("bound_diagnostic.json", "command.json", "outcome.json"):
    shutil.copy2(diagnostic / filename, logs / f"diagnostic_{filename}")
shutil.copy2(root / "pilot_gate.json", logs / "pilot_gate.json")
for original, destination in ((case / "solver.out", logs / "pilot_solver.out.gz"), (diagnostic / "solver.out", logs / "diagnostic_solver.out.gz"), (diagnostic / "ale_diagnostics.csv", logs / "diagnostic_ale_diagnostics.csv.gz")):
    with original.open("rb") as input_handle, gzip.open(destination, "wb", compresslevel=9) as output_handle:
        shutil.copyfileobj(input_handle, output_handle)
for filename in ("run_reservoir.py", "verify_reservoir.py", "diagnose_bound_from_pilot.py", "postprocess_h3p5_compare.py", "README.md", "RESERVOIR_FORMULATION.md"):
    shutil.copy2(root / "repo/model/p2_h_ale_coupling" / filename, source / filename)

checks = []
for stem in ("figure1_activation_comparison", "figure2_post_activation_comparison", "figure3_failure_profiles", "figure4_censored_endpoints"):
    image = Image.open(figures / f"{stem}.png").convert("RGB")
    checks.append({"stem": stem, "pixels": image.size, "nonblank": any(lo < hi for lo, hi in image.getextrema()), "formats": all((figures / f"{stem}.{ext}").stat().st_size > 1000 for ext in ("png", "pdf", "svg"))})
figure_checks = {"passed": all(row["nonblank"] and row["formats"] for row in checks), "checks": checks}
(out / "figure_checks.json").write_text(json.dumps(figure_checks, indent=2) + "\n")
verification = {
    "registered_activation_height_passed": summary["activation_height_m"] >= 3.5e-6,
    "activation_transfer_passed": transfer_relative <= 1e-12,
    "numerical_conservation_passed": summary["salt_relative_residual"] <= 1e-8 and summary["water_relative_residual"] <= 1e-8,
    "requested_horizon_passed": False,
    "same_upper_porosity_failure_reproduced": rejected_summary["epsilon_max"] >= 0.95 and rejected_summary["epsilon_max_cell"] == 0,
    "scientific_status": "FAIL_ALE_UPPER_POROSITY_BOUND",
    "figures_passed": figure_checks["passed"],
    "overall_passed": False,
}
(out / "verification.json").write_text(json.dumps(verification, indent=2) + "\n")
provenance = {
    "branch": "research/p2-h-reservoir-x008-h3p5",
    "commit": (root / "unit/COMMIT_SHA.txt").read_text().strip(),
    "hpcc_root": str(root),
    "solver_sha256": hashlib.sha256((root / "build/soil_solver").read_bytes()).hexdigest(),
    "jobs": {"pilot": "16438585", "bound_diagnostic": "16438590"},
    "comparison_source": str(old),
}
(out / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
with (out / "CHECKSUMS.sha256").open("w") as manifest:
    for path in sorted(item for item in out.rglob("*") if item.is_file() and item.name != "CHECKSUMS.sha256"):
        manifest.write(hashlib.sha256(path.read_bytes()).hexdigest() + "  " + str(path.relative_to(out)) + "\n")
print(json.dumps({"summary": summary, "comparison": comparison, "verification": verification}, indent=2))
raise SystemExit(0 if figure_checks["passed"] else 1)
