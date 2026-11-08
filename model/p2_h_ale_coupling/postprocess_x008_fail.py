"""Package and plot the censored x=0.08 reservoir pilot on a compute node."""
import csv
import gzip
import hashlib
import json
import pathlib
import shutil

import matplotlib
matplotlib.use("Agg")
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from PIL import Image


root = pathlib.Path.cwd()
case = root / "pilot/pilot_ref_ale64"
diagnostic = root / "diagnostic_bound"
out = root / "analysis_x008_fail"
results = out / "results"
figures = out / "figures"
logs = out / "logs"
source = out / "source"
for path in (results, figures, logs, source):
    path.mkdir(parents=True, exist_ok=True)


def floats(row):
    return {key: float(value) for key, value in row.items()}


def compact(path, destination, interval=10.0):
    selected = []
    first = last = None
    next_time = 0.0
    maxima = {"salt": 0.0, "water": 0.0, "solid": 0.0}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames
        for raw in reader:
            row = floats(raw)
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


first, last, history, maxima = compact(case / "coupled_history.csv", results / "pilot_history.csv")
_, boundary_last, boundary, _ = compact(case / "open_boundary_budget.csv", results / "pilot_boundary.csv")

with (case / "coupled_history.csv").open(newline="") as handle:
    first_precipitation = first_activation = None
    for raw in csv.DictReader(handle):
        row = floats(raw)
        if first_precipitation is None and row["top_solid_kg_m2"] > 0:
            first_precipitation = row["time_s"]
        if first_activation is None and row["active"] == 1:
            first_activation = row["activation_s"]
        if first_precipitation is not None and first_activation is not None:
            break

bound = json.loads((diagnostic / "bound_diagnostic.json").read_text())
last_accepted_time = bound["last_accepted"]["time_s"]
first_rejected_time = bound["first_rejected"]["time_s"]


def summarize_group(rows):
    eps = [float(row["epsilon"]) for row in rows]
    concentration = [float(row["c"]) for row in rows]
    min_eps_index = eps.index(min(eps))
    max_eps_index = eps.index(max(eps))
    min_c_index = concentration.index(min(concentration))
    max_c_index = concentration.index(max(concentration))
    return {
        "time_s": float(rows[0]["time_s"]),
        "status": rows[0]["status"],
        "epsilon_min": eps[min_eps_index],
        "epsilon_min_cell": int(float(rows[min_eps_index]["i"])),
        "epsilon_min_z_m": 0.5 * (float(rows[min_eps_index]["z_bottom_m"]) + float(rows[min_eps_index]["z_top_m"])),
        "epsilon_max": eps[max_eps_index],
        "epsilon_max_cell": int(float(rows[max_eps_index]["i"])),
        "epsilon_max_z_m": 0.5 * (float(rows[max_eps_index]["z_bottom_m"]) + float(rows[max_eps_index]["z_top_m"])),
        "concentration_min": concentration[min_c_index],
        "concentration_min_cell": int(float(rows[min_c_index]["i"])),
        "concentration_max": concentration[max_c_index],
        "concentration_max_cell": int(float(rows[max_c_index]["i"])),
    }


extrema = []
accepted_profile = rejected_profile = None
with (diagnostic / "ale_diagnostics.csv").open(newline="") as handle:
    current_key = None
    group = []
    for row in csv.DictReader(handle):
        key = (row["time_s"], row["status"])
        if group and key != current_key:
            summary = summarize_group(group)
            if summary["status"] == "accepted":
                extrema.append(summary)
            if abs(summary["time_s"] - last_accepted_time) < 1e-9 and summary["status"] == "accepted":
                accepted_profile = group
            if abs(summary["time_s"] - first_rejected_time) < 1e-9 and summary["status"] == "bound_rejected" and rejected_profile is None:
                rejected_profile = group
            group = []
        group.append(row)
        current_key = key
    if group:
        summary = summarize_group(group)
        if summary["status"] == "accepted":
            extrema.append(summary)
        if abs(summary["time_s"] - last_accepted_time) < 1e-9 and summary["status"] == "accepted":
            accepted_profile = group
        if abs(summary["time_s"] - first_rejected_time) < 1e-9 and summary["status"] == "bound_rejected" and rejected_profile is None:
            rejected_profile = group

if not accepted_profile or not rejected_profile:
    raise RuntimeError("missing accepted/rejected ALE diagnostic profiles")

with (results / "ale_extrema.csv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=list(extrema[0]))
    writer.writeheader()
    writer.writerows(extrema)

profile_fieldnames = list(accepted_profile[0])
with (results / "ale_last_accepted_and_first_rejected.csv").open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=profile_fieldnames)
    writer.writeheader()
    writer.writerows(accepted_profile)
    writer.writerows(rejected_profile)

attempts = list(csv.DictReader((case / "donor_attempts.csv").open(newline="")))
accepted_attempts = [row for row in attempts if int(row["sufficient"]) == 1]
transfer_relative = 0.0
for row in accepted_attempts:
    transfer_relative = max(
        transfer_relative,
        abs(float(row["salt_transfer_residual_kg"])) / max(abs(float(row["initial_salt_kg"])), 1e-30),
        abs(float(row["water_transfer_residual_kg"])) / max(abs(float(row["initial_water_kg"])), 1e-30),
    )

initial_salt = abs(first["soil_salt_kg"])
initial_water = abs(first["soil_water_kg"])
accepted_summary = summarize_group(accepted_profile)
rejected_summary = summarize_group(rejected_profile)
summary = {
    "status": "FAIL_ALE_UPPER_POROSITY_BOUND",
    "requested_end_s": 180000.0,
    "last_accepted_time_s": last["time_s"],
    "first_precipitation_time_s": first_precipitation,
    "activation_time_s": first_activation,
    "accepted_ale_window_s": last["time_s"] - first_activation,
    "initial_surface_mass_fraction": first["top_soil_mass_fraction"],
    "last_surface_mass_fraction": last["top_soil_mass_fraction"],
    "last_surface_sw": last["top_soil_sw"],
    "last_height_m": last["h_m"],
    "last_mass_height_m": last["h_mass_m"],
    "last_top_solid_kg_m2": last["top_solid_kg_m2"],
    "last_internal_solid_kg_m2": last["internal_solid_kg_m2"],
    "last_internal_positive_kg_m2": last["internal_positive_kg_m2"],
    "last_internal_dissolution_kg_m2": last["internal_dissolution_kg_m2"],
    "last_internal_fraction_net": last["internal_fraction_net"],
    "last_epsilon_min": accepted_summary["epsilon_min"],
    "last_epsilon_max": accepted_summary["epsilon_max"],
    "last_epsilon_max_cell": accepted_summary["epsilon_max_cell"],
    "last_epsilon_max_z_m": accepted_summary["epsilon_max_z_m"],
    "first_rejected_time_s": first_rejected_time,
    "first_rejected_epsilon_max": rejected_summary["epsilon_max"],
    "first_rejected_epsilon_max_cell": rejected_summary["epsilon_max_cell"],
    "first_rejected_epsilon_max_z_m": rejected_summary["epsilon_max_z_m"],
    "last_concentration_min": accepted_summary["concentration_min"],
    "last_concentration_max": accepted_summary["concentration_max"],
    "cumulative_evaporation_kg_m2": last["vapour_kg"] / 0.001,
    "cumulative_bottom_water_kg_m2": boundary_last["cumulative_bottom_water_kg"] / 0.001,
    "cumulative_bottom_salt_kg_m2": boundary_last["cumulative_bottom_salt_kg"] / 0.001,
    "salt_relative_residual": maxima["salt"] / initial_salt,
    "water_relative_residual": maxima["water"] / initial_water,
    "solid_ledger_abs_residual_kg": maxima["solid"],
    "activation_transfer_relative_residual": transfer_relative,
    "accepted_activation_attempts": len(accepted_attempts),
    "matrix_status": "NOT_RUN_PILOT_GATE_FAILED",
}
(results / "pilot_summary.json").write_text(json.dumps(summary, indent=2) + "\n")

with gzip.open(results / "soil_fields.csv.gz", "wb", compresslevel=9) as output, (case / "soil_fields.csv").open("rb") as source_handle:
    shutil.copyfileobj(source_handle, output)
shutil.copy2(case / "crust_fields.csv", results / "crust_fields.csv")

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8.5, "axes.labelsize": 9,
    "legend.fontsize": 8, "lines.linewidth": 1.5, "axes.linewidth": 0.8,
})
blue, orange, green, purple = "#0072B2", "#D55E00", "#009E73", "#CC79A7"


def save(fig, name):
    for extension in ("png", "pdf", "svg"):
        fig.savefig(figures / f"{name}.{extension}", dpi=300 if extension == "png" else None, bbox_inches="tight")
    plt.close(fig)


time_h = [row["time_s"] / 3600 for row in history]
activation_h = first_activation / 3600
stop_h = last["time_s"] / 3600

fig, axes = plt.subplots(2, 1, figsize=(7.2, 5.8), sharex=True, constrained_layout=True)
axes[0].plot(time_h, [row["top_soil_mass_fraction"] for row in history], color=blue, label=r"$x_{NaCl}=0.08$")
axes[0].axhline(0.264940870223385, color="0.25", ls=":", label="equilibrium")
axes[1].plot(time_h, [row["top_soil_sw"] for row in history], color=green)
for axis in axes:
    axis.axvline(activation_h, color=orange, ls="--", label="ALE activation" if axis is axes[0] else None)
    axis.axvline(stop_h, color=purple, ls=":", label="upper-porosity stop" if axis is axes[0] else None)
axes[0].set_ylabel("Surface NaCl mass fraction")
axes[1].set_ylabel(r"Surface $S_w$")
axes[1].set_xlabel("Time (h)")
axes[1].ticklabel_format(axis="y", style="plain", useOffset=False)
axes[0].legend(frameon=False, ncol=2)
save(fig, "figure1_surface_activation")

fig, axes = plt.subplots(2, 2, figsize=(8.5, 6.2), sharex=True, constrained_layout=True)
axes[0, 0].plot(time_h, [row["h_m"] * 1e6 for row in history], color=blue)
axes[0, 0].set_ylabel("Actual ALE height (µm)")
ext_t = [row["time_s"] / 3600 for row in extrema]
axes[0, 1].plot(ext_t, [row["epsilon_min"] for row in extrema], color=blue, label=r"$\epsilon_{min}$")
axes[0, 1].plot(ext_t, [row["epsilon_max"] for row in extrema], color=orange, ls="--", label=r"$\epsilon_{max}$")
axes[0, 1].axhline(0.95, color="0.25", ls=":")
axes[0, 1].set_ylabel("ALE porosity extrema")
axes[0, 1].legend(frameon=False)
axes[1, 0].plot(time_h, [row["top_solid_kg_m2"] for row in history], color=blue, label="top solid")
axes[1, 0].plot(time_h, [row["internal_solid_kg_m2"] for row in history], color=orange, ls="--", label="net internal solid")
axes[1, 0].set_ylabel(r"Solid inventory (kg m$^{-2}$)")
axes[1, 0].legend(frameon=False)
axes[1, 1].plot(time_h, [row["internal_positive_kg_m2"] for row in history], color=green, label="internal precipitation")
axes[1, 1].plot(time_h, [row["internal_dissolution_kg_m2"] for row in history], color=purple, ls="--", label="internal dissolution")
axes[1, 1].set_ylabel(r"Cumulative magnitude (kg m$^{-2}$)")
axes[1, 1].legend(frameon=False)
for axis in axes[1]:
    axis.set_xlabel("Time (h)")
for axis in axes.flat:
    axis.axvline(activation_h, color="0.4", lw=0.8, ls=":")
save(fig, "figure2_ale_evolution")

fig, axes = plt.subplots(1, 2, figsize=(8.4, 3.8), constrained_layout=True)
for rows, label, color, style in (
    (accepted_profile, "last accepted", blue, "-"),
    (rejected_profile, "first rejected trial", orange, "--"),
):
    z = [0.5 * (float(row["z_bottom_m"]) + float(row["z_top_m"])) * 1e6 for row in rows]
    axes[0].plot([float(row["epsilon"]) for row in rows], z, color=color, ls=style, label=label)
    axes[1].plot([(float(row["c"]) - 0.264940870223385) * 1e6 for row in rows], z, color=color, ls=style, label=label)
axes[0].axvline(0.95, color="0.25", ls=":", label="upper bound")
axes[1].axvline(0, color="0.25", ls=":")
axes[0].set_xlabel("Porosity")
axes[1].set_xlabel(r"$10^6(C-C_{eq})$")
for axis in axes:
    axis.set_ylabel("Height above soil (µm)")
axes[0].legend(frameon=False)
save(fig, "figure3_last_accepted_rejected_profiles")

boundary_t = [row["time_s"] / 3600 for row in boundary]
fig, axes = plt.subplots(2, 2, figsize=(8.5, 6.2), constrained_layout=True)
axes[0, 0].plot(boundary_t, [row["bottom_water_in_kg_m2_s"] for row in boundary], color=blue, label="bottom water influx")
axes[0, 0].plot(time_h, [row["evaporation_kg_m2_s"] for row in history], color=orange, ls="--", label="top evaporation")
axes[0, 0].set_ylabel(r"Flux (kg m$^{-2}$ s$^{-1}$)")
axes[0, 0].legend(frameon=False)
axes[0, 1].plot(boundary_t, [row["cumulative_bottom_water_kg"] / 0.001 for row in boundary], color=blue, label="bottom water")
axes[0, 1].plot(time_h, [row["vapour_kg"] / 0.001 for row in history], color=orange, ls="--", label="evaporated water")
axes[0, 1].set_ylabel(r"Cumulative water (kg m$^{-2}$)")
axes[0, 1].legend(frameon=False)
axes[1, 0].plot(time_h, [row["water_residual_kg"] / initial_water for row in history], color=blue)
axes[1, 0].set_ylabel("Relative water residual")
axes[1, 1].plot(time_h, [row["salt_residual_kg"] / initial_salt for row in history], color=orange)
axes[1, 1].set_ylabel("Relative salt residual")
for axis in axes[1]:
    axis.set_xlabel("Time (h)")
save(fig, "figure4_flux_and_conservation")

crust_groups = {}
with (case / "crust_fields.csv").open(newline="") as handle:
    for raw in csv.DictReader(handle):
        row = floats(raw)
        crust_groups.setdefault(row["time_s"], []).append(row)
target_times = [first_activation + offset for offset in (100, 400, 800, 1100)]
selected_times = []
for target in target_times:
    chosen = min(crust_groups, key=lambda value: abs(value - target))
    if chosen not in selected_times:
        selected_times.append(chosen)
norm = mpl.colors.Normalize(vmin=min(row["epsilon"] for rows in crust_groups.values() for row in rows), vmax=0.95)
cmap = mpl.colormaps["viridis"]
fig, axis = plt.subplots(figsize=(7.2, 4.2), constrained_layout=True)
for column, time in enumerate(selected_times):
    for row in crust_groups[time]:
        bottom = row["z_bottom_m"] * 1e6
        height = (row["z_top_m"] - row["z_bottom_m"]) * 1e6
        axis.add_patch(Rectangle((column - 0.34, bottom), 0.68, height, facecolor=cmap(norm(row["epsilon"])), edgecolor="none"))
axis.set_xlim(-0.7, len(selected_times) - 0.3)
axis.set_ylim(0, max(row["z_top_m"] for time in selected_times for row in crust_groups[time]) * 1e6 * 1.05)
axis.set_xticks(range(len(selected_times)), [f"+{(time-first_activation):.0f} s" for time in selected_times])
axis.set_xlabel("Time after ALE activation")
axis.set_ylabel("Actual height above soil (µm)")
fig.colorbar(mpl.cm.ScalarMappable(norm=norm, cmap=cmap), ax=axis, label="Local ALE porosity")
save(fig, "figure5_porosity_columns")

for name in (
    "command.json", "params.input", "outcome.json", "reservoir_config.csv",
    "boundary_faces.csv", "coupling_events.csv", "donor_attempts.csv",
):
    shutil.copy2(case / name, logs / name)
for name in ("bound_diagnostic.json", "command.json", "outcome.json"):
    shutil.copy2(diagnostic / name, logs / f"diagnostic_{name}")
for name in ("verification_checks.json",):
    shutil.copy2(root / name, logs / name)
for name in ("solver.out",):
    with (case / name).open("rb") as input_handle, gzip.open(logs / f"pilot_{name}.gz", "wb", compresslevel=9) as output_handle:
        shutil.copyfileobj(input_handle, output_handle)
    with (diagnostic / name).open("rb") as input_handle, gzip.open(logs / f"diagnostic_{name}.gz", "wb", compresslevel=9) as output_handle:
        shutil.copyfileobj(input_handle, output_handle)
with (diagnostic / "ale_diagnostics.csv").open("rb") as input_handle, gzip.open(logs / "diagnostic_ale_diagnostics.csv.gz", "wb", compresslevel=9) as output_handle:
    shutil.copyfileobj(input_handle, output_handle)
for name in (
    "README.md", "soil_methods.inc", "RESERVOIR_FORMULATION.md", "run_reservoir.py",
    "verify_reservoir.py", "diagnose_x008_bound.py", "postprocess_x008_fail.py",
):
    shutil.copy2(root / "repo/model/p2_h_ale_coupling" / name, source / name)

checks = []
for stem in (
    "figure1_surface_activation", "figure2_ale_evolution",
    "figure3_last_accepted_rejected_profiles", "figure4_flux_and_conservation",
    "figure5_porosity_columns",
):
    image = Image.open(figures / f"{stem}.png").convert("RGB")
    checks.append({
        "stem": stem,
        "pixels": image.size,
        "nonblank": any(low < high for low, high in image.getextrema()),
        "formats": all((figures / f"{stem}.{extension}").stat().st_size > 1000 for extension in ("png", "pdf", "svg")),
    })
figure_checks = {"passed": all(row["nonblank"] and row["formats"] for row in checks), "checks": checks}
(out / "figure_checks.json").write_text(json.dumps(figure_checks, indent=2) + "\n")

verification = {
    "numerical_conservation_passed": summary["salt_relative_residual"] <= 1e-8 and summary["water_relative_residual"] <= 1e-8,
    "activation_transfer_passed": transfer_relative <= 1e-12 and bool(accepted_attempts),
    "ale_activation_passed": first_activation is not None and summary["accepted_ale_window_s"] >= 100,
    "requested_horizon_passed": last["time_s"] >= 180000.0 - 1e-8,
    "scientific_status": summary["status"],
    "matrix_status": summary["matrix_status"],
    "figures_passed": figure_checks["passed"],
    "overall_passed": False,
}
(out / "verification.json").write_text(json.dumps(verification, indent=2) + "\n")
provenance = {
    "branch": "research/p2-h-kr-av-reservoir-x008",
    "commit": (root / "unit/COMMIT_SHA.txt").read_text().strip(),
    "hpcc_root": str(root),
    "solver_sha256": hashlib.sha256((root / "build/soil_solver").read_bytes()).hexdigest(),
    "jobs": {"verification": "16437446", "pilot": "16437448", "bound_diagnostic": "16438191"},
}
(out / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")
with (out / "CHECKSUMS.sha256").open("w") as manifest:
    for path in sorted(item for item in out.rglob("*") if item.is_file() and item.name != "CHECKSUMS.sha256"):
        manifest.write(hashlib.sha256(path.read_bytes()).hexdigest() + "  " + str(path.relative_to(out)) + "\n")
print(json.dumps({"summary": summary, "verification": verification}, indent=2))
raise SystemExit(0 if figure_checks["passed"] else 1)
