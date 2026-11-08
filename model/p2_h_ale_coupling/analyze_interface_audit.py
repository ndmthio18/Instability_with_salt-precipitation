"""Postprocess the registered interface audit without altering raw outputs."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import os
import pathlib
from collections import defaultdict

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


ROOT = pathlib.Path.cwd()
OUT = ROOT / "analysis_output"
FIG = OUT / "figures"
OUT.mkdir(exist_ok=True)
FIG.mkdir(exist_ok=True)
CEQ = 0.2649408697676898
RHO = 1202.1685774416314
DIFFUSION = 1.3e-9
AREA = 0.001
TARGETS = [0.5, 1.0, 5.0, 10.0, 30.0, 100.0, 300.0, 500.0]


def convert(value: str):
    try:
        return float(value)
    except (TypeError, ValueError):
        return value


def read_csv(path: pathlib.Path):
    with path.open(newline="") as handle:
        return [{key: convert(value) for key, value in row.items()}
                for row in csv.DictReader(handle)]


def write_csv(path: pathlib.Path, rows, fields=None):
    rows = list(rows)
    if fields is None:
        fields = list(rows[0]) if rows else []
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def closest(rows, target, tolerance=0.251):
    if not rows:
        return None
    row = min(rows, key=lambda r: abs(r["time_s"] - target))
    return row if abs(row["time_s"] - target) <= tolerance else None


def successful_donor(case_dir):
    rows = read_csv(case_dir / "donor_attempts.csv")
    rows = [r for r in rows if int(r["sufficient"]) == 1]
    if not rows:
        return None
    r = rows[0]
    r["c_initial"] = r["salt_transfer_kg"] / (r["salt_transfer_kg"] + r["water_transfer_kg"])
    return r


def diagnostic_groups(case_dir):
    groups = defaultdict(list)
    path = case_dir / "ale_diagnostics.csv"
    if not path.exists():
        return groups
    for row in read_csv(path):
        if row["status"] == "accepted":
            groups[row["time_s"]].append(row)
    for values in groups.values():
        values.sort(key=lambda r: int(r["i"]))
    return groups


case_dirs = sorted((ROOT / "coupled").glob("*"))
cases = {}
summaries = []
for directory in case_dirs:
    if not (directory / "case_spec.json").exists():
        continue
    spec = json.loads((directory / "case_spec.json").read_text())
    history = read_csv(directory / "coupled_history.csv")
    boundary = read_csv(directory / "open_boundary_budget.csv")
    diagnostics = diagnostic_groups(directory)
    donor = successful_donor(directory)
    active = [r for r in history if int(r["active"]) == 1]
    activation = active[0]["activation_s"] if active else math.nan
    final = history[-1]
    outcome = json.loads((directory / "outcome.json").read_text())
    last_diag = diagnostics[max(diagnostics)] if diagnostics else []
    eps_max = max((r["epsilon"] for r in last_diag), default=math.nan)
    c_bottom = last_diag[0]["bottom_c"] if last_diag else math.nan
    first_precip = next((r["time_s"] for r in history if r["top_solid_kg_m2"] > 0), math.nan)
    status = "COMPLETED" if outcome["exit_code"] == 0 else "STOPPED"
    if not active:
        status = "NO_ALE_ACTIVATION"
    summary = {
        "case": spec["name"], "status": status, "exit_code": outcome["exit_code"],
        "x0_mole_fraction": spec["salt_mole_fraction"], "bottom_mode": spec["bottom_mode"],
        "soil_cells": spec["soil_cells"], "ale_cells": spec["ale_cells"],
        "activation_threshold_m": spec["activation_height_m"],
        "first_precipitation_s": first_precip, "activation_s": activation,
        "initial_ale_mass_fraction": donor["c_initial"] if donor else math.nan,
        "initial_ale_minus_ceq": donor["c_initial"]-CEQ if donor else math.nan,
        "donor_count": donor["donor_count"] if donor else 0,
        "donor_depth_m": donor["used_depth_m"] if donor else math.nan,
        "accepted_end_s": final["time_s"], "accepted_post_activation_s": final["time_s"]-activation if active else math.nan,
        "final_height_m": final["h_m"], "final_internal_positive_kg_m2": final["internal_positive_kg_m2"],
        "final_internal_dissolution_kg_m2": final["internal_dissolution_kg_m2"],
        "final_internal_net_kg_m2": final["internal_solid_kg_m2"],
        "final_epsilon_min": final["epsilon_min"], "final_epsilon_max": eps_max,
        "final_interface_mass_fraction": c_bottom,
        "salt_relative_residual": abs(final["salt_residual_kg"]) / max(abs(history[0]["soil_salt_kg"]), 1e-30),
        "water_relative_residual": abs(final["water_residual_kg"]) / max(abs(history[0]["soil_water_kg"]), 1e-30),
        "slurm_job_id": outcome["slurm_job_id"], "hostname": outcome["hostname"],
    }
    summaries.append(summary)
    cases[spec["name"]] = {
        "dir": directory, "spec": spec, "history": history, "boundary": boundary,
        "diagnostics": diagnostics, "donor": donor, "summary": summary,
    }
write_csv(OUT / "CASE_COMPARISON.csv", summaries)


def concentration_fields(prefix, row, ceq=CEQ):
    value = row if row is not None else math.nan
    return {
        prefix: value,
        prefix + "_minus_ceq": value-ceq if isinstance(value, float) and math.isfinite(value) else math.nan,
        prefix + "_relative": (value-ceq)/ceq if isinstance(value, float) and math.isfinite(value) else math.nan,
    }


jump_rows = []
flux_rows = []
audit_names = ["x008_reservoir_h0p1_n80_a64", "x008_reservoir_n80_a64"]
for name in audit_names:
    case = cases[name]
    history = case["history"]
    activation = case["summary"]["activation_s"]
    if not math.isfinite(activation):
        continue
    active_index = next(i for i, row in enumerate(history) if int(row["active"]) == 1)
    selections = [("pre", history[active_index-1]), ("activation", history[active_index])]
    selections += [(f"plus_{delta:g}s", closest(history, activation+delta)) for delta in TARGETS]
    bmap = {row["time_s"]: row for row in case["boundary"]}
    for label, hrow in selections:
        if hrow is None:
            jump_rows.append({"case": name, "sample": label, "status": "CENSORED_AFTER_FAILURE"})
            continue
        time = hrow["time_s"]
        dgroup = case["diagnostics"].get(time, [])
        bottom = dgroup[0] if dgroup else None
        top = dgroup[-1] if dgroup else None
        cinit = case["donor"]["c_initial"] if case["donor"] else math.nan
        row = {
            "case": name, "sample": label, "status": "ACCEPTED", "time_s": time,
            "elapsed_after_activation_s": time-activation, "ceq": CEQ,
            "c_ale_initial": cinit, "soil_top_sw": hrow["top_soil_sw"],
            "ale_bottom_porosity": bottom["epsilon"] if bottom else math.nan,
            "ale_top_porosity": top["epsilon"] if top else math.nan,
            "ale_height_m": hrow["h_m"],
            "interface_total_liquid_flux_kg_m2_s": bottom["bottom_liquid_kg_m2_s"] if bottom else math.nan,
            "interface_water_flux_kg_m2_s": (bottom["bottom_liquid_kg_m2_s"]-bottom["bottom_salt_kg_m2_s"]) if bottom else math.nan,
            "interface_salt_flux_kg_m2_s": bottom["bottom_salt_kg_m2_s"] if bottom else math.nan,
            "interface_advective_salt_flux_kg_m2_s": bottom["bottom_adv_kg_m2_s"] if bottom else math.nan,
            "interface_diffusive_salt_flux_kg_m2_s": bottom["bottom_diff_kg_m2_s"] if bottom else math.nan,
            "ale_bottom_reaction_kg_m3_s": bottom["reaction_kg_m3_s"] if bottom else math.nan,
        }
        row.update(concentration_fields("c_soil_top_dof", hrow["top_soil_mass_fraction"]))
        row.update(concentration_fields("c_interface", bottom["bottom_c"] if bottom else math.nan))
        row.update(concentration_fields("c_ale_bottom_cell", bottom["c"] if bottom else math.nan))
        row.update(concentration_fields("c_ale_top_cell", top["c"] if top else math.nan))
        jump_rows.append(row)

        if label in ("pre", "activation") or not bottom:
            continue
        index = history.index(hrow)
        previous = history[index-1]
        dt = hrow["time_s"] - previous["time_s"]
        boundary = bmap[time]
        soil_salt_rate = (hrow["soil_salt_kg"]-previous["soil_salt_kg"])/dt
        soil_liquid_rate = ((hrow["soil_salt_kg"]+hrow["soil_water_kg"])
                            -(previous["soil_salt_kg"]+previous["soil_water_kg"]))/dt
        soil_salt_out = (boundary["bottom_salt_in_kg_s"]-soil_salt_rate)/AREA
        soil_liquid_out = (boundary["bottom_liquid_in_kg_s"]-soil_liquid_rate)/AREA
        # The two interface ledger columns are BDF-integrated conserved fluxes.
        # Their accepted-step differences, rather than the instantaneous flux
        # evaluated at the final state, are the quantities imposed on DuMuX.
        accepted_liquid = (hrow["interface_liquid_kg_m2"]
                           - previous["interface_liquid_kg_m2"])/dt
        accepted_salt = (hrow["interface_salt_kg_m2"]
                         - previous["interface_salt_kg_m2"])/dt
        distance = 0.5*(bottom["z_bottom_m"]+bottom["z_top_m"])
        deff = math.sqrt(bottom["epsilon"])*DIFFUSION
        independent_diff = -RHO*bottom["epsilon"]*deff*(bottom["c"]-bottom["bottom_c"])/distance
        independent_adv = bottom["bottom_c"]*bottom["bottom_liquid_kg_m2_s"]
        independent_salt = independent_adv + independent_diff
        flux_rows.append({
            "case": name, "sample": label, "time_s": time, "elapsed_after_activation_s": time-activation,
            "soil_storage_reconstructed_liquid_out_kg_m2_s": soil_liquid_out,
            "soil_storage_reconstructed_salt_out_kg_m2_s": soil_salt_out,
            "ale_accepted_integrated_liquid_in_kg_m2_s": accepted_liquid,
            "ale_accepted_integrated_salt_in_kg_m2_s": accepted_salt,
            "ale_recorded_liquid_in_kg_m2_s": bottom["bottom_liquid_kg_m2_s"],
            "ale_recorded_advective_salt_in_kg_m2_s": bottom["bottom_adv_kg_m2_s"],
            "ale_recorded_diffusive_salt_in_kg_m2_s": bottom["bottom_diff_kg_m2_s"],
            "ale_recorded_total_salt_in_kg_m2_s": bottom["bottom_salt_kg_m2_s"],
            "independent_ale_advective_salt_kg_m2_s": independent_adv,
            "independent_ale_diffusive_salt_kg_m2_s": independent_diff,
            "independent_ale_total_salt_kg_m2_s": independent_salt,
            "soil_minus_ale_liquid_kg_m2_s": soil_liquid_out-accepted_liquid,
            "soil_minus_ale_salt_kg_m2_s": soil_salt_out-accepted_salt,
            "instantaneous_minus_step_average_liquid_kg_m2_s": bottom["bottom_liquid_kg_m2_s"]-accepted_liquid,
            "instantaneous_minus_step_average_salt_kg_m2_s": bottom["bottom_salt_kg_m2_s"]-accepted_salt,
            "operator_rebuild_minus_recorded_salt_kg_m2_s": independent_salt-bottom["bottom_salt_kg_m2_s"],
            "net_salt_fraction": bottom["bottom_salt_kg_m2_s"]/bottom["bottom_liquid_kg_m2_s"],
            "interface_mass_fraction": bottom["bottom_c"],
            "ale_bottom_cell_mass_fraction": bottom["c"],
            "boundary_to_cell_distance_m": distance,
            "ale_bottom_deff_m2_s": deff,
        })
write_csv(OUT / "ACTIVATION_JUMP.csv", jump_rows)
write_csv(OUT / "INTERFACE_FLUX_AUDIT.csv", flux_rows)


counter_rows = []
for directory in sorted((ROOT / "counterfactual").glob("*")):
    if not (directory / "history.csv").exists():
        continue
    metadata = read_csv(directory / "metadata.csv")[0]
    history = read_csv(directory / "history.csv")
    status_lines = (directory / "status.txt").read_text().splitlines()
    for target in (0.0, 1.0, 10.0, 100.0, 500.0):
        row = closest(history, target, tolerance=0.251)
        if row is None:
            continue
        counter_rows.append({
            "case": directory.name, "sample_time_s": row["time_s"], "run_status": status_lines[0],
            "cells": metadata["cells"], "mode": metadata["mode"], "ceq": metadata["ceq"],
            "c_initial": metadata["c_initial"], "c_bottom": metadata["c_bottom"],
            "c_bottom_minus_ceq": metadata["c_bottom"]-metadata["ceq"],
            "solid_adjustment_kg_m2": metadata["solid_adjustment_kg_m2"],
            "soil_water_return_kg_m2": metadata["soil_water_return_kg_m2"],
            "height_m": row["h_m"], "epsilon_min": row["epsilon_min"], "epsilon_max": row["epsilon_max"],
            "c_bottom_cell": row["c_bottom_cell"], "c_top_cell": row["c_top_cell"],
            "cumulative_internal_positive_kg_m2": row["cumulative_internal_positive_kg_m2"],
            "cumulative_internal_dissolution_kg_m2": row["cumulative_internal_dissolution_kg_m2"],
            "cumulative_bottom_salt_kg_m2": row["cumulative_bottom_salt_kg_m2"],
            "bottom_liquid_kg_m2_s": row["bottom_liquid_kg_m2_s"],
            "bottom_advective_salt_kg_m2_s": row["bottom_advective_salt_kg_m2_s"],
            "bottom_diffusive_salt_kg_m2_s": row["bottom_diffusive_salt_kg_m2_s"],
            "bottom_salt_kg_m2_s": row["bottom_salt_kg_m2_s"],
        })
write_csv(OUT / "COUNTERFACTUAL_RESULTS.csv", counter_rows)


mpl.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8, "axes.labelsize": 8,
    "axes.titlesize": 9, "legend.fontsize": 7, "xtick.labelsize": 7,
    "ytick.labelsize": 7, "pdf.fonttype": 42, "ps.fonttype": 42,
    "svg.fonttype": "none", "axes.spines.top": False, "axes.spines.right": False,
})
COLORS = {"h0p1": "#0072B2", "h3p5": "#D55E00", "adv": "#009E73", "diff": "#CC79A7", "net": "#000000"}


def save(fig, stem):
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(FIG / f"{stem}.{suffix}", dpi=400 if suffix == "png" else None,
                    facecolor="white", transparent=False)
    plt.close(fig)


# Figure 1: accepted concentration histories, with failure censoring retained.
fig, axes = plt.subplots(1, 2, figsize=(7.1, 2.75), layout="constrained", sharey=True)
for ax, name, title, color in zip(axes, audit_names, ("0.1 µm activation", "3.5 µm activation"),
                                  (COLORS["h0p1"], COLORS["h3p5"])):
    case = cases[name]
    activation = case["summary"]["activation_s"]
    diag = case["diagnostics"]
    times = sorted(t for t in diag if 0 <= t-activation <= 500)
    interface = [diag[t][0]["bottom_c"] for t in times]
    bottom_cell = [diag[t][0]["c"] for t in times]
    history_map = {r["time_s"]: r for r in case["history"]}
    soil = [history_map[t]["top_soil_mass_fraction"] for t in times]
    elapsed = np.array(times)-activation
    ax.plot(elapsed, 1e6*(np.array(soil)-CEQ), color=color, label="DuMuX top DOF")
    ax.plot(elapsed, 1e6*(np.array(interface)-CEQ), color="black", ls="--", label="ALE boundary value")
    ax.plot(elapsed, 1e6*(np.array(bottom_cell)-CEQ), color="#009E73", ls=":", label="ALE cell 0")
    ax.axhline(0, color="0.65", lw=0.8)
    ax.set_yscale("symlog", linthresh=0.5, linscale=0.8)
    ax.set(xlabel="Time after activation (s)", title=title)
axes[0].set_ylabel(r"$(C-C_{eq})\times10^6$ (mass fraction)")
axes[0].legend(frameon=False)
save(fig, "figure1_activation_concentration")


# Figure 2: soil and crust concentration profiles for the 3.5 µm failure.
case = cases["x008_reservoir_n80_a64"]
activation = case["summary"]["activation_s"]
soil_rows = read_csv(case["dir"] / "soil_fields.csv")
available_soil_times = sorted(set(r["time_s"] for r in soil_rows if r["time_s"] >= activation))
profile_targets = [activation, activation+1, activation+100, min(case["diagnostics"])]
profile_targets[-1] = max(case["diagnostics"])
selected_soil = [min(available_soil_times, key=lambda t: abs(t-target)) for target in profile_targets]
fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.0), layout="constrained")
for target, soil_time in zip(profile_targets, selected_soil):
    label = f"+{target-activation:.1f} s"
    sr = [r for r in soil_rows if r["time_s"] == soil_time]
    top_y = max(r["y_m"] for r in sr)
    sr = [r for r in sr if top_y-r["y_m"] <= 0.001]
    axes[0].plot([r["mass_fraction"] for r in sr], 1e3*np.array([top_y-r["y_m"] for r in sr]), label=label)
    dt = min(case["diagnostics"], key=lambda t: abs(t-target))
    dr = case["diagnostics"][dt]
    z = [0.5*(r["z_bottom_m"]+r["z_top_m"])*1e6 for r in dr]
    axes[1].plot([r["c"] for r in dr], z, label=label)
axes[0].invert_yaxis()
axes[0].axvline(CEQ, color="0.6", lw=0.8, ls="--")
axes[1].axvline(CEQ, color="0.6", lw=0.8, ls="--")
axes[0].set(xlabel="NaCl mass fraction", ylabel="Depth below soil surface (mm)", title="Soil")
axes[1].set(xlabel="NaCl mass fraction", ylabel="Height above interface (µm)", title="Porous salt crust")
axes[1].legend(frameon=False)
save(fig, "figure2_soil_crust_concentration_profiles")


# Figure 3: interface salt-flux decomposition.
fig, axes = plt.subplots(1, 2, figsize=(7.1, 2.8), layout="constrained", sharey=True)
for ax, name, title in zip(axes, audit_names, ("0.1 µm", "3.5 µm")):
    case = cases[name]
    activation = case["summary"]["activation_s"]
    times = sorted(t for t in case["diagnostics"] if 0 < t-activation <= 500)
    bottom = [case["diagnostics"][t][0] for t in times]
    elapsed = np.array(times)-activation
    ax.plot(elapsed, [r["bottom_adv_kg_m2_s"] for r in bottom], color=COLORS["adv"], label="Advective")
    ax.plot(elapsed, [r["bottom_diff_kg_m2_s"] for r in bottom], color=COLORS["diff"], ls="--", label="Diffusive")
    ax.plot(elapsed, [r["bottom_salt_kg_m2_s"] for r in bottom], color=COLORS["net"], lw=1.2, label="Net salt")
    ax.axhline(0, color="0.65", lw=0.8)
    ax.set(xlabel="Time after activation (s)", title=title)
axes[0].set_ylabel(r"Interface salt flux (kg m$^{-2}$ s$^{-1}$)")
axes[0].legend(frameon=False)
save(fig, "figure3_interface_flux_decomposition")


# Figure 4: signed reaction profiles in the 3.5 µm failure.
fig, ax = plt.subplots(figsize=(3.65, 3.0), layout="constrained")
for delta in (1.0, 100.0, 300.0, 497.0):
    time = min(case["diagnostics"], key=lambda t: abs(t-(activation+delta)))
    rows = case["diagnostics"][time]
    z = [0.5*(r["z_bottom_m"]+r["z_top_m"])*1e6 for r in rows]
    ax.plot([r["reaction_kg_m3_s"] for r in rows], z, label=f"+{time-activation:.1f} s")
ax.axvline(0, color="0.65", lw=0.8)
ax.set(xlabel=r"Internal reaction (kg m$^{-3}$ s$^{-1}$)", ylabel="Height (µm)",
       title="3.5 µm case: precipitation (+), dissolution (−)")
ax.legend(frameon=False)
save(fig, "figure4_internal_reaction_profiles")


# Figure 5: prescribed-bottom-concentration sensitivity.
fig, axes = plt.subplots(1, 2, figsize=(7.1, 2.8), layout="constrained")
cf_dirs = sorted((ROOT / "counterfactual").glob("actual_mix_cb_*_a64"))
for directory in cf_dirs:
    rows = read_csv(directory / "history.csv")
    meta = read_csv(directory / "metadata.csv")[0]
    delta = meta["c_bottom"]-CEQ
    label = f"ΔC={delta:+.0e}"
    axes[0].plot([r["time_s"] for r in rows], [r["epsilon_max"] for r in rows], label=label)
    axes[1].plot([r["time_s"] for r in rows], 1e3*np.array([r["cumulative_internal_dissolution_kg_m2"] for r in rows]), label=label)
for ax in axes:
    ax.set_xlabel("Time (s)")
axes[0].axhline(0.95, color="0.5", ls="--", lw=0.8)
axes[0].set(ylabel=r"Maximum porosity, $\epsilon_{max}$", title="Porosity response")
axes[1].set(ylabel=r"Cumulative dissolution (g m$^{-2}$)", title="Dissolution amplification")
axes[0].legend(frameon=False, ncol=2)
save(fig, "figure5_bottom_concentration_sensitivity")


# Figure 6: matched old/new soil conditions at +300 s after activation.
labels, dissolution, positive, height, epsmax = [], [], [], [], []
for name in ("x09999_noflux_n80_a64", "x09999_reservoir_n80_a64", "x008_noflux_n80_a64", "x008_reservoir_n80_a64"):
    c = cases[name]
    act = c["summary"]["activation_s"]
    label = (("x=.09999" if name.startswith("x09999") else "x=.08")
             + "\n" + ("no_flux" if "_noflux_" in name else "reservoir"))
    labels.append(label)
    if not math.isfinite(act):
        dissolution.append(math.nan); positive.append(math.nan); height.append(math.nan); epsmax.append(math.nan)
        continue
    row = closest(c["history"], act+300)
    group = c["diagnostics"].get(row["time_s"], []) if row else []
    dissolution.append(1e3*row["internal_dissolution_kg_m2"] if row else math.nan)
    positive.append(1e3*row["internal_positive_kg_m2"] if row else math.nan)
    height.append(1e6*row["h_m"] if row else math.nan)
    epsmax.append(max((r["epsilon"] for r in group), default=math.nan))
fig, axes = plt.subplots(2, 2, figsize=(7.1, 5.0), layout="constrained")
x = np.arange(len(labels))
for ax, values, ylabel in zip(axes.flat, (dissolution, positive, height, epsmax),
                              ("Dissolution (g m$^{-2}$)", "Precipitation (g m$^{-2}$)",
                               "Crust height (µm)", r"$\epsilon_{max}$")):
    ax.bar(x, np.nan_to_num(values, nan=0.0), color=["#0072B2", "#56B4E9", "#E69F00", "#D55E00"])
    for xi, value in enumerate(values):
        if not math.isfinite(value):
            ax.text(xi, 0, "no activation", rotation=90, va="bottom", ha="center", fontsize=7)
    ax.set_ylabel(ylabel)
    ax.set_xticks(x, labels, rotation=0, ha="center")
axes[0, 0].set_yscale("log")
axes[0, 0].set_ylim(1e-4, max(v for v in dissolution if math.isfinite(v))*1.4)
fig.suptitle("Matched comparison at +300 s after accepted activation")
save(fig, "figure6_matched_case_comparison")


verification = {
    "ceq_mass_fraction": CEQ,
    "coupled_cases": len(summaries),
    "counterfactual_cases": len(set(r["case"] for r in counter_rows)),
    "activation_jump_rows": len(jump_rows),
    "interface_flux_rows": len(flux_rows),
    "max_operator_rebuild_abs_error_kg_m2_s": max(abs(r["operator_rebuild_minus_recorded_salt_kg_m2_s"]) for r in flux_rows),
    "max_soil_ale_salt_abs_error_kg_m2_s": max(abs(r["soil_minus_ale_salt_kg_m2_s"]) for r in flux_rows),
    "max_soil_ale_liquid_abs_error_kg_m2_s": max(abs(r["soil_minus_ale_liquid_kg_m2_s"]) for r in flux_rows),
    "figures": sorted(p.name for p in FIG.glob("*")),
    "matplotlib_version": mpl.__version__,
    "numpy_version": np.__version__,
    "slurm_job_id": os.environ.get("SLURM_JOB_ID"),
    "hostname": os.uname().nodename,
}
(OUT / "verification.json").write_text(json.dumps(verification, indent=2) + "\n")

raw_hashes = []
raw_sizes = []
raw_names = {
    "case_spec.json", "command.json", "outcome.json", "exit_status.txt", "solver.out",
    "coupled_history.csv", "ale_diagnostics.csv", "open_boundary_budget.csv",
    "grid_profiles.csv", "donor_attempts.csv", "donor_profiles.csv", "coupling_events.csv",
}
for directory in case_dirs:
    for path in sorted(directory.iterdir()):
        if path.is_file() and path.name in raw_names:
            relative = path.relative_to(ROOT).as_posix()
            raw_hashes.append(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {relative}")
            raw_sizes.append({"path": relative, "bytes": path.stat().st_size})
for path in sorted((ROOT / "counterfactual").rglob("*")):
    if path.is_file():
        relative = path.relative_to(ROOT).as_posix()
        raw_hashes.append(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {relative}")
        raw_sizes.append({"path": relative, "bytes": path.stat().st_size})
(OUT / "raw_output_sha256.txt").write_text("\n".join(raw_hashes) + "\n")
write_csv(OUT / "raw_output_sizes.csv", raw_sizes)

provenance = {
    "branch": "research/p2-h-interface-concentration-audit",
    "base_commit": "1fe9d1f7f9fee8094b671367992d463a8e8b96c9",
    "historical_commits": {
        "av_sensitivity": "63f267f9eae85bf434a2ad9ed0804a63723e9478",
        "x008_h0p1": "b2e1797a54a68f56ea6e8c49c1ad72360f847ba9",
        "x008_h3p5": "1fe9d1f7f9fee8094b671367992d463a8e8b96c9",
    },
    "hpcc_root": str(ROOT),
    "soil_solver_sha256": hashlib.sha256((ROOT / "build" / "soil_solver").read_bytes()).hexdigest(),
    "coupled_jobs": sorted({str(s["slurm_job_id"]) for s in summaries}),
    "counterfactual_array_job": "16438812",
    "counterfactual_build_job": "16438804",
    "coupled_array_jobs": ["16438805", "16438821"],
    "postprocess_job": os.environ.get("SLURM_JOB_ID"),
    "raw_fields_retained_on_hpcc": True,
    "raw_fields_committed_to_git": False,
}
(OUT / "provenance.json").write_text(json.dumps(provenance, indent=2) + "\n")

hashes = []
for path in sorted(OUT.rglob("*")):
    if path.is_file() and path.name != "CHECKSUMS.sha256":
        hashes.append(f"{hashlib.sha256(path.read_bytes()).hexdigest()}  {path.relative_to(OUT).as_posix()}")
(OUT / "CHECKSUMS.sha256").write_text("\n".join(hashes) + "\n")
print(json.dumps(verification, indent=2))
