"""Audit M1/M2/M3 accepted histories and create the registered Step-1 deliverables."""
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

if not os.environ.get("SLURM_JOB_ID") or os.uname().nodename.startswith(("ibcln", "login", "test")):
    raise RuntimeError("Slurm compute allocation required")

root = Path.cwd()
out = root / "deliverable"
figdir = out / "figures"
records = out / "records"
out.mkdir(exist_ok=False)
figdir.mkdir()
records.mkdir()

CEQ = 0.1*0.05844/(0.9*0.01801528 + 0.1*0.05844)
MODELS = ["M1", "M2", "M3"]
COLORS = {"M1": "#0072B2", "M2": "#D55E00", "M3": "#009E73"}
STYLES = {"M1": "-", "M2": "--", "M3": "-."}

def digest(path):
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024*1024), b""):
            h.update(block)
    return h.hexdigest()

def load_case(name):
    directory = root / "runs" / name
    data = np.genfromtxt(directory / "history.csv", delimiter=",", names=True, dtype=None, encoding="utf-8")
    outcome = json.loads((directory / "outcome.json").read_text())
    rejected = list(csv.DictReader((directory / "rejected_trials.csv").open()))
    events = list(csv.DictReader((directory / "events.csv").open()))
    return {"dir": directory, "data": data, "outcome": outcome, "rejected": rejected, "events": events}

def write_table(path, rows):
    if not rows:
        path.write_text("")
        return
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

def value_at(data, field, time):
    return float(np.interp(time, data["time_s"], data[field]))

cases = {model: load_case(model.lower()+"_dt0p5") for model in MODELS}
half = {}
for model in ["M2", "M3"]:
    path = root / "runs" / (model.lower()+"_dt0p25")
    if path.exists():
        half[model] = load_case(path.name)

# Exact pre-activation identity, excluding the intentionally different growth-law label.
pre_fields = ["time_s", "dt_s", "active", "activation_s", "h_m", "solid_kg_m2", "pore_salt_kg_m2",
              "pore_water_kg_m2", "j_kg_m2_s", "v_m_s", "c_soil", "sw_soil", "soil_salt_kg",
              "soil_water_kg", "vapour_kg", "bottom_salt_kg", "bottom_water_kg", "precip_flux", "fill_liquid_flux"]
activation = min(float(cases[m]["outcome"]["activation_time_s"]) for m in MODELS)
pre_regression = {}
base = cases["M1"]["data"]
for model in ["M2", "M3"]:
    other = cases[model]["data"]
    a = base[base["time_s"] <= activation]
    b = other[other["time_s"] <= activation]
    assert len(a) == len(b)
    pre_regression[model] = {field: float(np.max(np.abs(a[field].astype(float)-b[field].astype(float)))) for field in pre_fields}

common_last = min(float(cases[m]["data"]["time_s"][-1]) for m in MODELS)
relative_times = [0.0, 10.0, 100.0, 300.0, 500.0]
times = sorted(set([activation+x for x in relative_times if activation+x <= common_last] + [common_last]))
comparison_rows = []
fields = ["h_m", "v_m_s", "j_kg_m2_s", "c_top", "c_soil", "c_soil_relative_to_ceq",
          "total_out_kg_m2_s", "salt_out_kg_m2_s", "solid_rate_kg_m2_s", "pore_salt_rate_kg_m2_s",
          "pore_water_rate_kg_m2_s", "da", "da_star", "da_top", "efficiency_m1", "efficiency_m2",
          "efficiency_m3", "selected_efficiency", "m3_relative_error"]
for time in times:
    for model in MODELS:
        data = cases[model]["data"]
        row = {"comparison": "common_time", "model": model, "time_s": time, "time_since_activation_s": time-activation,
               "height_target_um": "", "interpolated": int(not np.any(np.isclose(data["time_s"], time, rtol=0, atol=1e-10)))}
        row.update({field: value_at(data, field, time) for field in fields})
        comparison_rows.append(row)

min_final_height = min(float(cases[m]["data"]["h_m"][-1]) for m in MODELS)
height_targets = np.linspace(3.5e-6, min_final_height, 6)
for height in height_targets:
    for model in MODELS:
        data = cases[model]["data"]
        active = data[data["active"] == 1]
        row = {"comparison": "common_height", "model": model,
               "time_s": float(np.interp(height, active["h_m"], active["time_s"])), "time_since_activation_s": "",
               "height_target_um": height*1e6, "interpolated": int(not np.any(np.isclose(active["h_m"], height, rtol=0, atol=1e-14)))}
        row.update({field: float(np.interp(height, active["h_m"], active[field])) for field in fields})
        comparison_rows.append(row)
write_table(out / "MODEL_COMPARISON.csv", comparison_rows)

activation_rows = []
for model in MODELS:
    data = cases[model]["data"]
    idx = int(np.where(data["active"] == 1)[0][0])
    post = idx+1
    activation_rows.append({
        "model": model,
        "first_precipitation_time_s": next(float(e["time_s"]) for e in cases[model]["events"] if e["event"] == "accepted_first_precipitation"),
        "activation_time_s": float(data["time_s"][idx]),
        "activation_height_um": float(data["h_m"][idx])*1e6,
        "activation_solid_kg_m2": float(data["solid_kg_m2"][idx]),
        "activation_pore_salt_kg_m2": float(data["pore_salt_kg_m2"][idx]),
        "activation_pore_water_kg_m2": float(data["pore_water_kg_m2"][idx]),
        "activation_c_soil": float(data["c_soil"][idx]),
        "salt_flux_before_kg_m2_s": float(data["salt_out_kg_m2_s"][idx]),
        "salt_flux_first_growth_kg_m2_s": float(data["salt_out_kg_m2_s"][post]),
        "salt_flux_after_over_before": float(data["salt_out_kg_m2_s"][post]/data["salt_out_kg_m2_s"][idx]),
        "liquid_flux_first_growth_kg_m2_s": float(data["total_out_kg_m2_s"][post]),
        "solid_rate_first_growth_kg_m2_s": float(data["solid_rate_kg_m2_s"][post]),
        "pore_salt_rate_first_growth_kg_m2_s": float(data["pore_salt_rate_kg_m2_s"][post]),
        "pore_water_rate_first_growth_kg_m2_s": float(data["pore_water_rate_kg_m2_s"][post]),
        "da_activation": float(data["da"][idx]), "da_star_activation": float(data["da_star"][idx]),
        "da_top_activation": float(data["da_top"][idx]), "selected_efficiency_activation": float(data["selected_efficiency"][idx]),
    })
write_table(out / "ACTIVATION_FLUX_COMPARISON.csv", activation_rows)

budget_rows = []
validity_rows = []
for model in MODELS:
    case = cases[model]
    data = case["data"]
    outcome = case["outcome"]
    active = data[data["active"] == 1]
    final = data[-1]
    budget_rows.append({
        "model": model, "dt_s": 0.5, "requested_time_s": 9000, "last_accepted_time_s": float(final["time_s"]),
        "solver_exit_code": outcome["solver_exit_code"], "completed_requested_horizon": outcome["completed_requested_horizon"],
        "final_height_um": float(final["h_m"])*1e6, "final_solid_g_m2": float(final["solid_kg_m2"])*1e3,
        "final_pore_salt_g_m2": float(final["pore_salt_kg_m2"])*1e3, "final_pore_water_g_m2": float(final["pore_water_kg_m2"])*1e3,
        "max_salt_residual_relative": outcome["max_salt_residual_relative"],
        "max_water_residual_relative": outcome["max_water_residual_relative"],
        "job_id": outcome["slurm_job_id"], "node": outcome["hostname"], "history_sha256": outcome["history_sha256"],
    })
    crossing = {}
    deficit = -active["c_soil_relative_to_ceq"]
    for threshold in [0.005, 0.01, 0.02]:
        hits = np.where(deficit >= threshold)[0]
        crossing[str(threshold)] = float(active["time_s"][hits[0]]) if len(hits) else None
    rejected = case["rejected"]
    validity_rows.append({
        "model": model, "minimum_da_active": float(np.min(active["da"])), "maximum_da_active": float(np.max(active["da"])),
        "da_star": float(active["da_star"][0]), "all_m3_steps_da_ge_10": bool(np.all(active["m3_valid"] == 1)),
        "maximum_abs_m3_relative_error": float(np.max(np.abs(active["m3_relative_error"]))),
        "maximum_abs_soil_relative_deviation": float(np.max(np.abs(active["c_soil_relative_to_ceq"]))),
        "undersaturation_0p5pct_time_s": crossing["0.005"], "undersaturation_1pct_time_s": crossing["0.01"],
        "undersaturation_2pct_accepted_time_s": crossing["0.02"],
        "first_rejected_time_s": float(rejected[0]["time_s"]) if rejected else None,
        "first_rejected_reason": rejected[0]["reason"] if rejected else None,
        "first_rejected_c_soil": float(rejected[0]["c_soil"]) if rejected else None,
        "completed_requested_horizon": outcome["completed_requested_horizon"],
    })
for model, case in half.items():
    final = case["data"][-1]
    outcome = case["outcome"]
    budget_rows.append({
        "model": model, "dt_s": 0.25, "requested_time_s": 9000, "last_accepted_time_s": float(final["time_s"]),
        "solver_exit_code": outcome["solver_exit_code"], "completed_requested_horizon": outcome["completed_requested_horizon"],
        "final_height_um": float(final["h_m"])*1e6, "final_solid_g_m2": float(final["solid_kg_m2"])*1e3,
        "final_pore_salt_g_m2": float(final["pore_salt_kg_m2"])*1e3, "final_pore_water_g_m2": float(final["pore_water_kg_m2"])*1e3,
        "max_salt_residual_relative": outcome["max_salt_residual_relative"], "max_water_residual_relative": outcome["max_water_residual_relative"],
        "job_id": outcome["slurm_job_id"], "node": outcome["hostname"], "history_sha256": outcome["history_sha256"],
    })
write_table(out / "MASS_BUDGET.csv", budget_rows)
write_table(out / "VALIDITY_CHECKS.csv", validity_rows)

# Numerical and timestep verification.
m1_regression = json.loads((root / "unit/m1_regression.json").read_text())
timestep = {}
for model, case in half.items():
    coarse = cases[model]["data"]
    fine = case["data"]
    t = min(float(coarse["time_s"][-1]), float(fine["time_s"][-1]))
    timestep[model] = {field: abs(value_at(coarse, field, t)-value_at(fine, field, t))/max(abs(value_at(fine, field, t)), 1e-300)
                       for field in ["h_m", "v_m_s", "j_kg_m2_s", "c_soil", "salt_out_kg_m2_s"]}
verification = {
    "numerical_verification": "PASS" if m1_regression["status"] == "PASS" and
        max(max(row["max_salt_residual_relative"], row["max_water_residual_relative"]) for row in budget_rows) <= 1e-8 and
        max(max(values.values()) for values in pre_regression.values()) == 0 else "FAIL",
    "m1_exact_regression": m1_regression,
    "preactivation_exact_regression": pre_regression,
    "timestep_relative_differences_at_common_last_time": timestep,
    "all_m3_steps_da_ge_10": all(row["all_m3_steps_da_ge_10"] for row in validity_rows),
    "formula_replacement_only": True,
    "does_not_implement_internal_reaction_mass_partition": True,
    "postprocess_job_id": os.environ["SLURM_JOB_ID"],
    "postprocess_node": os.uname().nodename,
}
(out / "verification.json").write_text(json.dumps(verification, indent=2))

# Figures: full accepted histories only; stopped curves end at actual accepted time.
plt.rcParams.update({"font.size": 9, "pdf.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False})
def save(fig, name):
    fig.savefig(figdir/(name+".png"), dpi=300, facecolor="white")
    fig.savefig(figdir/(name+".pdf"), facecolor="white")
    plt.close(fig)

heights = np.linspace(3.5e-6, max(float(cases[m]["data"]["h_m"][-1]) for m in MODELS), 400)
epsilon, ds, kr, beta, dsub = 0.43, 1.3e-9, 0.001, 2.0, 25e-6
deff = np.sqrt(epsilon)*ds
av = 36*(1-epsilon)/dsub
da = heights**2*(1-CEQ)*av*kr/(epsilon*deff)
da_star = beta**2*(1-epsilon)**2*(1-CEQ)*kr/(epsilon*deff*av)
da_top = beta*(1-epsilon)*(1-CEQ)*heights*kr/(epsilon*deff)
efficiencies = {
    "M1": da_top/(1+da_top),
    "M2": np.sqrt(da_star)*np.tanh(np.sqrt(da))/(1+np.sqrt(da_star)*np.tanh(np.sqrt(da))),
    "M3": np.full_like(heights, np.sqrt(da_star)/(1+np.sqrt(da_star))),
}
fig, ax = plt.subplots(figsize=(7, 4.2), layout="constrained")
for model in MODELS:
    ax.plot(heights*1e6, efficiencies[model], color=COLORS[model], ls=STYLES[model], label=model)
ax.set(xlabel="Crust height (μm)", ylabel="Growth efficiency (-)")
ax.grid(alpha=.2); ax.legend()
save(fig, "figure1_efficiency_vs_height")

fig, axes = plt.subplots(2, 1, figsize=(7, 6), sharex=True, layout="constrained")
for model in MODELS:
    z = cases[model]["data"]
    axes[0].plot(z["time_s"], z["h_m"]*1e6, color=COLORS[model], ls=STYLES[model], label=model)
    axes[1].plot(z["time_s"], z["v_m_s"]*1e9, color=COLORS[model], ls=STYLES[model])
    axes[0].plot(z["time_s"][-1], z["h_m"][-1]*1e6, "o", color=COLORS[model], ms=4)
axes[0].set_ylabel("Crust height (μm)"); axes[0].legend(); axes[1].set(xlabel="Accepted physical time (s)", ylabel="Top growth velocity (nm/s)")
for ax in axes: ax.grid(alpha=.2)
save(fig, "figure2_height_and_velocity")

fig, axes = plt.subplots(2, 1, figsize=(7, 6), sharex=True, layout="constrained")
for model in MODELS:
    z = cases[model]["data"]
    axes[0].plot(z["time_s"], 100*z["c_soil"], color=COLORS[model], ls=STYLES[model], label=model)
    axes[1].plot(z["time_s"], -100*z["c_soil_relative_to_ceq"], color=COLORS[model], ls=STYLES[model])
axes[0].axhline(100*CEQ, color="black", ls=":", lw=.8, label="Ceq")
axes[1].axhline(2, color="black", ls=":", lw=.8, label="2% validity limit")
axes[0].set_ylabel("Soil-interface NaCl mass fraction (%)"); axes[0].legend()
axes[1].set(xlabel="Accepted physical time (s)", ylabel="Relative undersaturation (%)"); axes[1].legend()
for ax in axes: ax.grid(alpha=.2)
save(fig, "figure3_interface_concentration")

fig, axes = plt.subplots(2, 1, figsize=(7, 6), sharex=True, layout="constrained")
for model in MODELS:
    z = cases[model]["data"]
    axes[0].plot(z["time_s"], z["salt_out_kg_m2_s"]*1e6, color=COLORS[model], ls=STYLES[model], label=model)
    axes[1].plot(z["time_s"], z["solid_rate_kg_m2_s"]*1e6, color=COLORS[model], ls=STYLES[model], label=model+" solid")
    axes[1].plot(z["time_s"], z["pore_salt_rate_kg_m2_s"]*1e6, color=COLORS[model], ls=":", alpha=.8, label=model+" pore salt")
axes[0].set_ylabel("Soil-to-crust NaCl flux (mg m⁻² s⁻¹)"); axes[0].legend(ncol=3)
axes[1].set(xlabel="Accepted physical time (s)", ylabel="New-volume NaCl demand (mg m⁻² s⁻¹)"); axes[1].legend(ncol=2, fontsize=8)
for ax in axes: ax.grid(alpha=.2)
save(fig, "figure4_salt_supply_and_demand")

fig, axes = plt.subplots(2, 1, figsize=(7, 6), sharex=True, layout="constrained")
for model in MODELS:
    z = cases[model]["data"]
    axes[0].plot(z["time_s"], z["j_kg_m2_s"]*1e6, color=COLORS[model], ls=STYLES[model], label=model)
    axes[1].plot(z["time_s"], 100*z["c_top"], color=COLORS[model], ls=STYLES[model])
axes[0].set_ylabel("Evaporation flux (mg m⁻² s⁻¹)"); axes[0].legend()
axes[1].set(xlabel="Accepted physical time (s)", ylabel="Crust-top NaCl mass fraction (%)")
for ax in axes: ax.grid(alpha=.2)
save(fig, "figure5_evaporation_and_top_concentration")

fig, axes = plt.subplots(2, 1, figsize=(7, 6), sharex=True, layout="constrained")
for model in MODELS:
    z = cases[model]["data"]
    axes[0].semilogy(z["time_s"], np.maximum(np.abs(z["salt_residual_rel"]), 1e-18), color=COLORS[model], ls=STYLES[model], label=model)
    axes[1].semilogy(z["time_s"], np.maximum(np.abs(z["water_residual_rel"]), 1e-18), color=COLORS[model], ls=STYLES[model])
for ax in axes: ax.axhline(1e-8, color="black", ls=":", lw=.8); ax.grid(alpha=.2)
axes[0].set_ylabel("Absolute relative salt residual"); axes[0].legend()
axes[1].set(xlabel="Accepted physical time (s)", ylabel="Absolute relative water residual")
save(fig, "figure6_conservation")

fig, axes = plt.subplots(2, 1, figsize=(7, 6), sharex=True, layout="constrained")
for model in MODELS:
    z = cases[model]["data"]; active = z[z["active"] == 1]
    axes[0].plot(active["time_s"], active["da"], color=COLORS[model], ls=STYLES[model], label=model+" Da")
    axes[0].plot(active["time_s"], active["da_top"], color=COLORS[model], ls=":", alpha=.75)
    axes[1].plot(active["time_s"], 100*active["m3_relative_error"], color=COLORS[model], ls=STYLES[model], label=model)
axes[0].axhline(10, color="black", ls="--", lw=.8, label="M3 Da threshold")
axes[0].set_ylabel("Da (solid) and Da_top (dotted)"); axes[0].legend(ncol=2, fontsize=8)
axes[1].set(xlabel="Accepted physical time (s)", ylabel="(M3−M2)/M2 (%)"); axes[1].legend()
for ax in axes: ax.grid(alpha=.2)
save(fig, "figure7_damkohler_and_m3_error")

# Retain lightweight evidence and hashes. Full histories remain on HPCC.
hashes = {}
for model, case in {**cases, **{m+"_dt_half": c for m, c in half.items()}}.items():
    target = records / case["dir"].name
    target.mkdir()
    for name in ["outcome.json", "command.json", "case_spec.json", "events.csv", "rejected_trials.csv", "params.input"]:
        shutil.copy2(case["dir"]/name, target/name)
    for path in case["dir"].iterdir():
        if path.is_file(): hashes[str(path.relative_to(root))] = digest(path)
for path in (root / "unit").iterdir():
    if path.is_file(): shutil.copy2(path, records/path.name)
for path in [root/"build/soil_solver", root/"build/compiled_main.cc"]:
    hashes[str(path.relative_to(root))] = digest(path)
(out / "provenance.json").write_text(json.dumps({
    "base_commit": "c6f24b979999a9ae87602cfa1c6135597c82cbcb",
    "branch": "research/prat-analytical-growth-law-comparison",
    "remote_root": str(root),
    "postprocess_job_id": os.environ["SLURM_JOB_ID"],
    "postprocess_node": os.uname().nodename,
    "raw_history_policy": "Full accepted histories and solver logs retained on HPCC; no failed curve is extrapolated.",
    "comparison_interpolation": "Linear interpolation only within each accepted time/height range; every interpolated row is flagged.",
    "figure_policy": "General manuscript figures; target journal and final width not specified. Color is redundant with line style.",
    "sha256": hashes,
}, indent=2))
manifest = {str(path.relative_to(out)): digest(path) for path in out.rglob("*") if path.is_file() and path.name != "checksums.json"}
(out / "checksums.json").write_text(json.dumps(manifest, indent=2))
print(json.dumps(verification, indent=2))
