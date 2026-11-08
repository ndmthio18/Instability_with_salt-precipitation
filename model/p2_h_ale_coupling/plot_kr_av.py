"""Publication-style P2-H k_r/a_v comparison figures."""
import csv
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

root = pathlib.Path.cwd()
production = root / "production_corrected"
analysis = root / "analysis"
figures = analysis / "figures"
results = analysis / "results"
figures.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans",
    "font.size": 8.5,
    "axes.labelsize": 9,
    "axes.titlesize": 9,
    "legend.fontsize": 8,
    "axes.linewidth": 0.8,
    "lines.linewidth": 1.6,
    "xtick.major.width": 0.8,
    "ytick.major.width": 0.8,
})

groups = {
    "kr": [
        ("kr0p1_av1", r"$0.1k_{r,\mathrm{ref}}$"),
        ("kr1_av1_ref", r"$k_{r,\mathrm{ref}}$"),
        ("kr10_av1", r"$10k_{r,\mathrm{ref}}$"),
    ],
    "av": [
        ("kr1_av0p5", r"$0.5a_{v,\mathrm{ref}}$"),
        ("kr1_av1_ref", r"$a_{v,\mathrm{ref}}$"),
        ("kr1_av2", r"$2a_{v,\mathrm{ref}}$"),
    ],
}
colors = ["#0072B2", "#009E73", "#D55E00"]

def rows(path):
    with path.open(newline="") as handle:
        return [{key: float(value) for key, value in row.items()} for row in csv.DictReader(handle)]

def grouped(path):
    data = {}
    with path.open(newline="") as handle:
        for raw in csv.DictReader(handle):
            time = float(raw["time_s"])
            data.setdefault(time, []).append({key: float(value) for key, value in raw.items()})
    return data

summary = {row["case"]: row for row in csv.DictReader((results / "summary.csv").open(newline=""))}
all_completed = all(int(row["completed_50h"]) == 1 for row in summary.values())
maximum_end_h = max(float(row["actual_end_s"]) for row in summary.values()) / 3600
plot_end_h = 50.0 if all_completed else max(0.5, np.ceil(maximum_end_h * 2.1) / 2)

def save(fig, stem):
    fig.savefig(figures / f"{stem}.png", dpi=300, bbox_inches="tight")
    fig.savefig(figures / f"{stem}.pdf", bbox_inches="tight")
    fig.savefig(figures / f"{stem}.svg", bbox_inches="tight")
    plt.close(fig)

def soil_cloud(case, depth_grid):
    profiles = grouped(production / case / "soil_fields.csv")
    times = np.array(sorted(profiles))
    values = np.full((len(depth_grid), len(times)), np.nan)
    for column, time in enumerate(times):
        profile = profiles[time]
        depth = np.array([0.2 - row["y_m"] for row in profile])
        concentration = np.array([row["mass_fraction"] for row in profile])
        order = np.argsort(depth)
        values[:, column] = np.interp(depth_grid, depth[order], concentration[order])
    return times, values

def crust_cloud(case, height_grid):
    profiles = grouped(production / case / "crust_fields.csv")
    if not profiles:
        end_time = float(summary[case]["actual_end_s"])
        return np.array([0.0, end_time]), np.full((len(height_grid), 2), np.nan)
    times = np.array(sorted(profiles))
    values = np.full((len(height_grid), len(times)), np.nan)
    for column, time in enumerate(times):
        profile = profiles[time]
        center = np.array([0.5 * (row["z_bottom_m"] + row["z_top_m"]) for row in profile])
        porosity = np.array([row["epsilon"] for row in profile])
        height = max(row["z_top_m"] for row in profile)
        inside = height_grid <= height
        values[inside, column] = np.interp(height_grid[inside], center, porosity, left=porosity[0], right=porosity[-1])
    return times, values

all_cases = sorted({name for group in groups.values() for name, _ in group})
soil_raw = {name: soil_cloud(name, np.linspace(0, 0.01, 220)) for name in all_cases}
maximum_height = max(float(summary[name]["final_crust_height_m"]) for name in all_cases)
height_grid = np.linspace(0, max(maximum_height, 1e-7), 240)
crust_raw = {name: crust_cloud(name, height_grid) for name in all_cases}
soil_min = min(np.nanmin(value) for _, value in soil_raw.values())
soil_max = max(np.nanmax(value) for _, value in soil_raw.values())

def evolution_figure(group_key, stem):
    fig, axes = plt.subplots(2, 3, figsize=(11.3, 5.9), sharex=True, constrained_layout=True)
    crust_image = plt.cm.ScalarMappable(norm=plt.Normalize(0.02, 0.43), cmap="viridis_r")
    soil_image = None
    for column, (name, label) in enumerate(groups[group_key]):
        crust_time, crust_value = crust_raw[name]
        soil_time, soil_value = soil_raw[name]
        if not np.isfinite(crust_value).any():
            axes[0, column].text(0.5, 0.5, "No ALE activation", transform=axes[0, column].transAxes,
                                 ha="center", va="center", color="0.35")
        elif crust_time.size == 1:
            valid = np.isfinite(crust_value[:, 0])
            axes[0, column].scatter(
                np.full(np.count_nonzero(valid), crust_time[0] / 3600),
                height_grid[valid] * 1e6, c=crust_value[valid, 0], s=4,
                cmap="viridis_r", vmin=0.02, vmax=0.43, linewidths=0,
                rasterized=True,
            )
            axes[0, column].text(0.98, 0.96, "single accepted\nactivation state",
                                 transform=axes[0, column].transAxes,
                                 ha="right", va="top", color="0.25", fontsize=7)
        else:
            axes[0, column].pcolormesh(
                crust_time / 3600, height_grid * 1e6, crust_value,
                shading="auto", cmap="viridis_r", vmin=0.02, vmax=0.43,
                rasterized=True,
            )
        soil_image = axes[1, column].pcolormesh(
            soil_time / 3600, np.linspace(0, 10, soil_value.shape[0]), soil_value,
            shading="auto", cmap="magma", vmin=soil_min, vmax=soil_max,
            rasterized=True,
        )
        end_h = float(summary[name]["actual_end_s"]) / 3600
        if int(summary[name]["completed_50h"]) == 0:
            for row_index, axis in enumerate(axes[:, column]):
                marker_color = "0.2" if row_index == 0 else "white"
                axis.axvline(end_h, color=marker_color, ls="--", lw=1.0)
                axis.text(end_h, 0.98, "stop", transform=axis.get_xaxis_transform(),
                          ha="right", va="top", color=marker_color, fontsize=7)
        axes[0, column].set_title(label)
        axes[0, column].set_xlim(0, plot_end_h)
        axes[0, column].set_ylim(0, height_grid[-1] * 1e6)
        axes[1, column].set_xlabel("Time (h)")
        axes[1, column].invert_yaxis()
    axes[0, 0].set_ylabel("Height above soil (µm)")
    axes[1, 0].set_ylabel("Depth below surface (mm)")
    for index, axis in enumerate(axes.flat):
        label_color = "0.15" if index < 3 else "white"
        axis.text(0.02, 0.96, f"({chr(97 + index)})", transform=axis.transAxes,
                  ha="left", va="top", fontweight="bold", color=label_color)
    fig.colorbar(crust_image, ax=axes[0, :], label="Crust porosity, ε", pad=0.015)
    fig.colorbar(soil_image, ax=axes[1, :], label="Soil dissolved NaCl mass fraction", pad=0.015)
    save(fig, stem)

evolution_figure("kr", "figure1_fixed_av_kr_evolution")
evolution_figure("av", "figure2_fixed_kr_av_evolution")

fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.6), sharey=True, constrained_layout=True)
all_rates = []
for axis, group_key, title in zip(
    axes, ("kr", "av"), (r"Fixed $a_v=a_{v,\mathrm{ref}}$", r"Fixed $k_r=k_{r,\mathrm{ref}}$")
):
    for color, (name, label) in zip(colors, groups[group_key]):
        history = rows(results / "time_series" / f"{name}.csv")
        all_rates.extend(row["evaporation_kg_m2_s"] for row in history)
        axis.plot([row["time_s"] / 3600 for row in history],
                  [row["evaporation_kg_m2_s"] for row in history], color=color, label=label)
        if int(summary[name]["completed_50h"]) == 0:
            axis.plot(history[-1]["time_s"] / 3600, history[-1]["evaporation_kg_m2_s"],
                      marker="x", ms=6, mew=1.2, color=color)
    axis.set(xlabel="Time (h)", title=title, xlim=(0, plot_end_h))
    axis.legend(frameon=False)
    axis.text(0.02, 0.96, f"({chr(97 + list(axes).index(axis))})", transform=axis.transAxes,
              ha="left", va="top", fontweight="bold")
axes[0].set_ylabel(r"Evaporation rate (kg m$^{-2}$ s$^{-1}$)")
rate_pad = max(0.1 * (max(all_rates) - min(all_rates)), 0.005 * np.mean(all_rates))
axes[0].set_ylim(min(all_rates) - rate_pad, max(all_rates) + rate_pad)
save(fig, "figure3_evaporation_rate")

fig, axes = plt.subplots(2, 2, figsize=(9.0, 6.3), sharex="col", constrained_layout=True)
for column, group_key in enumerate(("kr", "av")):
    for color, (name, label) in zip(colors, groups[group_key]):
        history = rows(results / "time_series" / f"{name}.csv")
        time = [row["time_s"] / 3600 for row in history]
        axes[0, column].plot(time, [row["h_mass_m"] * 1e6 for row in history],
                             color=color, ls="--", label=label)
        axes[0, column].plot(time, [row["h_m"] * 1e6 for row in history], color=color)
        active = [row for row in history if row["active"] == 1]
        if active:
            axes[1, column].plot([row["time_s"] / 3600 for row in active],
                                 [row["epsilon_min"] for row in active], color=color)
    axes[0, column].legend(frameon=False)
    axes[0, column].set_title(
        r"Fixed $a_v=a_{v,\mathrm{ref}}$" if group_key == "kr"
        else r"Fixed $k_r=k_{r,\mathrm{ref}}$"
    )
    axes[1, column].set_xlabel("Time (h)")
    axes[0, column].set_xlim(0, plot_end_h)
    axes[0, column].text(0.98, 0.96, "solid: ALE height\ndashed: solid-equivalent height",
                         transform=axes[0, column].transAxes, ha="right", va="top", fontsize=7)
    if not any(int(summary[name]["ale_activated"]) for name, _ in groups[group_key]):
        axes[1, column].text(0.5, 0.5, "No accepted ALE state",
                             transform=axes[1, column].transAxes,
                             ha="center", va="center", color="0.35")
axes[0, 0].set_ylabel("Height metric (µm)")
axes[1, 0].set_ylabel("Minimum crust porosity")
axes[1, 0].set_ylim(0.015, 0.44)
axes[1, 1].set_ylim(0.015, 0.44)
for index, axis in enumerate(axes.flat):
    axis.text(0.02, 0.96, f"({chr(97 + index)})", transform=axis.transAxes,
              ha="left", va="top", fontweight="bold")
save(fig, "figure4_height_and_minimum_porosity")
