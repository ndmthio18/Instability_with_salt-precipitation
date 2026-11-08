"""Publication-style reservoir-fed P2-H sensitivity figures."""
import csv
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

root = pathlib.Path.cwd()
raw = root / "matrix"
analysis = root / "analysis"
results = analysis / "results"
figures = analysis / "figures"
figures.mkdir(parents=True, exist_ok=True)

plt.rcParams.update({
    "font.family": "DejaVu Sans", "font.size": 8.5, "axes.labelsize": 9,
    "axes.titlesize": 9, "legend.fontsize": 7.5, "axes.linewidth": 0.8,
    "lines.linewidth": 1.5, "xtick.major.width": 0.8, "ytick.major.width": 0.8,
})
colors = ("#0072B2", "#009E73", "#D55E00")
styles = ("--", "-", ":")
groups = {
    "kr": (("kr0p1_av1", r"$0.1k_{r,ref}$"), ("kr1_av1_ref", r"$k_{r,ref}$"), ("kr10_av1", r"$10k_{r,ref}$")),
    "av": (("kr1_av0p5", r"$0.5a_{v,ref}$"), ("kr1_av1_ref", r"$a_{v,ref}$"), ("kr1_av2", r"$2a_{v,ref}$")),
}


def read(path):
    with path.open(newline="") as handle:
        return [{key: float(value) for key, value in row.items()} for row in csv.DictReader(handle)]


def group_rows(path):
    grouped = {}
    for row in read(path):
        grouped.setdefault(row["time_s"], []).append(row)
    return grouped


def save(fig, name):
    for ext in ("png", "pdf", "svg"):
        fig.savefig(figures / f"{name}.{ext}", dpi=300 if ext == "png" else None, bbox_inches="tight")
    plt.close(fig)


summary = {row["case"]: row for row in csv.DictReader((results / "summary.csv").open(newline=""))}
histories = {name: read(results / "time_series" / f"{name}.csv") for name in summary}
boundaries = {name: read(results / "time_series" / f"{name}_boundary.csv") for name in summary}
max_height = max(float(row["final_height_m"]) for row in summary.values())
z_grid = np.linspace(0, max_height, 260)


def crust_cloud(name):
    grouped = group_rows(raw / name / "crust_fields.csv")
    times = np.array(sorted(grouped))
    cloud = np.full((len(z_grid), len(times)), np.nan)
    for column, time in enumerate(times):
        profile = grouped[time]
        center = np.array([0.5 * (row["z_bottom_m"] + row["z_top_m"]) for row in profile])
        epsilon = np.array([row["epsilon"] for row in profile])
        height = max(row["z_top_m"] for row in profile)
        inside = z_grid <= height
        cloud[inside, column] = np.interp(z_grid[inside], center, epsilon, left=epsilon[0], right=epsilon[-1])
    return times, cloud


clouds = {name: crust_cloud(name) for name in summary}


soil_depth = np.linspace(0, 0.01, 240)
def soil_cloud(name):
    grouped = group_rows(raw / name / "soil_fields.csv")
    times = np.array(sorted(grouped))
    cloud = np.full((len(soil_depth), len(times)), np.nan)
    for column, time in enumerate(times):
        profile = grouped[time]
        depth = np.array([0.2 - row["y_m"] for row in profile])
        concentration = np.array([row["mass_fraction"] for row in profile])
        order = np.argsort(depth)
        cloud[:, column] = np.interp(soil_depth, depth[order], concentration[order])
    return times, cloud


soil_clouds = {name: soil_cloud(name) for name in summary}
soil_min = min(np.nanmin(cloud) for _, cloud in soil_clouds.values())
soil_max = max(np.nanmax(cloud) for _, cloud in soil_clouds.values())


def cloud_figure(key, filename, title):
    fig, axes = plt.subplots(2, 3, figsize=(10.8, 6.0), sharex=True, constrained_layout=True)
    crust_image = soil_image = None
    for index, (name, label) in enumerate(groups[key]):
        axis = axes[0, index]
        times, cloud = clouds[name]
        crust_image = axis.pcolormesh(times / 3600, z_grid * 1e6, cloud, shading="auto", cmap="viridis_r", vmin=0.02, vmax=0.43, rasterized=True)
        activation = float(summary[name]["activation_s"]) / 3600
        axis.axvline(activation, color="white", lw=0.9, ls="--")
        axis.set(title=label, xlim=(0, 50), ylim=(0, max_height * 1e6))
        axis.text(0.02, 0.97, f"({chr(97+index)})", transform=axis.transAxes, ha="left", va="top", fontweight="bold", color="white")
        soil_time, soil = soil_clouds[name]
        lower = axes[1, index]
        soil_image = lower.pcolormesh(soil_time / 3600, soil_depth * 1e3, soil, shading="auto", cmap="magma", vmin=soil_min, vmax=soil_max, rasterized=True)
        lower.axvline(activation, color="white", lw=0.9, ls="--")
        lower.set(xlabel="Time (h)", xlim=(0, 50), ylim=(0, 10))
        lower.invert_yaxis()
        lower.text(0.02, 0.97, f"({chr(100+index)})", transform=lower.transAxes, ha="left", va="top", fontweight="bold", color="white")
        end_h = float(summary[name]["actual_end_s"]) / 3600
        if int(summary[name]["completed_50h"]) == 0:
            for target_axis in (axis, lower):
                target_axis.axvline(end_h, color="white", lw=1.1, ls=":")
                target_axis.text(end_h, 0.03, "censored", transform=target_axis.get_xaxis_transform(), rotation=90, ha="right", va="bottom", color="white", fontsize=6.5)
    axes[0, 0].set_ylabel("Actual crust height above soil (µm)")
    axes[1, 0].set_ylabel("Depth below surface (mm)")
    fig.suptitle(title, fontsize=10)
    fig.colorbar(crust_image, ax=axes[0, :], label="Local ALE porosity, ε", pad=0.015)
    fig.colorbar(soil_image, ax=axes[1, :], label="Dissolved NaCl mass fraction", pad=0.015)
    save(fig, filename)


cloud_figure("kr", "figure1_porosity_cloud_kr", r"Fixed $a_v=a_{v,ref}$")
cloud_figure("av", "figure2_porosity_cloud_av", r"Fixed $k_r=k_{r,ref}$")

fig, axes = plt.subplots(1, 2, figsize=(9.0, 3.45), sharex=True, sharey=True, constrained_layout=True)
for column, key in enumerate(("kr", "av")):
    for color, style, (name, label) in zip(colors, styles, groups[key]):
        h = histories[name]
        axes[column].plot([row["time_s"] / 3600 for row in h], [row["evaporation_kg_m2_s"] for row in h], color=color, ls=style, label=label)
        if int(summary[name]["completed_50h"]) == 0:
            axes[column].plot(h[-1]["time_s"] / 3600, h[-1]["evaporation_kg_m2_s"], marker="x", color=color, ms=6)
    axes[column].set(xlabel="Time (h)", xlim=(0, 50), title=(r"Fixed $a_v$" if key == "kr" else r"Fixed $k_r$"))
    axes[column].legend(frameon=False)
    axes[column].text(0.02, 0.97, f"({chr(97+column)})", transform=axes[column].transAxes, ha="left", va="top", fontweight="bold")
axes[0].set_ylabel(r"Top film evaporation (kg m$^{-2}$ s$^{-1}$)")
save(fig, "figure3_evaporation")

fig, axes = plt.subplots(3, 2, figsize=(9.0, 8.0), sharex="col", constrained_layout=True)
for column, key in enumerate(("kr", "av")):
    for color, style, (name, label) in zip(colors, styles, groups[key]):
        h = histories[name]
        time = [row["time_s"] / 3600 for row in h]
        axes[0, column].plot(time, [row["h_m"] * 1e6 for row in h], color=color, ls=style, label=label)
        axes[1, column].plot(time, [row["epsilon_min"] for row in h], color=color, ls=style)
        axes[2, column].plot(time, [row["internal_fraction_net"] for row in h], color=color, ls=style)
        if int(summary[name]["completed_50h"]) == 0:
            for row_index, metric in enumerate(("h_m", "epsilon_min", "internal_fraction_net")):
                value = h[-1][metric] * (1e6 if metric == "h_m" else 1)
                axes[row_index, column].plot(time[-1], value, marker="x", color=color, ms=6)
    axes[0, column].set_title(r"Fixed $a_v$" if key == "kr" else r"Fixed $k_r$")
    axes[0, column].legend(frameon=False)
    axes[2, column].set(xlabel="Time (h)", xlim=(0, 50))
axes[0, 0].set_ylabel("Actual ALE height (µm)")
axes[1, 0].set_ylabel("Minimum porosity")
axes[2, 0].set_ylabel("Internal solid fraction")
for index, axis in enumerate(axes.flat):
    axis.text(0.02, 0.97, f"({chr(97+index)})", transform=axis.transAxes, ha="left", va="top", fontweight="bold")
save(fig, "figure4_height_porosity_internal_fraction")

fig, axes = plt.subplots(4, 2, figsize=(9.0, 10.0), sharex="col", constrained_layout=True)
for column, key in enumerate(("kr", "av")):
    for color, style, (name, label) in zip(colors, styles, groups[key]):
        h = histories[name]; b = boundaries[name]
        ht = [row["time_s"] / 3600 for row in h]; bt = [row["time_s"] / 3600 for row in b]
        axes[0, column].plot(bt, [row["bottom_water_in_kg_m2_s"] for row in b], color=color, ls=style, label=label)
        axes[1, column].plot(ht, [row["evaporation_kg_m2_s"] for row in h], color=color, ls=style)
        axes[2, column].plot(bt, [row["cumulative_bottom_water_kg"] / 0.001 for row in b], color=color, ls=style)
        axes[3, column].plot(ht, [row["top_soil_sw"] for row in h], color=color, ls=style)
        if int(summary[name]["completed_50h"]) == 0:
            endpoint = h[-1]["time_s"] / 3600
            for row_index, value in enumerate((
                b[-1]["bottom_water_in_kg_m2_s"], h[-1]["evaporation_kg_m2_s"], b[-1]["cumulative_bottom_water_kg"] / 0.001, h[-1]["top_soil_sw"]
            )):
                axes[row_index, column].plot(endpoint, value, marker="x", color=color, ms=6)
    axes[0, column].set_title(r"Fixed $a_v$" if key == "kr" else r"Fixed $k_r$")
    axes[0, column].legend(frameon=False)
    axes[3, column].set(xlabel="Time (h)", xlim=(0, 50))
axes[0, 0].set_ylabel(r"Bottom water influx (kg m$^{-2}$ s$^{-1}$)")
axes[1, 0].set_ylabel(r"Top evaporation (kg m$^{-2}$ s$^{-1}$)")
axes[2, 0].set_ylabel(r"Cumulative bottom water (kg m$^{-2}$)")
axes[3, 0].set_ylabel(r"Surface soil $S_w$")
for index, axis in enumerate(axes.flat):
    axis.text(0.02, 0.97, f"({chr(97+index)})", transform=axis.transAxes, ha="left", va="top", fontweight="bold")
save(fig, "figure5_reservoir_supply")
