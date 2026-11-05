"""Create provisional, source-data-backed WP11 comparison figures.

Usage:
  python scripts/wp11_figures.py OUT LABEL=REDUCED_DIR [LABEL=REDUCED_DIR ...]
      [--dt2 REDUCED_DIR]

Every figure is exported as PNG and PDF. Lines connect accepted/saved samples
without smoothing; color and linestyle both encode the grid.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

COLORS = {"R24": "#0072B2", "R48": "#D55E00", "R96": "#009E73", "R96-dt/2": "#CC79A7"}
STYLES = {"R24": "-", "R48": "--", "R96": "-.", "R96-dt/2": ":"}
MARKERS = {"R24": "o", "R48": "s", "R96": "^", "R96-dt/2": "D"}
DEPTHS = (1.0, 2.0, 5.0, 10.0, 20.0)
MODES = (1, 2, 4, 6, 10)


def sha256(path):
    digest = hashlib.sha256(Path(path).read_bytes()).hexdigest()
    return digest


def records(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle))


def numeric(rows, key):
    return np.array([float(row[key]) for row in rows])


def finish(fig, out, stem, registry, alt):
    if fig.get_layout_engine() is None:
        fig.tight_layout()
    paths = []
    for suffix, dpi in (("png", 220), ("pdf", None)):
        path = out / f"{stem}.{suffix}"
        fig.savefig(path, dpi=dpi, bbox_inches="tight")
        paths.append({"path": path.name, "sha256": sha256(path), "bytes": path.stat().st_size})
    plt.close(fig)
    registry.append({"figure": stem, "files": paths, "description": alt})


def load_grid(label, path):
    path = Path(path)
    return {
        "label": label,
        "path": path,
        "surface": records(path / "surface_series.csv"),
        "flux": records(path / "flux_series.csv"),
        "scalar": json.loads((path / "scalar_metrics.json").read_text()),
        "modes": records(path / "tracked_mode_series.csv"),
        "field": json.loads((path / "field_metrics.json").read_text()),
        "panels": np.load(path / "field_panels.npz"),
    }


def figure_surface(grids, out, registry):
    fig, axes = plt.subplots(2, 2, figsize=(10, 6), sharex=True)
    quantities = (("surface_mean", "Surface mean"), ("surface_max", "Surface maximum"),
                  ("domain_max", "Domain maximum"), ("surface_sigma", "Surface standard deviation"))
    for axis, (key, title) in zip(axes.ravel(), quantities):
        for grid in grids:
            r = grid["surface"]
            axis.plot(numeric(r, "time_over_T"), numeric(r, key), color=COLORS[grid["label"]],
                      ls=STYLES[grid["label"]], lw=1.25, label=grid["label"])
        axis.axhline(0.1, color="0.35", ls=":", lw=0.8)
        axis.set(title=title, ylabel="salt mole fraction")
        axis.grid(alpha=0.18)
    for axis in axes[-1]:
        axis.set_xlabel("t / T")
    axes[0, 0].legend(frameon=False, ncol=len(grids))
    finish(fig, out, "Figure_01_surface_concentration", registry,
           "Four panels compare accepted-step surface mean, surface maximum, domain maximum and surface standard deviation. A dotted horizontal line marks x=0.1; curves are unsmoothed.")


def figure_flux(grids, out, registry):
    fig, axes = plt.subplots(len(DEPTHS), 2, figsize=(11, 12), sharex=True)
    for row_index, depth in enumerate(DEPTHS):
        for col_index, (key, title) in enumerate((("advective_down_m_s", "Advective"), ("total_down_m_s", "Total"))):
            axis = axes[row_index, col_index]
            for grid in grids:
                rr = [r for r in grid["flux"] if float(r["depth_mm"]) == depth]
                axis.plot(numeric(rr, "time_over_T"), numeric(rr, key), color=COLORS[grid["label"]],
                          ls=STYLES[grid["label"]], lw=1.0, label=grid["label"])
            axis.axhline(0, color="0.2", lw=0.6)
            axis.axhline(1e-15, color="0.65", ls=":", lw=0.6)
            axis.set_yscale("symlog", linthresh=1e-15)
            axis.set_ylabel(f"{depth:g} mm\nflux (m s$^{{-1}}$)")
            axis.set_title(f"{title} signed downward flux")
            axis.grid(alpha=0.16)
    for axis in axes[-1]:
        axis.set_xlabel("t / T")
    axes[0, 0].legend(frameon=False, ncol=len(grids), fontsize=8)
    finish(fig, out, "Figure_02_signed_conservative_flux", registry,
           "Rows show the 1, 2, 5, 10 and 20 mm planes; columns show signed plane-integrated advective and total salt flux. Positive is downward. The dotted line is the registered 1e-15 m/s resolvability screen, not a physical threshold.")


def panel_for_time(grid, target):
    arrays = grid["panels"]
    frames = sorted({int(key.split("_")[1]) for key in arrays.files if key.endswith("_actual_time_over_T")})
    frame = min(frames, key=lambda value: abs(float(arrays[f"frame_{value}_actual_time_over_T"]) - target))
    prefix = f"frame_{frame}"
    return {
        "frame": frame,
        "actual": float(arrays[f"{prefix}_actual_time_over_T"]),
        "px": arrays[f"{prefix}_point_x_m"],
        "pz": arrays[f"{prefix}_point_z_m"],
        "concentration": arrays[f"{prefix}_concentration"],
        "cx": arrays[f"{prefix}_cell_x_m"],
        "cz": arrays[f"{prefix}_cell_z_m"],
        "down": arrays[f"{prefix}_downward_velocity_m_s"],
    }


def figure_fields(grids, out, registry, target):
    panels = [panel_for_time(grid, target) for grid in grids]
    top_depth_mm = 20.0
    c_values, v_values = [], []
    for panel in panels:
        point_depth = (panel["pz"][-1] - panel["pz"]) * 1000
        cell_depth = (panel["pz"][-1] - panel["cz"]) * 1000
        c_values.append(panel["concentration"][point_depth <= top_depth_mm])
        v_values.append(panel["down"][cell_depth <= top_depth_mm])
    cmin, cmax = min(v.min() for v in c_values), max(v.max() for v in c_values)
    vmax = max(abs(v).max() for v in v_values)
    fig, axes = plt.subplots(2, len(grids), figsize=(4.1 * len(grids), 6.7), squeeze=False,
                             layout="constrained")
    concentration_image = velocity_image = None
    for column, (grid, panel) in enumerate(zip(grids, panels)):
        x_mm = panel["px"] * 1000
        depth_mm = (panel["pz"][-1] - panel["pz"]) * 1000
        concentration_image = axes[0, column].pcolormesh(x_mm, depth_mm, panel["concentration"],
                                                         shading="auto", cmap="viridis", vmin=cmin, vmax=cmax)
        cx_mm = panel["cx"] * 1000
        cdepth_mm = (panel["pz"][-1] - panel["cz"]) * 1000
        velocity_image = axes[1, column].pcolormesh(cx_mm, cdepth_mm, panel["down"], shading="auto",
                                                    cmap="coolwarm", vmin=-vmax, vmax=vmax)
        axes[0, column].set_title(f"{grid['label']}, {panel['actual']:.3f}T")
        for axis in axes[:, column]:
            axis.set_ylim(top_depth_mm, 0)
            axis.set_xlabel("x (mm)")
        axes[0, column].set_ylabel("depth (mm)")
        axes[1, column].set_ylabel("depth (mm)")
    fig.colorbar(concentration_image, ax=axes[0].tolist(), label="salt mole fraction", shrink=0.88)
    fig.colorbar(velocity_image, ax=axes[1].tolist(), label="qualitative downward cell velocity (m/s)", shrink=0.88)
    finish(fig, out, f"Figure_03_fields_{target:g}T", registry,
           f"Matched saved frames nearest {target:g}T show the top 20 mm. Concentration uses one common sequential scale and downward cell-centre velocity one common zero-centred diverging scale. Velocity is qualitative and is mapped through VTK cell connectivity, not interpreted as a conservative face flux.")


def figure_modes(grids, out, registry):
    fig, axes = plt.subplots(len(MODES), 1, figsize=(9, 11), sharex=True)
    for axis, mode in zip(axes, MODES):
        for grid in grids:
            rr = [r for r in grid["modes"] if int(r["mode"]) == mode]
            axis.plot(numeric(rr, "time_over_T"), numeric(rr, "energy_fraction"),
                      color=COLORS[grid["label"]], ls=STYLES[grid["label"]], lw=1.1, label=grid["label"])
        axis.set_yscale("log")
        axis.set_ylabel(f"m={mode}\nenergy fraction")
        axis.grid(alpha=0.16)
    axes[-1].set_xlabel("t / T")
    axes[0].legend(frameon=False, ncol=len(grids))
    finish(fig, out, "Figure_04_modal_energy", registry,
           "Unsmoothed Neumann-mode energy fractions for m=1, 2, 4, 6 and 10 versus physical time. Modes use k_m=m*pi/W; the same physical basis is used for every grid.")


def figure_inventory(grids, out, registry):
    fig, axes = plt.subplots(len(DEPTHS), 2, figsize=(11, 12), sharex=True)
    for row_index, depth in enumerate(DEPTHS):
        for grid in grids:
            rr = [r for r in grid["flux"] if float(r["depth_mm"]) == depth]
            time = numeric(rr, "time_over_T")
            cumulative = numeric(rr, "cumulative_advective_m") + numeric(rr, "cumulative_diffusive_m")
            inventory = numeric(rr, "inventory_salt_m")
            axes[row_index, 0].plot(time, cumulative, color=COLORS[grid["label"]],
                                    ls=STYLES[grid["label"]], lw=1.1, label=grid["label"])
            axes[row_index, 1].plot(time, inventory - inventory[0], color=COLORS[grid["label"]],
                                    ls=STYLES[grid["label"]], lw=1.1, label=grid["label"])
        axes[row_index, 0].set_ylabel(f"{depth:g} mm\ncumulative (m)")
        axes[row_index, 1].set_ylabel(f"{depth:g} mm\ninventory change (m)")
        axes[row_index, 0].grid(alpha=0.16)
        axes[row_index, 1].grid(alpha=0.16)
    axes[0, 0].set_title("Cumulative signed conservative transport")
    axes[0, 1].set_title("Upper-region salt inventory change")
    for axis in axes[-1]:
        axis.set_xlabel("t / T")
    axes[0, 0].legend(frameon=False, ncol=len(grids), fontsize=8)
    finish(fig, out, "Figure_05_cumulative_transport_inventory", registry,
           "Rows show five fixed physical depths. Left: cumulative signed advective plus diffusive transport; right: independently recorded upper-region salt inventory change from each run's initial state.")


def figure_flux_summary(grids, out, registry):
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.2), sharex=True)
    for grid in grids:
        first, peak, duration = [], [], []
        for depth in DEPTHS:
            summary = grid["scalar"]["planes"][str(depth)]["advective"]["above_1e-15"]
            first.append(np.nan if summary["first_positive_step"] is None else summary["first_positive_time_over_T"])
            peak.append(np.nan if summary["peak_positive_flux_m_s"] is None else summary["peak_positive_flux_m_s"])
            longest = summary["longest_continuous_positive_interval"]
            duration.append(np.nan if longest is None else longest["accepted_step_support_duration_over_T"])
        style = dict(color=COLORS[grid["label"]], ls=STYLES[grid["label"]], marker=MARKERS[grid["label"]], label=grid["label"])
        axes[0].plot(DEPTHS, first, **style)
        axes[1].plot(DEPTHS, peak, **style)
        axes[2].plot(DEPTHS, duration, **style)
    axes[0].set_ylabel("first accepted positive step (t/T)")
    axes[1].set_ylabel("peak positive advective flux (m/s)")
    axes[1].set_yscale("log")
    axes[2].set_ylabel("longest continuous support (T)")
    for axis in axes:
        axis.set_xlabel("depth below surface (mm)")
        axis.grid(alpha=0.2)
    axes[0].legend(frameon=False)
    finish(fig, out, "Figure_06_downward_flux_summary", registry,
           "Across grids and depths, panels report the first accepted advective flux above 1e-15 m/s, its peak, and the accepted-step support duration of the longest continuous positive interval. Missing markers mean no resolved positive interval.")


def figure_threshold(grids, out, registry):
    fig, axis = plt.subplots(figsize=(7, 4.5))
    for index, grid in enumerate(grids):
        stop = grid["scalar"]["stop"]
        color = COLORS[grid["label"]]
        if "threshold" in stop["reason"]:
            low = stop.get("global_threshold_low_s")
            high = stop.get("global_threshold_high_s")
            low = (low / (545692.2550893936 / 4.0)) if low is not None else stop["time_over_T"]
            high = (high / (545692.2550893936 / 4.0)) if high is not None else stop["time_over_T"]
            axis.plot([index, index], [low, high], color=color, lw=4)
            axis.scatter([index], [high], color=color, marker=MARKERS[grid["label"]], s=55)
        else:
            bound = stop["time_over_T"]
            axis.scatter([index], [bound], facecolors="none", edgecolors=color, marker=MARKERS[grid["label"]], s=70)
            axis.annotate("right-censored", (index, bound), xytext=(0, 8), textcoords="offset points", ha="center", fontsize=8)
    axis.axhline(25, color="0.5", ls=":", lw=0.8)
    axis.set_xticks(range(len(grids)), [g["label"] for g in grids])
    axis.set_ylabel("local x=0.1 event time or censoring bound (t/T)")
    axis.grid(axis="y", alpha=0.2)
    finish(fig, out, "Figure_07_solubility_event", registry,
           "Filled markers and vertical brackets denote the accepted-step bracket of the first local x>=0.1 event. Open markers denote right-censoring at the completed horizon; the dotted line marks 25T.")


def figure_timestep(new_grid, old_grid, out, registry):
    fig, axes = plt.subplots(2, 2, figsize=(10, 6), sharex=True)
    for grid in (new_grid, old_grid):
        rr = [r for r in grid["surface"] if float(r["time_over_T"]) <= 4.0 + 1e-8]
        style = dict(color=COLORS[grid["label"]], ls=STYLES[grid["label"]], lw=1.2, label=grid["label"])
        axes[0, 0].plot(numeric(rr, "time_over_T"), numeric(rr, "surface_mean"), **style)
        axes[0, 1].plot(numeric(rr, "time_over_T"), numeric(rr, "surface_sigma"), **style)
        flux = [r for r in grid["flux"] if float(r["depth_mm"]) == 1.0 and float(r["time_over_T"]) <= 4.0 + 1e-8]
        axes[1, 0].plot(numeric(flux, "time_over_T"), numeric(flux, "advective_down_m_s"), **style)
        axes[1, 1].plot(numeric(flux, "time_over_T"), numeric(flux, "total_down_m_s"), **style)
    axes[0, 0].set_ylabel("surface mean")
    axes[0, 1].set_ylabel("surface standard deviation")
    axes[1, 0].set_ylabel("1 mm advective flux (m/s)")
    axes[1, 1].set_ylabel("1 mm total flux (m/s)")
    for axis in axes[1]:
        axis.set_xlabel("t / T")
        axis.set_yscale("symlog", linthresh=1e-15)
        axis.axhline(0, color="0.25", lw=0.6)
    axes[0, 0].legend(frameon=False)
    for axis in axes.ravel():
        axis.grid(alpha=0.18)
    finish(fig, out, "Figure_08_R96_timestep_sensitivity", registry,
           "New full-dt R96 and independent WP9 dt/2 R96 are compared only over their common 0-4T window. Panels show surface mean, surface standard deviation and signed conservative 1 mm advective/total flux; no bitwise equality is expected.")


def main(out, specifications, dt2=None):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    grids = [load_grid(label, path) for label, path in specifications]
    plt.rcParams.update({
        "font.family": "DejaVu Sans", "font.size": 9, "axes.spines.top": False,
        "axes.spines.right": False, "pdf.fonttype": 42, "savefig.facecolor": "white",
    })
    registry = []
    figure_surface(grids, out, registry)
    figure_flux(grids, out, registry)
    for target in (5.0, 6.0):
        figure_fields(grids, out, registry, target)
    figure_modes(grids, out, registry)
    figure_inventory(grids, out, registry)
    figure_flux_summary(grids, out, registry)
    figure_threshold(grids, out, registry)
    if dt2 is not None:
        old = load_grid("R96-dt/2", dt2)
        new = next(grid for grid in grids if grid["label"] == "R96")
        figure_timestep(new, old, out, registry)
    manifest = {
        "schema_version": 1,
        "status": "provisional publication figures; target-journal compliance not claimed",
        "rendering": "accepted/saved samples connected without smoothing; color plus line style encode grid",
        "figures": registry,
        "source_directories": {grid["label"]: str(grid["path"]) for grid in grids},
    }
    (out / "figure_manifest.json").write_text(json.dumps(manifest, indent=2))
    lines = ["# Figure descriptions", "", "These are provisional publication figures; target-journal compliance is not claimed.", ""]
    for item in registry:
        lines.extend((f"## {item['figure']}", "", item["description"], ""))
    (out / "FIGURE_DESCRIPTIONS.md").write_text("\n".join(lines))
    print(json.dumps({"figures": len(registry), "out": str(out)}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("out")
    parser.add_argument("grids", nargs="+")
    parser.add_argument("--dt2")
    args = parser.parse_args()
    specs = []
    for specification in args.grids:
        label, path = specification.split("=", 1)
        specs.append((label, path))
    main(args.out, specs, args.dt2)
