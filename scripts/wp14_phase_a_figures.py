"""Render WP14 Phase A publication figures from frozen reduced products."""
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


PAPER2 = Path(__file__).resolve().parents[1]
BASE = PAPER2 / "analysis" / "wp14_nonlinear_timestep_20261106" / "phase_a"
FIGURES = BASE / "figures"
WP12 = PAPER2 / "analysis" / "wp12_amplitude_convergence_20261106"
WP13 = PAPER2 / "analysis" / "wp13_amplitude_1e5_20261106"
GRIDS = ("R24", "R48", "R96")
AMPLITUDES = (1e-6, 1e-5, 1e-4)
GRID_COLORS = {"R24": "#5B8DB8", "R48": "#D18F4E", "R96": "#3F8C6A"}
MODE_COLORS = {"m2": "#376EA6", "m4": "#70A7C7", "m6": "#D18F4E", "m10": "#8F6A9F", "other": "#C8C8C8"}

audit_dir = Path.home() / ".agents" / "skills" / "nature-figure" / "scripts"
sys.path.insert(0, str(audit_dir))
from audit_panel_alignment import require_matplotlib_panel_alignment  # noqa: E402


mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "font.size": 7,
    "axes.labelsize": 7,
    "axes.titlesize": 7,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "legend.fontsize": 6,
    "legend.handletextpad": 2.0,
    "legend.handlelength": 1.8,
    "axes.titlepad": 8,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 0.7,
    "legend.frameon": False,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
})


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def case_dir(grid: str, amplitude: float) -> Path:
    if amplitude == 1e-6:
        return WP12 / "reference_metrics" / grid
    if amplitude == 1e-5:
        return WP13 / "metrics" / grid
    return WP12 / "enhanced_metrics" / grid


def panel_labels(axes) -> None:
    for letter, axis in zip("abcdefghijklmnopqrstuvwxyz", np.ravel(axes)):
        axis.annotate(letter, xy=(0, 1), xycoords="axes fraction", xytext=(-18, 7),
                      textcoords="offset points", fontsize=8, fontweight="bold", va="bottom")


def save_figure(fig, name: str, exclude_axes=None) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    for axis in fig.axes:
        axis.tick_params(axis="both", pad=5)
    fig.canvas.draw()
    require_matplotlib_panel_alignment(
        fig,
        json_out=str(FIGURES / f"{name}.alignment.json"),
        overlay_svg=str(FIGURES / f"{name}.alignment.svg"),
        tolerance_pt=1.5,
        gutter_tolerance_pt=1.5,
        require_panel_labels=True,
        strict=True,
        exclude_axes=exclude_axes or [],
    )
    fig.savefig(FIGURES / f"{name}.svg", bbox_inches="tight", pad_inches=0.02)
    fig.savefig(FIGURES / f"{name}.pdf", bbox_inches="tight", pad_inches=0.02)
    fig.savefig(FIGURES / f"{name}.png", dpi=300, bbox_inches="tight", pad_inches=0.02)
    fig.savefig(FIGURES / f"{name}.tiff", dpi=600, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)


def figure_1() -> None:
    series = read_csv(BASE / "normalized_harmonic_timeseries.csv")
    milestones = read_csv(BASE / "normalized_milestones.csv")
    fig = plt.figure(figsize=(180 / 25.4, 110 / 25.4), constrained_layout=True)
    grid_spec = fig.add_gridspec(3, 3, height_ratios=(1.0, 0.10, 1.0))
    axes = np.empty((2, 3), dtype=object)
    for col in range(3):
        axes[0, col] = fig.add_subplot(grid_spec[0, col], sharey=axes[0, 0] if col else None)
        axes[1, col] = fig.add_subplot(grid_spec[2, col], sharey=axes[1, 0] if col else None)
    legend_axis = fig.add_subplot(grid_spec[1, :])
    legend_axis.set_axis_off()
    for col, amplitude in enumerate(AMPLITUDES):
        axis = axes[0, col]
        for grid in GRIDS:
            rows = [r for r in series if r["grid"] == grid and float(r["amplitude"]) == amplitude]
            x = [float(r["time_over_T"]) for r in rows]
            y = [max(float(r["R_harm"]), 1e-8) for r in rows]
            axis.plot(x, y, color=GRID_COLORS[grid], lw=1.2, label=grid)
        axis.set_yscale("log", nonpositive="clip")
        axis.set_yticks((1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1, 10),
                        ("1e-5", "1e-4", "1e-3", "1e-2", "1e-1", "1", "10"))
        axis.set_title(f"Initial amplitude {amplitude:.0e}")
        axis.set_xlim(0, 6)
        axis.set_xticks((1, 2, 3, 4, 5, 6))
        axis.set_ylim(1e-5, 30)
        axis.set_xlabel("Physical time, t/T")
        if col == 0:
            axis.set_ylabel("Normalized nonseed/seed ratio, R_harm")

        axis = axes[1, col]
        x_grid = np.arange(len(GRIDS))
        bottom_colors = {0.01: "#376EA6", 0.1: "#E37B25", 1.0: "#3A9D3F"}
        for threshold, marker in zip((0.01, 0.1), ("o", "s")):
            times = []
            for grid in GRIDS:
                row = next(r for r in milestones if r["grid"] == grid
                           and float(r["amplitude"]) == amplitude and float(r["threshold"]) == threshold)
                times.append(float(row["crossing_T"]))
            axis.plot(x_grid, times, marker=marker, lw=1.0, ms=3.5, color=bottom_colors[threshold],
                      label=f"R_harm={threshold:g}")
        absolute = []
        for grid in GRIDS:
            row = next(r for r in milestones if r["grid"] == grid
                       and float(r["amplitude"]) == amplitude and float(r["threshold"]) == 0.1)
            absolute.append(float(row["absolute_H23_milestone_T"]))
        axis.plot(x_grid, absolute, color="#7A4E3A", marker="D", lw=1.0, ms=3.2, label="Absolute H23")
        axis.set_xticks(x_grid, GRIDS)
        axis.set_ylim(0.8, 4.5)
        axis.set_xlabel("Horizontal grid")
        if col == 0:
            axis.set_ylabel("Registered milestone time, t/T")
    legend_items = (
        ("Grid: R24", GRID_COLORS["R24"]), ("R48", GRID_COLORS["R48"]),
        ("R96", GRID_COLORS["R96"]), ("R_harm=0.01", bottom_colors[0.01]),
        ("R_harm=0.1", bottom_colors[0.1]), ("Absolute H23", "#7A4E3A"),
    )
    legend_positions = ((0.10, 0.72), (0.22, 0.72), (0.34, 0.72),
                        (0.51, 0.28), (0.67, 0.28), (0.85, 0.28))
    for (x, y), (label, color) in zip(legend_positions, legend_items):
        legend_axis.text(x, y, label, transform=legend_axis.transAxes, ha="center", va="center",
                         fontsize=6.2, fontweight="bold", color=color)
    panel_labels(axes)
    save_figure(fig, "Figure_1_normalized_harmonic_development", exclude_axes=[legend_axis])


def stacked_modes(axis, rows: list[dict], title: str) -> None:
    labels = [f"{float(r['amplitude']):.0e}\n{r['grid']}" for r in rows]
    x = np.arange(len(rows))
    bottom = np.zeros(len(rows))
    for key, display in (("m2_fraction", "m=2"), ("m4_fraction", "m=4"),
                         ("m6_fraction", "m=6"), ("m10_fraction", "m=10")):
        values = np.array([float(r[key]) for r in rows])
        axis.bar(x, values, bottom=bottom, width=0.78, color=MODE_COLORS[key.split("_")[0]], label=display)
        bottom += values
    other = np.maximum(0.0, 1.0 - bottom)
    axis.bar(x, other, bottom=bottom, width=0.78, color=MODE_COLORS["other"], label="Other modes")
    axis.set_xticks(x, labels)
    axis.set_ylim(0, 1)
    axis.set_ylabel("Modal energy fraction")
    axis.set_title(title)


def figure_2() -> None:
    rows = read_csv(BASE / "mode_checkpoints.csv")
    series = read_csv(BASE / "normalized_harmonic_timeseries.csv")
    fig, axes = plt.subplots(2, 2, figsize=(180 / 25.4, 118 / 25.4), constrained_layout=True)
    order = [(amplitude, grid) for amplitude in AMPLITUDES for grid in GRIDS]
    ylabels = [f"{amplitude:.0e} {grid}" for amplitude, grid in order]
    mode_palette = {2: MODE_COLORS["m2"], 4: MODE_COLORS["m4"], 6: MODE_COLORS["m6"],
                    10: MODE_COLORS["m10"]}
    axis = axes[0, 0]
    for y, (amplitude, grid) in enumerate(order):
        selected = [r for r in series if r["grid"] == grid and float(r["amplitude"]) == amplitude]
        for mode in sorted({int(r["dominant_mode"]) for r in selected}):
            points = [r for r in selected if int(r["dominant_mode"]) == mode]
            axis.scatter([float(r["time_over_T"]) for r in points], [y] * len(points), s=4,
                         color=mode_palette.get(mode, "#B8B8B8"), linewidths=0,
                         label=f"m={mode}" if y == 0 else None)
    axis.set_yticks(range(len(order)), ylabels)
    axis.set_xlim(0, 6)
    axis.set_xlabel("Physical time, t/T")
    axis.set_title("Dominant-mode evolution")

    checkpoint_order = [(amplitude, grid) for amplitude in AMPLITUDES for grid in GRIDS]
    for axis, coordinate, title in (
        (axes[0, 1], "physical_6T", "Identical physical time: 6T"),
        (axes[1, 0], "normalized_Rharm_0.1", "Identical normalized milestone: R_harm=0.1"),
        (axes[1, 1], "absolute_H23_plus_1T_auxiliary", "Auxiliary: absolute t_nl + 1T"),
    ):
        selected = []
        for amplitude, grid in checkpoint_order:
            selected.append(next(r for r in rows if r["grid"] == grid
                                 and float(r["amplitude"]) == amplitude and r["coordinate"] == coordinate))
        stacked_modes(axis, selected, title)
    handles, labels = axes[0, 1].get_legend_handles_labels()
    axes[0, 1].legend(handles, labels, loc="upper center", ncol=5, bbox_to_anchor=(0.5, 1.22))
    panel_labels(axes)
    save_figure(fig, "Figure_2_physical_and_normalized_mode_selection")


def figure_3() -> None:
    pulse = read_csv(BASE / "transport_pulse_metrics.csv")
    fig, axes = plt.subplots(3, 4, figsize=(180 / 25.4, 130 / 25.4), constrained_layout=True,
                             sharex=True, sharey=True)
    for row_index, amplitude in enumerate(AMPLITUDES):
        for col_index, depth in enumerate((1.0, 2.0, 5.0, 10.0)):
            axis = axes[row_index, col_index]
            for grid in GRIDS:
                flux = read_csv(case_dir(grid, amplitude) / "flux_series.csv")
                selected = [r for r in flux if float(r["depth_mm"]) == depth]
                t = np.array([float(r["time_over_T"]) for r in selected])
                j = np.array([float(r["advective_down_m_s"]) for r in selected])
                resolved = np.where(j > 1e-15, j * 1e11, np.nan)
                axis.plot(t, resolved, color=GRID_COLORS[grid], lw=0.8, label=grid)
                metric = next(r for r in pulse if r["grid"] == grid
                              and float(r["amplitude"]) == amplitude and float(r["depth_mm"]) == depth)
                if metric["principal_pulse_start_T"]:
                    start, end = float(metric["principal_pulse_start_T"]), float(metric["principal_pulse_end_T"])
                    mask = (t >= start) & (t <= end) & np.isfinite(resolved)
                    axis.plot(t[mask], resolved[mask], color=GRID_COLORS[grid], lw=2.0)
            axis.set_yscale("log", nonpositive="clip")
            axis.set_yticks((1e-4, 1e-3, 1e-2, 1e-1, 1, 10),
                            ("1e-4", "1e-3", "1e-2", "1e-1", "1", "10"))
            axis.set_xlim(0, 6)
            axis.set_ylim(1e-4, 30)
            if row_index == 0:
                axis.set_title(f"{depth:g} mm below surface")
            if col_index == 0:
                axis.set_ylabel(f"A={amplitude:.0e}\nResolved positive J_adv\n(x 10^-11 m s^-1)")
            if row_index == 2:
                axis.set_xlabel("Physical time, t/T")
    axes[0, 0].legend(loc="upper left")
    panel_labels(axes)
    save_figure(fig, "Figure_3_principal_positive_advective_pulses")


def figure_4() -> None:
    rows = read_csv(BASE / "grid_difference_by_amplitude.csv")
    metrics = [
        ("normalized_milestone_R0.1_T", "R_harm milestone time"),
        ("surface_mean_peak", "Surface mean peak"),
        ("surface_peak_T", "Surface peak time"),
        ("surface_decline", "Surface decline"),
        ("5mm_principal_pulse_peak_T", "5 mm main-pulse time"),
        ("5mm_peak_advective_m_s", "5 mm peak J_adv"),
        ("5mm_positive_advective_integral_m", "5 mm positive integral"),
        ("10mm_principal_pulse_peak_T", "10 mm main-pulse time"),
        ("10mm_peak_advective_m_s", "10 mm peak J_adv"),
        ("10mm_positive_advective_integral_m", "10 mm positive integral"),
        ("maximum_mirror_error", "Maximum mirror error"),
        ("m2_fraction_6T", "m=2 fraction at 6T"),
        ("m6_fraction_6T", "m=6 fraction at 6T"),
    ]
    fig, axes = plt.subplots(1, 3, figsize=(180 / 25.4, 112 / 25.4), constrained_layout=True,
                             sharey=True)
    image = None
    for axis, amplitude in zip(axes, AMPLITUDES):
        matrix = np.full((len(metrics), 2), np.nan)
        annotations = [["" for _ in range(2)] for _ in metrics]
        for i, (metric, _) in enumerate(metrics):
            row = next(r for r in rows if float(r["amplitude"]) == amplitude and r["metric"] == metric)
            for j, key in enumerate(("relative_difference_R24_R48", "relative_difference_R48_R96")):
                if row[key] != "":
                    percent = 100 * float(row[key])
                    matrix[i, j] = math.log10(1 + percent)
                    annotations[i][j] = f"{percent:.1f}%" if percent < 1000 else ">999%"
                else:
                    annotations[i][j] = "N/A"
        image = axis.imshow(matrix, aspect="auto", cmap="viridis", vmin=0, vmax=3,
                            interpolation="nearest")
        axis.set_xticks((0, 1), ("R24-R48", "R48-R96"))
        axis.set_yticks(range(len(metrics)), [label for _, label in metrics])
        axis.set_title(f"Initial amplitude {amplitude:.0e}")
        for i in range(len(metrics)):
            for j in range(2):
                value = matrix[i, j]
                color = "white" if np.isfinite(value) and value > 1.35 else "black"
                axis.text(j, i, annotations[i][j], ha="center", va="center", fontsize=5.6, color=color)
    colorbar = fig.colorbar(image, ax=axes, location="bottom", shrink=0.55, pad=0.08,
                            ticks=(0, 1, 2, 3))
    colorbar.ax.set_xticklabels(("0%", "9%", "99%", "999%"))
    colorbar.set_label("Absolute adjacent-grid difference relative to the finer grid")
    panel_labels(axes)
    save_figure(fig, "Figure_4_grid_sensitivity_by_amplitude", exclude_axes=[colorbar.ax])


def main() -> None:
    figure_1()
    figure_2()
    figure_3()
    figure_4()
    print(f"rendered 4 figures in {FIGURES}")


if __name__ == "__main__":
    main()
