"""Render the registered WP12 comparison figures from reduced source data."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

from audit_panel_alignment import require_matplotlib_panel_alignment


PAPER2 = Path(__file__).resolve().parent.parent
BASE = PAPER2 / "analysis" / "wp12_amplitude_convergence_20261106"
FIGURES = BASE / "figures"
GRIDS = ("R24", "R48", "R96")
GRID_COLORS = {"R24": "#5B8DB8", "R48": "#D18F4E", "R96": "#4C9A73"}
AMPLITUDES = (1e-6, 1e-5, 1e-4)
AMP_COLORS = {1e-6: "#7D7D7D", 1e-5: "#D18F4E", 1e-4: "#376EA6"}

mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "font.size": 7,
    "axes.labelsize": 7,
    "axes.titlesize": 7,
    "xtick.labelsize": 6,
    "ytick.labelsize": 6,
    "legend.fontsize": 6,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 0.7,
    "legend.frameon": False,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
})


def read_csv(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def load_case(grid: str, amplitude: float):
    if amplitude == 1e-6:
        path = BASE / "reference_metrics" / grid
    elif amplitude == 1e-5:
        path = PAPER2 / "analysis" / "wp13_amplitude_1e5_20261106" / "metrics" / grid
    else:
        path = BASE / "enhanced_metrics" / grid
    return {
        "path": path,
        "scalar": json.loads((path / "wp12_scalar_metrics.json").read_text()),
        "field": json.loads((path / "wp12_field_metrics.json").read_text()),
        "harmonic": read_csv(path / "harmonic_signal.csv"),
        "surface": read_csv(path / "surface_series.csv"),
        "flux": read_csv(path / "flux_series.csv"),
        "field_series": read_csv(path / "field_state_series.csv"),
    }


def panel_labels(axes):
    for letter, axis in zip("abcdefghijklmnopqrstuvwxyz", np.ravel(axes)):
        axis.annotate(letter, xy=(0, 1), xycoords="axes fraction", xytext=(-18, 7),
                      textcoords="offset points", fontsize=8, fontweight="bold", va="bottom")


def save_figure(fig, axes, name: str) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig.canvas.draw()
    require_matplotlib_panel_alignment(
        fig,
        json_out=str(FIGURES / f"{name}.alignment.json"),
        overlay_svg=str(FIGURES / f"{name}.alignment.svg"),
        tolerance_pt=1.5,
        gutter_tolerance_pt=1.5,
        require_panel_labels=True,
        strict=True,
    )
    fig.savefig(FIGURES / f"{name}.svg", bbox_inches="tight", pad_inches=0.02)
    fig.savefig(FIGURES / f"{name}.pdf", bbox_inches="tight", pad_inches=0.02)
    fig.savefig(FIGURES / f"{name}.png", dpi=300, bbox_inches="tight", pad_inches=0.02)
    fig.savefig(FIGURES / f"{name}.tiff", dpi=600, bbox_inches="tight", pad_inches=0.02)
    plt.close(fig)


def nearest(rows, target):
    return min(rows, key=lambda row: abs(float(row["time_over_T"]) - target))


def timing_figure(cases):
    fig, axes = plt.subplots(2, 3, figsize=(180 / 25.4, 105 / 25.4), constrained_layout=True, sharey=True)
    for column, grid in enumerate(GRIDS):
        for amplitude, style in ((1e-6, "--"), (1e-5, "-."), (1e-4, "-")):
            case = cases[(grid, amplitude)]
            physical = np.array([float(row["time_over_T"]) for row in case["harmonic"]])
            signal = np.maximum(np.array([float(row["H23"]) for row in case["harmonic"]]), 1e-16)
            milestone = case["scalar"]["nonlinear_milestone"]["time_over_T"]
            label = f"Amplitude {amplitude:.0e}"
            axes[0, column].plot(physical, signal, style, color=AMP_COLORS[amplitude], lw=1.2, label=label)
            axes[1, column].plot(physical - milestone, signal, style, color=AMP_COLORS[amplitude], lw=1.2, label=label)
        axes[0, column].axhline(1e-5, color="#A33F3F", lw=0.7, ls=":")
        axes[1, column].axhline(1e-5, color="#A33F3F", lw=0.7, ls=":")
        axes[0, column].set_title(grid)
        axes[0, column].set_yscale("log", nonpositive="clip")
        axes[1, column].set_yscale("log", nonpositive="clip")
        axes[0, column].tick_params(axis="y", labelsize=8)
        axes[1, column].tick_params(axis="y", labelsize=8)
        axes[0, column].set_xlim(0, 6)
        axes[1, column].set_xlim(-4, 4.1)
        axes[0, column].set_xlabel("Physical time, t/T")
        axes[1, column].set_xlabel("Aligned time, (t-t_nl)/T")
    axes[0, 0].set_ylabel("Nonlinear harmonic signal H23")
    axes[1, 0].set_ylabel("Nonlinear harmonic signal H23")
    axes[0, 2].legend(loc="lower right")
    panel_labels(axes)
    save_figure(fig, axes, "Figure_1_timing_and_registered_alignment")


def transport_figure(cases, amplitude: float, name: str):
    fig, axes = plt.subplots(2, 2, figsize=(180 / 25.4, 120 / 25.4), constrained_layout=True)
    depths = np.array([1.0, 2.0, 5.0, 10.0])
    for grid in GRIDS:
        case = cases[(grid, amplitude)]
        color = GRID_COLORS[grid]
        axes[0, 0].plot([float(row["time_over_T"]) for row in case["surface"]],
                        [float(row["surface_mean"]) for row in case["surface"]], color=color, lw=1.2, label=grid)
        flux = [row for row in case["flux"] if float(row["depth_mm"]) == 2.0]
        axes[0, 1].plot([float(row["time_over_T"]) for row in flux],
                        [float(row["advective_down_m_s"]) * 1e11 for row in flux], color=color, lw=1.0, label=grid)
        planes = case["scalar"]["planes"]
        axes[1, 0].plot(depths, [planes[str(depth)]["advective"]["peak_positive_m_s"] * 1e11 for depth in depths],
                        marker="o", ms=3, color=color, lw=1.0, label=grid)
        axes[1, 1].plot(depths, [planes[str(depth)]["advective"]["integrated_positive_m"] * 1e6 for depth in depths],
                        marker="o", ms=3, color=color, lw=1.0, label=grid)
    axes[0, 0].set(xlabel="Physical time, t/T", ylabel="Surface mean mole fraction", xlim=(0, 6))
    axes[0, 1].axhline(0, color="#444444", lw=0.6)
    axes[0, 1].set(xlabel="Physical time, t/T", ylabel="Advective salt flux at 2 mm\n(x 10^-11 m s^-1)", xlim=(0, 6))
    axes[1, 0].set(xlabel="Depth below surface (mm)", ylabel="Peak positive advective flux\n(x 10^-11 m s^-1)", xticks=depths)
    axes[1, 1].set(xlabel="Depth below surface (mm)", ylabel="Integrated positive advective flux\n(x 10^-6 m)", xticks=depths)
    axes[0, 0].legend(loc="best")
    panel_labels(axes)
    save_figure(fig, axes, name)


def mode_figure(cases):
    fig, axes = plt.subplots(2, 3, figsize=(180 / 25.4, 105 / 25.4), constrained_layout=True, sharey=True)
    modes = (2, 4, 6, 10)
    x = np.arange(len(modes))
    width = 0.24
    for column, grid in enumerate(GRIDS):
        for row_index, target_kind in enumerate(("physical", "aligned")):
            for offset, amplitude in ((-width, 1e-6), (0.0, 1e-5), (width, 1e-4)):
                case = cases[(grid, amplitude)]
                target = 6.0 if target_kind == "physical" else case["scalar"]["nonlinear_milestone"]["time_over_T"] + 1.0
                state = nearest(case["field_series"], target)
                fractions = [float(state[f"mode_{mode}_energy_fraction"]) for mode in modes]
                axes[row_index, column].bar(x + offset, fractions, width=width, color=AMP_COLORS[amplitude],
                                            label=f"Amplitude {amplitude:.0e}")
            axes[row_index, column].set_xticks(x, [str(mode) for mode in modes])
            axes[row_index, column].set_ylim(0, 1)
            axes[row_index, column].set_xlabel("Horizontal mode m")
        axes[0, column].set_title(grid)
    axes[0, 0].set_ylabel("Energy fraction at 6T")
    axes[1, 0].set_ylabel("Energy fraction at t_nl+1T")
    handles, labels = axes[0, 2].get_legend_handles_labels()
    fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 1.035), ncol=3)
    panel_labels(axes)
    save_figure(fig, axes, "Figure_3_modal_energy_physical_and_aligned")


def symmetry_timing_figure(cases):
    fig, axes = plt.subplots(2, 2, figsize=(180 / 25.4, 95 / 25.4), constrained_layout=True)
    x = np.arange(3)
    width = 0.24
    for offset, amplitude in ((-width, 1e-6), (0.0, 1e-5), (width, 1e-4)):
        odd = [max(cases[(grid, amplitude)]["field"]["maximum_odd_mode_energy_fraction"], 1e-16) for grid in GRIDS]
        mirror = [max(cases[(grid, amplitude)]["field"]["maximum_mirror_symmetry_error"], 1e-16) for grid in GRIDS]
        axes[0, 0].bar(x + offset, odd, width=width, color=AMP_COLORS[amplitude], label=f"Amplitude {amplitude:.0e}")
        axes[0, 1].bar(x + offset, mirror, width=width, color=AMP_COLORS[amplitude], label=f"Amplitude {amplitude:.0e}")
    for axis in axes[0]:
        axis.set_xticks(x, GRIDS)
        axis.set_yscale("log", nonpositive="clip")
        axis.tick_params(axis="y", labelsize=8)
    axes[0, 0].set_ylabel("Maximum odd-mode energy fraction")
    axes[0, 1].set_ylabel("Maximum mirror-symmetry error")
    axes[0, 1].legend(loc="best")
    for axis, metric, ylabel in ((axes[1, 0], "nonlinear_milestone", "Nonlinear milestone, t/T"),
                                 (axes[1, 1], "downflow", "Resolved downflow onset, t/T")):
        for grid_index, grid in enumerate(GRIDS):
            values = []
            for amplitude in AMPLITUDES:
                case = cases[(grid, amplitude)]
                if metric == "nonlinear_milestone":
                    values.append(case["scalar"]["nonlinear_milestone"]["time_over_T"])
                else:
                    values.append(case["field"]["resolved_downflow_onset"]["time_over_T"])
            positions = [grid_index - 0.16, grid_index, grid_index + 0.16]
            axis.plot(positions, values, color=GRID_COLORS[grid], lw=0.9)
            axis.scatter(positions, values, color=[AMP_COLORS[amplitude] for amplitude in AMPLITUDES], s=13, zorder=3)
        axis.set_xticks(x, GRIDS)
        axis.set_ylabel(ylabel)
        axis.set_ylim(0, 6)
    panel_labels(axes)
    save_figure(fig, axes, "Figure_4_symmetry_and_timing")


def main() -> None:
    cases = {(grid, amplitude): load_case(grid, amplitude) for grid in GRIDS for amplitude in AMPLITUDES}
    timing_figure(cases)
    transport_figure(cases, 1e-4, "Figure_2_A1e-4_grid_transport")
    mode_figure(cases)
    symmetry_timing_figure(cases)
    transport_figure(cases, 1e-5, "Figure_5_A1e-5_grid_transport")
    print(FIGURES)


if __name__ == "__main__":
    main()
