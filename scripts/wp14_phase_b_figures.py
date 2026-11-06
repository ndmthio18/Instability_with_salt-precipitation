"""Render WP14 Phase B timestep-refinement figures from reduced comparison tables."""
from __future__ import annotations

import csv
import math
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np


PAPER2 = Path(__file__).resolve().parents[1]
BASE = PAPER2 / "analysis" / "wp14_nonlinear_timestep_20261106" / "phase_b"
DATA = BASE / "comparison"
FIGURES = BASE / "figures"
DT_COLOR = "#376EA6"
DT2_COLOR = "#D66B3D"
MODE_COLORS = {"m2": "#376EA6", "m4": "#70A7C7", "m6": "#D18F4E", "m10": "#8F6A9F"}

audit_dir = Path.home() / ".agents" / "skills" / "nature-figure" / "scripts"
sys.path.insert(0, str(audit_dir))
from audit_panel_alignment import require_matplotlib_panel_alignment  # noqa: E402


mpl.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans", "sans-serif"],
    "font.size": 8,
    "axes.labelsize": 8,
    "axes.titlesize": 8,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "legend.fontsize": 7,
    "legend.handletextpad": 2.5,
    "legend.handlelength": 1.8,
    "axes.titlepad": 8,
    "axes.spines.right": False,
    "axes.spines.top": False,
    "axes.linewidth": 0.7,
    "legend.frameon": False,
    "pdf.fonttype": 42,
    "svg.fonttype": "none",
})


def read_csv(name: str) -> list[dict[str, str]]:
    with (DATA / name).open(newline="") as handle:
        return list(csv.DictReader(handle))


def as_float(row: dict[str, str], key: str) -> float:
    return float(row[key])


def require_strictly_positive(values) -> np.ndarray:
    """Fail before log plotting instead of silently dropping nonpositive data."""
    array = np.asarray(values, dtype=float)
    if np.any(array <= 0):
        raise ValueError("log-scale values must be strictly positive")
    return array


def panel_labels(axes) -> None:
    for letter, axis in zip("abcdefghijklmnopqrstuvwxyz", np.ravel(axes)):
        axis.annotate(letter, xy=(0, 1), xycoords="axes fraction", xytext=(-18, 7),
                      textcoords="offset points", fontsize=9, fontweight="bold", va="bottom")


def save_figure(fig, name: str, exclude_axes=None) -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    for axis in fig.axes:
        axis.tick_params(axis="both", pad=7)
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
    fig.savefig(FIGURES / f"{name}.svg", bbox_inches="tight", pad_inches=0.08)
    fig.savefig(FIGURES / f"{name}.pdf", bbox_inches="tight", pad_inches=0.08)
    fig.savefig(FIGURES / f"{name}.png", dpi=300, bbox_inches="tight", pad_inches=0.08)
    fig.savefig(FIGURES / f"{name}.tiff", dpi=600, bbox_inches="tight", pad_inches=0.08)
    plt.close(fig)


def figure_5() -> None:
    surface = read_csv("matched_surface_and_H23.csv")
    modes = read_csv("matched_field_modal_series.csv")
    t_surface = np.array([as_float(row, "time_over_T") for row in surface])
    t_modes = np.array([as_float(row, "time_over_T") for row in modes])
    fig, axes = plt.subplots(2, 2, figsize=(180 / 25.4, 112 / 25.4), constrained_layout=True)

    axis = axes[0, 0]
    axis.plot(t_surface, [as_float(row, "surface_mean_dt") for row in surface],
              color=DT_COLOR, lw=1.1)
    axis.plot(t_surface, [as_float(row, "surface_mean_dt2") for row in surface],
              color=DT2_COLOR, lw=0.9, ls="--")
    axis.set_ylabel("Surface mean mole fraction")
    axis.set_title("Surface concentration at matched physical times")
    axis.text(0.03, 0.94, "WP12 R48, dt", transform=axis.transAxes, color=DT_COLOR, va="top")
    axis.text(0.03, 0.82, "WP14 R48, dt/2", transform=axis.transAxes, color=DT2_COLOR, va="top")

    axis = axes[0, 1]
    for key, color, style, label in (
        ("R_harm_dt", DT_COLOR, "-", "WP12 R48, dt"),
        ("R_harm_dt2", DT2_COLOR, "--", "WP14 R48, dt/2"),
    ):
        values = require_strictly_positive([max(as_float(row, key), 1e-8) for row in modes])
        axis.plot(t_modes, values, color=color, lw=1.1, ls=style)
    axis.axhline(0.1, color="#777777", lw=0.7, ls=":")
    axis.set_yscale("log")
    axis.set_ylim(1e-5, 30)
    axis.set_ylabel("Normalized nonseed/seed ratio, R_harm")
    axis.set_title("Normalized nonlinear development")
    axis.text(0.03, 0.94, "WP12 R48, dt", transform=axis.transAxes, color=DT_COLOR, va="top")
    axis.text(0.03, 0.82, "WP14 R48, dt/2", transform=axis.transAxes, color=DT2_COLOR, va="top")
    axis.text(0.12, 0.14, "Registered R_harm=0.1", transform=axis.transAxes, color="#666666")

    for axis, suffix, title in (
        (axes[1, 0], "dt", "WP12 R48 modal fractions, dt"),
        (axes[1, 1], "dt2", "WP14 R48 modal fractions, dt/2"),
    ):
        for mode, color in MODE_COLORS.items():
            axis.plot(t_modes, [as_float(row, f"{mode}_{suffix}") for row in modes],
                      color=color, lw=1.0)
        axis.set_ylim(0, 1.03)
        axis.set_ylabel("Modal energy fraction")
        axis.set_title(title, pad=24)
        for xpos, (mode, color) in zip((0.18, 0.40, 0.62, 0.84), MODE_COLORS.items()):
            axis.text(xpos, 1.015, mode.replace("m", "m="), transform=axis.transAxes,
                      color=color, ha="center", va="bottom")
    for axis in axes.flat:
        axis.set_xlim(-0.2, 6)
        axis.set_xticks(range(1, 7))
        axis.set_xlabel("Physical time, t/T")
    panel_labels(axes)
    save_figure(fig, "Figure_5_R48_dt_surface_and_modal_evolution")


def cumulative_positive(time_s: np.ndarray, flux: np.ndarray) -> np.ndarray:
    dt = np.r_[0.0, np.diff(time_s)]
    return np.cumsum(np.maximum(flux, 0.0) * dt)


def figure_6() -> None:
    rows = read_csv("matched_conservative_flux_series.csv")
    fig, axes = plt.subplots(2, 4, figsize=(180 / 25.4, 88 / 25.4), constrained_layout=False,
                             sharex=True)
    fig.subplots_adjust(left=0.09, right=0.99, bottom=0.13, top=0.91, wspace=0.40, hspace=0.40)
    for col, depth in enumerate((1.0, 2.0, 5.0, 10.0)):
        selected = [row for row in rows if as_float(row, "depth_mm") == depth]
        time_s = np.array([as_float(row, "time_s") for row in selected])
        time_t = np.array([as_float(row, "time_over_T") for row in selected])
        dt_flux = np.array([as_float(row, "advective_dt") for row in selected])
        dt2_flux = np.array([as_float(row, "advective_dt2") for row in selected])

        axis = axes[0, col]
        axis.plot(time_t, dt_flux * 1e11, color=DT_COLOR, lw=0.8)
        axis.plot(time_t, dt2_flux * 1e11, color=DT2_COLOR, lw=0.8, ls="--")
        axis.axhline(0, color="#777777", lw=0.5)
        axis.set_title(f"{depth:g} mm below surface")
        if col == 0:
            axis.set_ylabel("Advective salt flux\n(x 10^-11 m s^-1)")
            axis.text(0.04, 0.93, "WP12 dt", transform=axis.transAxes, color=DT_COLOR, va="top")
            axis.text(0.04, 0.81, "WP14 dt/2", transform=axis.transAxes, color=DT2_COLOR, va="top")

        axis = axes[1, col]
        axis.plot(time_t, cumulative_positive(time_s, dt_flux) * 1e6,
                  color=DT_COLOR, lw=1.0, label="WP12 dt")
        axis.plot(time_t, cumulative_positive(time_s, dt2_flux) * 1e6,
                  color=DT2_COLOR, lw=1.0, ls="--", label="WP14 dt/2")
        if col == 0:
            axis.set_ylabel("Cumulative positive\nadvective transport (x 10^-6 m)")
        axis.set_xlabel("Physical time, t/T")
    for axis in axes.flat:
        axis.set_xlim(0, 6)
    panel_labels(axes)
    save_figure(fig, "Figure_6_R48_dt_conservative_transport")


def relative_difference(delta: float, reference: float) -> float | None:
    scale = abs(reference)
    if not math.isfinite(delta) or not math.isfinite(scale) or scale <= 1e-30:
        return None
    return delta / scale


def metric_group(metric: str) -> str:
    if metric.startswith(("1mm", "2mm", "5mm", "10mm")):
        return "transport"
    if metric.startswith("surface"):
        return "surface"
    if metric.startswith(("m2", "m4", "m6", "m10", "odd", "maximum_mirror", "normalized", "resolved")):
        return "nonlinear/modal"
    return "other"


def figure_7() -> None:
    rows = read_csv("finite_resolution_error_budget.csv")
    colors = {"transport": "#3F8C6A", "surface": "#D18F4E", "nonlinear/modal": "#6B61A8", "other": "#777777"}
    robust = []
    for row in rows:
        if not row["delta_x"] or not row["delta_t"]:
            continue
        dx = relative_difference(as_float(row, "delta_x"), as_float(row, "R96_dt"))
        dt = relative_difference(as_float(row, "delta_t"), as_float(row, "R48_dt2"))
        if dx is None or dt is None:
            continue
        robust.append((row, max(dx, 1e-8), max(dt, 1e-8), metric_group(row["metric"])))

    fig, axes = plt.subplots(1, 2, figsize=(180 / 25.4, 88 / 25.4), constrained_layout=True)
    axis = axes[0]
    for group in ("transport", "surface", "nonlinear/modal", "other"):
        selected = [item for item in robust if item[3] == group]
        if selected:
            axis.scatter([item[1] * 100 for item in selected], [item[2] * 100 for item in selected],
                         s=17, alpha=0.8, color=colors[group], edgecolor="white", linewidth=0.3,
                         label=group.capitalize())
    limits = (1e-4, 1e4)
    require_strictly_positive(limits)
    diagonal = (2e-4, 5e3)
    axis.plot(diagonal, diagonal, color="#555555", lw=0.8, ls="--", label="Equal relative differences")
    axis.fill_between(limits, [value * 0.25 for value in limits], [value * 4 for value in limits],
                      color="#BBBBBB", alpha=0.15, edgecolor="none", linewidth=0, zorder=-1)
    axis.set_xscale("log")
    axis.set_yscale("log")
    axis.set_xlim(limits)
    axis.set_ylim(limits)
    axis.set_xlabel("Horizontal difference, |R96 dt - R48 dt| / |R96 dt| (%)")
    axis.set_ylabel("Timestep difference, |R48 dt/2 - R48 dt| / |R48 dt/2| (%)")
    axis.set_title("Finite-resolution differences across all resolved QoIs")
    axis.legend(loc="lower right")

    selected_names = [
        ("1mm_positive_integral_m", "1 mm cumulative"),
        ("surface_mean_peak", "Surface peak"),
        ("surface_peak_T", "Peak time"),
        ("2mm_principal_integral_m", "2 mm pulse"),
        ("5mm_principal_integral_m", "5 mm pulse"),
        ("10mm_principal_integral_m", "10 mm pulse"),
        ("5mm_signed_total_integral_m", "5 mm total"),
        ("m2_fraction_6T", "m=2 at 6T"),
        ("m6_fraction_6T", "m=6 at 6T"),
    ]
    axis = axes[1]
    values, labels, bar_colors = [], [], []
    lookup = {row["metric"]: row for row in rows}
    for name, label in selected_names:
        row = lookup[name]
        ratio = float(row["delta_t_over_delta_x"]) if row["delta_t_over_delta_x"] else math.nan
        values.append(max(ratio, 1e-4) if math.isfinite(ratio) else math.nan)
        labels.append(label)
        bar_colors.append(colors[metric_group(name)])
    y = np.arange(len(values))
    require_strictly_positive([value for value in values if math.isfinite(value)])
    axis.barh(y, values, color=bar_colors, height=0.72)
    axis.axvspan(0.25, 4, facecolor="#BBBBBB", alpha=0.15, edgecolor="none", linewidth=0, zorder=-1)
    axis.axvline(1, color="#555555", lw=0.8, ls="--")
    axis.set_xscale("log")
    axis.set_xlim(1e-3, 10)
    axis.set_yticks(y, labels)
    axis.invert_yaxis()
    axis.set_xlabel("Absolute difference ratio, Delta_t / Delta_x")
    axis.set_title("Selected QoIs: temporal versus horizontal variation")
    panel_labels(axes)
    save_figure(fig, "Figure_7_spatial_vs_temporal_differences")


def main() -> None:
    figure_5()
    figure_6()
    figure_7()
    print(f"rendered 3 Phase B figures in {FIGURES}")


if __name__ == "__main__":
    main()
