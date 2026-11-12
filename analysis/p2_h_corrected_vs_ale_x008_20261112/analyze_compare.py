"""Compare the corrected net-growth law with the prescribed-flux ALE run.

The comparison observable is the crystal-inventory-equivalent net height,
not the fixed-bottom ALE top coordinate.  The prescribed forcing is linearly
interpolated exactly as in the C++ driver; its piecewise-linear integral is
evaluated analytically at every ALE output time.
"""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
FIGURES = HERE / "figures"
FIGURES.mkdir(exist_ok=True)

# Frozen common parameters used by the C++ driver.
C_EQ = 0.264940869768
RHO_CR = 2165.0
EPS_BIRTH = 0.43
D_MOLECULAR = 1.3e-9
K_R = 1.0e-3
BETA = 2.0
GRAIN = 25.0e-6
H0 = 3.5e-6


def piecewise_linear_integral(t_nodes: np.ndarray, y_nodes: np.ndarray, query: np.ndarray) -> np.ndarray:
    """Integral from t_nodes[0] to query for a linearly interpolated history."""
    if np.any(np.diff(t_nodes) <= 0) or query.min() < t_nodes[0] or query.max() > t_nodes[-1]:
        raise ValueError("query must lie inside a strictly increasing forcing support")
    dt = np.diff(t_nodes)
    segment = 0.5 * (y_nodes[:-1] + y_nodes[1:]) * dt
    cumulative = np.r_[0.0, np.cumsum(segment)]
    index = np.searchsorted(t_nodes, query, side="right") - 1
    index = np.minimum(index, len(t_nodes) - 2)
    local = query - t_nodes[index]
    slope = (y_nodes[index + 1] - y_nodes[index]) / dt[index]
    return cumulative[index] + y_nodes[index] * local + 0.5 * slope * local**2


def main() -> None:
    forcing = pd.read_csv(RESULTS / "x008_forcing.csv")
    ale = pd.read_csv(RESULTS / "ale_history.csv")

    deff = np.sqrt(EPS_BIRTH) * D_MOLECULAR
    av = 36.0 * (1.0 - EPS_BIRTH) / GRAIN
    da_star = BETA**2 * (1.0 - EPS_BIRTH) ** 2 * (1.0 - C_EQ) * K_R / (
        EPS_BIRTH * deff * av
    )
    root_da_star = np.sqrt(da_star)
    phi_core = (1.0 - EPS_BIRTH) * (1.0 + root_da_star) / root_da_star

    tf = forcing["elapsed_time_s"].to_numpy()
    j = forcing["evaporation_water_kg_m2_s"].to_numpy()
    b_soil = forcing["b_soil_kg_m2_s"].to_numpy()
    p_net = (C_EQ * j - b_soil) / (1.0 - C_EQ)
    exported_soil_net = forcing["soil_net_nacl_kg_m2_s"].to_numpy()
    supply_discrepancy = p_net - exported_soil_net

    ta = ale["time_s"].to_numpy()
    net_supply = piecewise_linear_integral(tf, p_net, ta)
    corrected_height = H0 + net_supply / (RHO_CR * phi_core)
    ale_net_height = H0 + (
        ale["solid_kg_m2"].to_numpy() - ale["solid_kg_m2"].iloc[0]
    ) / (RHO_CR * phi_core)
    ale_geometric_height = ale["height_m"].to_numpy()
    difference = ale_net_height - corrected_height
    corrected_increment = corrected_height - H0

    comparison = pd.DataFrame(
        {
            "elapsed_time_s": ta,
            "elapsed_time_h": ta / 3600.0,
            "corrected_net_height_m": corrected_height,
            "ale_net_height_m": ale_net_height,
            "ale_fixed_bottom_geometric_top_m": ale_geometric_height,
            "ale_minus_corrected_m": difference,
            "ale_epsilon_bottom": ale["epsilon_bottom"].to_numpy(),
            "ale_epsilon_min": ale["epsilon_min"].to_numpy(),
            "ale_solid_kg_m2": ale["solid_kg_m2"].to_numpy(),
        }
    )
    comparison.to_csv(RESULTS / "comparison_timeseries.csv", index=False)

    dissolved = np.flatnonzero(ale["epsilon_bottom"].to_numpy() >= 1.0 - 1e-12)
    dissolution_time = float(ta[dissolved[0]]) if len(dissolved) else None
    denom = max(float(corrected_increment[-1]), np.finfo(float).tiny)
    summary = {
        "status": "completed_full_forcing_support",
        "forcing_end_s": float(tf[-1]),
        "ale_last_time_s": float(ta[-1]),
        "ale_rows": int(len(ale)),
        "da_star": float(da_star),
        "sqrt_da_star": float(root_da_star),
        "phi_core": float(phi_core),
        "initial_height_m": H0,
        "corrected_final_net_height_m": float(corrected_height[-1]),
        "ale_final_net_height_m": float(ale_net_height[-1]),
        "ale_final_geometric_top_m": float(ale_geometric_height[-1]),
        "corrected_final_increment_m": float(corrected_increment[-1]),
        "ale_final_net_increment_m": float(ale_net_height[-1] - H0),
        "final_ale_minus_corrected_m": float(difference[-1]),
        "final_relative_difference_of_increment": float(difference[-1] / denom),
        "max_abs_ale_minus_corrected_m": float(np.max(np.abs(difference))),
        "integrated_corrected_supply_kg_m2": float(net_supply[-1]),
        "integrated_exported_soil_net_kg_m2": float(
            piecewise_linear_integral(tf, exported_soil_net, np.array([tf[-1]]))[0]
        ),
        "ale_solid_gain_kg_m2": float(ale["solid_kg_m2"].iloc[-1] - ale["solid_kg_m2"].iloc[0]),
        "bottom_complete_dissolution_time_s": dissolution_time,
        "final_epsilon_bottom": float(ale["epsilon_bottom"].iloc[-1]),
        "final_epsilon_top": float(ale["epsilon_top"].iloc[-1]),
        "final_epsilon_min": float(ale["epsilon_min"].iloc[-1]),
        "max_abs_salt_residual_kg_m2": float(ale["salt_residual_kg_m2"].abs().max()),
        "max_abs_water_residual_kg_m2": float(ale["water_residual_kg_m2"].abs().max()),
        "max_abs_corrected_supply_minus_exported_soil_net_kg_m2_s": float(
            np.max(np.abs(supply_discrepancy))
        ),
        "observable_definition": "h0 + (Mcr(t)-Mcr(0))/(rho_cr*phi_core)",
        "scheduler_job": 16569032,
    }
    (RESULTS / "comparison_summary.json").write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )

    mpl.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "font.size": 10,
            "axes.labelsize": 11,
            "legend.fontsize": 9,
            "xtick.labelsize": 9,
            "ytick.labelsize": 9,
            "svg.fonttype": "none",
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
        }
    )
    colors = {"corrected": "#0072B2", "ale": "#D55E00", "geometry": "#666666"}
    fig, axes = plt.subplots(2, 1, figsize=(6.5, 5.7), sharex=True, height_ratios=(2.15, 1.0))
    time_h = ta / 3600.0
    axes[0].plot(time_h, (corrected_height - H0) * 1e6, color=colors["corrected"], lw=2.1,
                 label="Corrected law: net height")
    axes[0].plot(time_h, (ale_net_height - H0) * 1e6, color=colors["ale"], lw=1.5, ls="--",
                 label="ALE: crystal-inventory net height")
    axes[0].plot(time_h, (ale_geometric_height - H0) * 1e6, color=colors["geometry"], lw=1.25,
                 ls=":", label="ALE: fixed-bottom geometric top")
    axes[0].set_ylabel("Height increase (µm)")
    axes[0].legend(frameon=False, loc="upper left")
    axes[0].text(0.015, 0.94, "(a)", transform=axes[0].transAxes, va="top")

    axes[1].plot(time_h, difference * 1e9, color=colors["ale"], lw=1.7)
    axes[1].axhline(0.0, color="black", lw=0.7)
    axes[1].set_xlabel("Time after precipitation onset (h)")
    axes[1].set_ylabel("ALE net − corrected (nm)")
    axes[1].text(0.015, 0.92, "(b)", transform=axes[1].transAxes, va="top")
    if dissolution_time is not None:
        x = dissolution_time / 3600.0
        for axis in axes:
            axis.axvline(x, color="#009E73", lw=1.0, ls="-.")
        axes[1].annotate(
            r"bottom $\varepsilon=1$",
            xy=(x, np.interp(x, time_h, difference * 1e9)),
            xytext=(x + 4.0, 0.62 * axes[1].get_ylim()[1]),
            arrowprops={"arrowstyle": "->", "color": "#009E73", "lw": 0.8},
            color="#007A58",
        )
    for axis in axes:
        axis.grid(True, color="#D9D9D9", lw=0.6)
        axis.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(FIGURES / f"net_growth_comparison.{suffix}", dpi=300, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(6.5, 3.2))
    ax.plot(time_h, ale["epsilon_bottom"], color="#009E73", lw=1.8, label="Bottom cell")
    ax.plot(time_h, ale["epsilon_top"], color="#CC79A7", lw=1.6, ls="--", label="Top cell")
    ax.plot(time_h, ale["epsilon_min"], color="#666666", lw=1.4, ls=":", label="Minimum")
    ax.set(xlabel="Time after precipitation onset (h)", ylabel="ALE porosity, ε", ylim=(0, 1.04))
    ax.grid(True, color="#D9D9D9", lw=0.6)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(frameon=False, ncol=3, loc="lower right")
    fig.tight_layout()
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(FIGURES / f"porosity_diagnostic.{suffix}", dpi=300, bbox_inches="tight")
    plt.close(fig)

    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
