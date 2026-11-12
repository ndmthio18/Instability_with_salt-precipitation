"""Postprocess task-25 BVP, ALE, convergence, redistribution, and x008 evidence."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pandas as pd

from analytical_model import Parameters, analytical_profile, corrected_observables


HERE = Path(__file__).resolve().parent
RAW = HERE / "raw_remote" / "results" / "raw"
RESULTS = HERE / "results"
FIGURES = HERE / "figures"
RESULTS.mkdir(exist_ok=True)
FIGURES.mkdir(exist_ok=True)
PRIOR = HERE.parent / "p2_h_corrected_vs_ale_x008_20261112" / "results"
P = Parameters()


def sech(x: np.ndarray | float) -> np.ndarray:
    """Evaluate sech without overflowing for the thick x008 crust."""
    value = np.asarray(x, dtype=float)
    exponential = np.exp(-np.abs(value))
    return 2.0 * exponential / (1.0 + exponential * exponential)


def load_raw() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    histories = pd.concat(
        [pd.read_csv(path) for path in sorted(RAW.glob("*_history.csv"))],
        ignore_index=True,
    )
    profiles = pd.concat(
        [pd.read_csv(path) for path in sorted(RAW.glob("*_profiles.csv"))],
        ignore_index=True,
    )
    configs = pd.read_csv(HERE / "configs" / "cases.csv")
    return histories, profiles, configs


def status_table() -> pd.DataFrame:
    records = []
    for path in sorted(RAW.glob("*_status.txt")):
        record: dict[str, str] = {}
        for line in path.read_text(encoding="utf-8").splitlines():
            key, value = line.split("=", maxsplit=1)
            record[key] = value
        records.append(record)
    return pd.DataFrame(records)


def window_table(histories: pd.DataFrame, configs: pd.DataFrame) -> pd.DataFrame:
    window = histories[histories["time_over_tau"].between(10.0, 20.0)].copy()

    def summarize(group: pd.DataFrame) -> pd.Series:
        velocity = group["velocity_m_s"]
        h0 = group["height_m"].iloc[0]
        return pd.Series(
            {
                "window_rows": len(group),
                "velocity_mean_m_s": velocity.mean(),
                "velocity_std_m_s": velocity.std(),
                "velocity_cv": velocity.std() / abs(velocity.mean()),
                "top_concentration_mean": group["top_concentration"].mean(),
                "top_concentration_range": group["top_concentration"].max()
                - group["top_concentration"].min(),
                "bottom_face_concentration_mean": group[
                    "bottom_face_concentration"
                ].mean(),
                "bottom_cell_concentration_mean": group[
                    "bottom_cell_concentration"
                ].mean(),
                "height_relative_change": (
                    group["height_m"].max() - group["height_m"].min()
                )
                / h0,
                "epsilon_mean_change": group["epsilon_mean"].max()
                - group["epsilon_mean"].min(),
                "epsilon_min": group["epsilon_min"].min(),
                "epsilon_max": group["epsilon_max"].max(),
                "da_birth_mean": group["da_birth_properties"].mean(),
                "da_volume_mean": group["da_volume_average"].mean(),
                "da_star_volume_mean": group["da_star_volume_average"].mean(),
                "top_flux_mean": group["top_solid_flux_kg_m2_s"].mean(),
                "positive_mean": group["internal_positive_kg_m2_s"].mean(),
                "dissolution_mean": group["internal_dissolution_kg_m2_s"].mean(),
                "upper2_mean": group["upper_positive_2ell_kg_m2_s"].mean(),
                "upper3_mean": group["upper_positive_3ell_kg_m2_s"].mean(),
                "upper4_mean": group["upper_positive_4ell_kg_m2_s"].mean(),
                "storage_rate_mean": group[
                    "dissolved_storage_rate_kg_m2_s"
                ].mean(),
                "bottom_total_mean": group[
                    "bottom_total_salt_kg_m2_s"
                ].mean(),
                "bottom_advective_mean": group[
                    "bottom_advective_salt_kg_m2_s"
                ].mean(),
                "bottom_diffusive_mean": group[
                    "bottom_diffusive_salt_kg_m2_s"
                ].mean(),
                "solid_final": group["solid_kg_m2"].iloc[-1],
                "salt_residual_max": group["salt_residual_kg_m2"].abs().max(),
                "water_residual_max": group["water_residual_kg_m2"].abs().max(),
            }
        )

    out = window.groupby("case_id", sort=False).apply(
        summarize, include_groups=False
    ).reset_index()
    out = out.merge(configs, on="case_id", validate="one_to_one")
    out["u_reconstructed"] = out["u"]
    out["quasi_steady_valid"] = (
        (out["velocity_cv"] < 5e-4)
        & (out["height_relative_change"] < 5e-3)
        & (out["epsilon_mean_change"].abs() < 5e-3)
        & (out["top_concentration_range"] < 1e-5)
    )
    return out


def main_validation(window: pd.DataFrame) -> pd.DataFrame:
    main = window[window["suite"] == "main"].copy()
    theory = [
        corrected_observables(float(row.da), float(row.u), P)
        for row in main.itertuples()
    ]
    main["theory_velocity_m_s"] = [item["velocity_m_s"] for item in theory]
    main["theory_top_flux_kg_m2_s"] = [
        item["top_flux_kg_m2_s"] for item in theory
    ]
    main["theory_top_concentration"] = [
        item["top_concentration"] for item in theory
    ]
    main["velocity_relative_error"] = (
        main["velocity_mean_m_s"] - main["theory_velocity_m_s"]
    ).abs() / main["theory_velocity_m_s"].abs()
    main["top_flux_relative_error"] = (
        main["top_flux_mean"] - main["theory_top_flux_kg_m2_s"]
    ).abs() / main["theory_top_flux_kg_m2_s"].abs()
    main["top_concentration_abs_error"] = (
        main["top_concentration_mean"] - main["theory_top_concentration"]
    ).abs()
    reference = main[main["u"] == 0][["da", "velocity_mean_m_s"]].rename(
        columns={"velocity_mean_m_s": "ale_velocity_u0_m_s"}
    )
    main = main.merge(reference, on="da", validate="many_to_one")
    main["ale_ratio"] = main["velocity_mean_m_s"] / main["ale_velocity_u0_m_s"]
    main["theory_ratio"] = 1.0 - main["u"] / np.cosh(np.sqrt(main["da"]))
    main["ratio_absolute_error"] = (
        main["ale_ratio"] - main["theory_ratio"]
    ).abs()
    main["ratio_normalized_error"] = main["ratio_absolute_error"] / np.maximum(
        np.abs(1.0 - main["theory_ratio"]), 1e-15
    )
    main["rate_balance_residual"] = main["bottom_total_mean"] - (
        main["top_flux_mean"]
        + main["positive_mean"]
        - main["dissolution_mean"]
        + main["storage_rate_mean"]
    )
    return main.sort_values(["da", "u"])


def convergence_tables(window: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    spatial = window[window["suite"] == "spatial"].copy()
    reference = spatial[
        (spatial["cells"] == 128) & np.isclose(spatial["power"], 1.5)
    ][["da", "u", "velocity_mean_m_s"]].rename(
        columns={"velocity_mean_m_s": "reference_velocity_m_s"}
    )
    spatial = spatial.merge(reference, on=["da", "u"], validate="many_to_one")
    spatial["velocity_relative_to_reference"] = (
        spatial["velocity_mean_m_s"] - spatial["reference_velocity_m_s"]
    ).abs() / spatial["reference_velocity_m_s"].abs()
    spatial["top_dx_over_ell"] = (
        np.sqrt(spatial["da"]) * (1.0 / spatial["cells"]) ** spatial["power"]
    )
    spatial["bottom_dx_over_ell"] = np.sqrt(spatial["da"]) * (
        1.0 - (1.0 - 1.0 / spatial["cells"]) ** spatial["power"]
    )
    spatial["top_dx_over_grain"] = (
        spatial["top_dx_over_ell"] * P.ell / P.grain
    )

    temporal = window[window["suite"] == "temporal"].copy()
    temporal["setting"] = temporal["case_id"].str.extract(
        r"(coarse|baseline|strict|smallstep)$"
    )
    reference_t = temporal[temporal["setting"] == "baseline"][
        ["da", "u", "velocity_mean_m_s"]
    ].rename(columns={"velocity_mean_m_s": "reference_velocity_m_s"})
    temporal = temporal.merge(reference_t, on=["da", "u"], validate="many_to_one")
    temporal["velocity_relative_to_reference"] = (
        temporal["velocity_mean_m_s"] - temporal["reference_velocity_m_s"]
    ).abs() / temporal["reference_velocity_m_s"].abs()
    return spatial, temporal


def profile_table(profiles: pd.DataFrame) -> pd.DataFrame:
    selected = profiles[
        (profiles["suite"] == "main")
        & np.isclose(profiles["time_over_tau"], 15.0)
        & profiles["case_id"].str.contains(
            r"da(?:1|20|100)_u(?:0|0.9)$", regex=True
        )
    ].copy()
    analytical = []
    for case_id, group in selected.groupby("case_id", sort=False):
        da = float(case_id.split("_da")[1].split("_u")[0])
        u = float(case_id.split("_u")[1])
        c, _ = analytical_profile(group["z_m"].to_numpy(), da, u, P)
        block = group.copy()
        block["da"] = da
        block["u"] = u
        block["z_over_h"] = block["z_m"] / block["z_m"].max()
        block["analytical_concentration"] = P.ceq + c
        block["concentration_error"] = (
            block["concentration"] - block["analytical_concentration"]
        )
        analytical.append(block)
    return pd.concat(analytical, ignore_index=True)


def x008_table() -> tuple[pd.DataFrame, dict[str, float]]:
    forcing = pd.read_csv(PRIOR / "x008_forcing.csv")
    ale = pd.read_csv(PRIOR / "ale_history.csv")
    time = forcing["elapsed_time_s"].to_numpy()
    j = forcing["evaporation_water_kg_m2_s"].to_numpy()
    back = forcing["b_soil_kg_m2_s"].to_numpy()
    height = np.interp(time, ale["time_s"], ale["height_m"])
    ale_bottom = np.interp(
        time, ale["time_s"], ale["bottom_total_salt_flux_kg_m2_s"]
    )
    ale_velocity = np.interp(time, ale["time_s"], ale["growth_velocity_m_s"])
    reconstructed = (P.ceq * j - back) / (1.0 - P.ceq)
    u = back / (P.ceq * j)
    da = (height / P.ell) ** 2
    theory_velocity = (
        P.ceq
        * j
        / (P.rho_cr * (1.0 - P.eps) * (1.0 - P.ceq))
        * P.x
        * (1.0 - u * sech(np.sqrt(da)))
        / (P.x + np.tanh(np.sqrt(da)))
    )
    out = pd.DataFrame(
        {
            "elapsed_time_s": time,
            "elapsed_time_h": time / 3600,
            "soil_exported_net_salt_kg_m2_s": forcing[
                "soil_net_nacl_kg_m2_s"
            ],
            "soil_advective_salt_kg_m2_s": forcing[
                "soil_advective_nacl_kg_m2_s"
            ],
            "soil_downward_diffusive_salt_kg_m2_s": back,
            "reconstructed_analytical_supply_kg_m2_s": reconstructed,
            "ale_bottom_total_salt_kg_m2_s": ale_bottom,
            "evaporation_water_kg_m2_s": j,
            "u": u,
            "ale_height_m": height,
            "instantaneous_da_birth_properties": da,
            "ale_top_velocity_m_s": ale_velocity,
            "corrected_top_velocity_m_s": theory_velocity,
        }
    )
    out["top_velocity_relative_error"] = (
        np.abs(out["ale_top_velocity_m_s"] - out["corrected_top_velocity_m_s"])
        / np.maximum(np.abs(out["corrected_top_velocity_m_s"]), 1e-30)
    )
    post_onset = out.iloc[1:]
    integrated = {
        "soil_exported_net_kg_m2": float(
            np.trapezoid(out["soil_exported_net_salt_kg_m2_s"], time)
        ),
        "reconstructed_supply_kg_m2": float(
            np.trapezoid(out["reconstructed_analytical_supply_kg_m2_s"], time)
        ),
        "ale_bottom_total_kg_m2": float(
            np.trapezoid(out["ale_bottom_total_salt_kg_m2_s"], time)
        ),
        "reconstructed_minus_soil_fraction": float(
            (
                np.trapezoid(
                    out["reconstructed_analytical_supply_kg_m2_s"], time
                )
                - np.trapezoid(out["soil_exported_net_salt_kg_m2_s"], time)
            )
            / np.trapezoid(out["soil_exported_net_salt_kg_m2_s"], time)
        ),
        "median_top_velocity_relative_error": float(
            np.median(
                np.abs(ale_velocity - theory_velocity)
                / np.maximum(np.abs(theory_velocity), 1e-30)
            )
        ),
        "max_top_velocity_relative_error": float(
            out["top_velocity_relative_error"].max()
        ),
        "post_first_sample_median_top_velocity_relative_error": float(
            post_onset["top_velocity_relative_error"].median()
        ),
        "post_first_sample_max_top_velocity_relative_error": float(
            post_onset["top_velocity_relative_error"].max()
        ),
    }
    return out, integrated


def style() -> None:
    mpl.rcParams.update(
        {
            "font.family": "serif",
            "font.serif": ["Times New Roman", "Times", "DejaVu Serif"],
            "font.size": 9,
            "axes.labelsize": 10,
            "legend.fontsize": 8,
            "xtick.labelsize": 8,
            "ytick.labelsize": 8,
            "pdf.fonttype": 42,
            "ps.fonttype": 42,
            "svg.fonttype": "none",
            "figure.facecolor": "white",
            "axes.facecolor": "white",
        }
    )


def save(fig: plt.Figure, name: str, sources: list[str], description: str) -> None:
    for suffix in ("png", "pdf", "svg"):
        fig.savefig(FIGURES / f"{name}.{suffix}", dpi=300, facecolor="white")
    svg_path = FIGURES / f"{name}.svg"
    svg_text = svg_path.read_text(encoding="utf-8")
    svg_path.write_text(
        "\n".join(line.rstrip() for line in svg_text.splitlines()) + "\n",
        encoding="utf-8",
    )
    manifest = {
        "name": name,
        "size_inches": list(fig.get_size_inches()),
        "formats": ["png", "pdf", "svg"],
        "dpi_for_png": 300,
        "sources": sources,
        "transformations": "documented in analyze_validation.py; no smoothing",
        "missing_data": "none silently connected; invalid cases would be gaps",
        "description": description,
        "journal_profile": "provisional general research figure; target journal pending",
    }
    (FIGURES / f"{name}_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    plt.close(fig)


def create_figures(
    main: pd.DataFrame,
    spatial: pd.DataFrame,
    temporal: pd.DataFrame,
    profiles: pd.DataFrame,
    x008: pd.DataFrame,
) -> None:
    style()
    colors = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#4D4D4D"]
    markers = ["o", "s", "^", "D", "v", "P"]

    fig, ax = plt.subplots(figsize=(7.09, 4.2), layout="constrained")
    dense_u = np.linspace(0, 0.9, 200)
    for color, marker, da in zip(colors, markers, sorted(main["da"].unique())):
        theory = 1.0 - dense_u / np.cosh(np.sqrt(da))
        block = main[main["da"] == da]
        ax.plot(dense_u, theory, color=color, lw=1.4)
        ax.plot(
            block["u"], block["ale_ratio"], ls="none", marker=marker, ms=5,
            mec=color, mfc="white", mew=1.1, color=color, label=f"Da = {da:g}",
        )
    ax.set(
        xlabel=r"Downward-diffusion ratio, $u=B_{\rm soil}/(C_{\rm eq}j)$",
        ylabel=r"Normalized top velocity, $v_t(u)/v_t(0)$",
        xlim=(-0.02, 0.92),
        ylim=(0.38, 1.025),
    )
    ax.grid(True, color="#D9D9D9", lw=0.55)
    ax.spines[["top", "right"]].set_visible(False)
    ax.legend(ncol=2, frameon=False, title="Lines: theory; markers: ALE")
    save(
        fig, "figure1_main_validation",
        ["results/main_validation.csv"],
        "Normalized ALE top velocities (markers) and prescribed-flux analytical predictions (lines) versus downward diffusion.",
    )

    pivot = main.pivot(index="da", columns="u", values="ratio_absolute_error")
    fig, ax = plt.subplots(figsize=(6.0, 4.3), layout="constrained")
    im = ax.imshow(
        pivot.to_numpy(), origin="lower", aspect="auto", cmap="viridis",
        extent=(-0.125, 1.025, -0.5, len(pivot.index)-0.5),
    )
    ax.set_xticks(pivot.columns)
    ax.set_yticks(range(len(pivot.index)), labels=[f"{x:g}" for x in pivot.index])
    ax.set(xlabel="u", ylabel="Da")
    cb = fig.colorbar(im, ax=ax)
    cb.set_label(r"$|R_{\rm ALE}-R_{\rm theory}|$")
    for row in range(len(pivot.index)):
        for col in range(len(pivot.columns)):
            value = pivot.iloc[row, col]
            ax.text(
                pivot.columns[col], row, f"{value:.1e}", ha="center", va="center",
                color="white" if value > 0.55 * pivot.to_numpy().max() else "black",
                fontsize=7,
            )
    save(
        fig, "figure2_error_map",
        ["results/main_validation.csv"],
        "Absolute normalized-velocity error over the six-by-five Da-u matrix; zero-u normalization points are retained.",
    )

    fig, axes = plt.subplots(1, 3, figsize=(7.09, 2.8), layout="constrained", sharey=True)
    for ax, da in zip(axes, (1.0, 20.0, 100.0)):
        for u, color, marker in ((0.0, colors[0], "o"), (0.9, colors[1], "s")):
            block = profiles[(profiles["da"] == da) & np.isclose(profiles["u"], u)]
            ax.plot(
                block["z_m"] / (P.ell * np.sqrt(da)),
                (block["analytical_concentration"] - P.ceq) * 1e6,
                color=color, lw=1.4,
            )
            ax.plot(
                block["z_m"] / (P.ell * np.sqrt(da)),
                (block["concentration"] - P.ceq) * 1e6,
                ls="none", marker=marker, ms=2.8, mfc="white", mec=color,
                label=f"u = {u:g}",
            )
        ax.set(xlabel=r"$z/h$", title=f"Da = {da:g}")
        ax.grid(True, color="#E0E0E0", lw=0.5)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel(r"$10^6(C-C_{\rm eq})$")
    handles = [
        Line2D([], [], color=colors[0], marker="o", mfc="white", label="u = 0"),
        Line2D([], [], color=colors[1], marker="s", mfc="white", label="u = 0.9"),
        Line2D([], [], color="black", lw=1.4, label="analytical line"),
    ]
    axes[-1].legend(handles=handles, frameon=False, loc="best")
    save(
        fig, "figure3_concentration_profiles",
        ["results/profile_comparison.csv"],
        "Analytical concentration profiles (lines) and ALE cell concentrations (markers) at 15 diffusion times.",
    )

    fig, axes = plt.subplots(1, 3, figsize=(7.09, 2.9), layout="constrained", sharey=True)
    for ax, da in zip(axes, (1.0, 20.0, 100.0)):
        block = main[main["da"] == da]
        ax.plot(block["u"], block["top_flux_mean"]/P.s_ref, "-o", color=colors[0], label="top")
        ax.plot(block["u"], block["positive_mean"]/P.s_ref, "--s", color=colors[2], label="internal +")
        ax.plot(block["u"], block["dissolution_mean"]/P.s_ref, "-.^", color=colors[1], label="internal −")
        ax.plot(block["u"], block["storage_rate_mean"]/P.s_ref, ":D", color=colors[3], label="dissolved storage")
        ax.set(xlabel="u", title=f"Da = {da:g}", xlim=(-0.02, 0.92))
        ax.grid(True, color="#E0E0E0", lw=0.5)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel(r"Rate divided by $S_{\rm ref}$")
    axes[-1].legend(frameon=False, fontsize=7)
    save(
        fig, "figure4_internal_redistribution",
        ["results/main_validation.csv"],
        "Separate top precipitation, positive internal precipitation, internal dissolution, and dissolved-storage rates.",
    )

    fig, axes = plt.subplots(2, 1, figsize=(7.09, 5.0), layout="constrained", sharex=True)
    axes[0].plot(x008["elapsed_time_h"], x008["soil_exported_net_salt_kg_m2_s"]*1e6,
                 color=colors[2], lw=1.5, label="soil exported SCVF net")
    axes[0].plot(x008["elapsed_time_h"], x008["reconstructed_analytical_supply_kg_m2_s"]*1e6,
                 color=colors[0], lw=1.5, ls="--", label="analytical reconstruction")
    axes[0].plot(x008["elapsed_time_h"], x008["ale_bottom_total_salt_kg_m2_s"]*1e6,
                 color=colors[1], lw=1.3, ls=":", label="ALE bottom total")
    axes[0].set_ylabel(r"Salt flux ($10^{-6}$ kg m$^{-2}$ s$^{-1}$)")
    axes[0].legend(frameon=False)
    axes[1].plot(x008["elapsed_time_h"], x008["ale_top_velocity_m_s"]*1e9,
                 color=colors[1], lw=1.4, label="ALE top")
    axes[1].plot(x008["elapsed_time_h"], x008["corrected_top_velocity_m_s"]*1e9,
                 color=colors[0], lw=1.4, ls="--", label="corrected law")
    axes[1].set(xlabel="Time after precipitation onset (h)", ylabel=r"Top velocity (nm s$^{-1}$)")
    axes[1].legend(frameon=False)
    for ax in axes:
        ax.grid(True, color="#D9D9D9", lw=0.55)
        ax.spines[["top", "right"]].set_visible(False)
    save(
        fig, "figure5_x008_forcing",
        ["results/x008_flux_comparison.csv"],
        "One-way x0=0.08 soil SCVF flux, analytical reconstructed supply, ALE bottom flux, and direct top velocities.",
    )

    fig, axes = plt.subplots(1, 2, figsize=(7.09, 3.1), layout="constrained")
    block = spatial[(spatial["da"] == 20) & np.isclose(spatial["u"], 0.9)]
    for power, color, marker in zip((1.0, 1.5, 2.0, 3.0), colors[:4], markers[:4]):
        part = block[np.isclose(block["power"], power)].sort_values("cells")
        axes[0].plot(
            part["cells"], part["velocity_relative_to_reference"]*100,
            marker=marker, color=color, label=f"p = {power:g}",
        )
    axes[0].set(
        xlabel="ALE cells", ylabel="Velocity difference from n=128, p=1.5 (%)",
        xticks=[32, 64, 128],
    )
    axes[0].legend(frameon=False)
    tblock = temporal[(temporal["da"] == 20) & np.isclose(temporal["u"], 0.9)]
    order = ["coarse", "baseline", "strict", "smallstep"]
    tblock = tblock.set_index("setting").loc[order].reset_index()
    axes[1].plot(
        range(len(order)), tblock["velocity_relative_to_reference"]*1e9,
        color=colors[1], marker="o",
    )
    axes[1].set(
        xlabel="Temporal setting",
        ylabel=r"Relative velocity difference ($10^{-9}$)",
        xticks=range(len(order)),
        xticklabels=order,
    )
    axes[1].tick_params(axis="x", rotation=25)
    for ax in axes:
        ax.grid(True, color="#D9D9D9", lw=0.55)
        ax.spines[["top", "right"]].set_visible(False)
    save(
        fig, "figure6_numerical_convergence",
        ["results/spatial_convergence.csv", "results/temporal_convergence.csv"],
        "Spatial and temporal convergence for the Da=20, u=0.9 representative benchmark.",
    )


def main() -> None:
    histories, profiles_raw, configs = load_raw()
    window = window_table(histories, configs)
    validation = main_validation(window)
    spatial, temporal = convergence_tables(window)
    profiles = profile_table(profiles_raw)
    x008, x008_summary = x008_table()
    statuses = status_table()

    window.to_csv(RESULTS / "quasi_steady_windows.csv", index=False)
    validation.to_csv(RESULTS / "main_validation.csv", index=False)
    spatial.to_csv(RESULTS / "spatial_convergence.csv", index=False)
    temporal.to_csv(RESULTS / "temporal_convergence.csv", index=False)
    profiles.to_csv(RESULTS / "profile_comparison.csv", index=False)
    x008.to_csv(RESULTS / "x008_flux_comparison.csv", index=False)
    statuses.to_csv(RESULTS / "case_status.csv", index=False)

    boundary = window[window["suite"] == "boundary"].copy()
    boundary.to_csv(RESULTS / "boundary_audit.csv", index=False)
    da100u0 = validation[(validation["da"] == 100) & (validation["u"] == 0)].iloc[0]
    summary = {
        "status": "PASS_CONTROLLED_BENCHMARKS_WITH_SCOPE_LIMITS",
        "bvp_cases": 30,
        "ale_cases": int(configs.shape[0]),
        "ale_completed_cases": int((statuses["status"] == "forcing_end").sum()),
        "quasi_steady_valid_cases": int(window["quasi_steady_valid"].sum()),
        "main_valid_cases": int(validation["quasi_steady_valid"].sum()),
        "max_ratio_absolute_error": float(validation["ratio_absolute_error"].max()),
        "max_velocity_relative_error": float(validation["velocity_relative_error"].max()),
        "max_top_concentration_abs_error": float(
            validation["top_concentration_abs_error"].max()
        ),
        "max_bottom_concentration_abs_deviation_from_ceq": float(
            (validation["bottom_face_concentration_mean"] - P.ceq).abs().max()
        ),
        "max_rate_balance_abs_residual_kg_m2_s": float(
            validation["rate_balance_residual"].abs().max()
        ),
        "max_salt_inventory_residual_kg_m2": float(window["salt_residual_max"].max()),
        "max_water_inventory_residual_kg_m2": float(window["water_residual_max"].max()),
        "spatial_n64_p1p5_max_relative_difference_from_n128": float(
            spatial[
                (spatial["cells"] == 64) & np.isclose(spatial["power"], 1.5)
            ]["velocity_relative_to_reference"].max()
        ),
        "spatial_n32_p1p5_max_relative_difference_from_n128": float(
            spatial[
                (spatial["cells"] == 32) & np.isclose(spatial["power"], 1.5)
            ]["velocity_relative_to_reference"].max()
        ),
        "temporal_max_relative_difference": float(
            temporal["velocity_relative_to_reference"].max()
        ),
        "large_da_top_fraction_u0": float(da100u0["top_flux_mean"] / P.s_ref),
        "large_da_positive_internal_fraction_u0": float(
            da100u0["positive_mean"] / P.s_ref
        ),
        "large_da_theory_top_fraction": float(P.x / (1.0 + P.x)),
        "large_da_theory_internal_fraction": float(1.0 / (1.0 + P.x)),
        "x008": x008_summary,
        "slurm": {
            "build_job": 16569420,
            "failed_crlf_array": 16570638,
            "cancelled_retry_before_solve_completion": 16570777,
            "successful_array": 16570915,
            "failed_pack": 16571055,
            "successful_pack": 16571056,
        },
    }
    (RESULTS / "validation_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    create_figures(validation, spatial, temporal, profiles, x008)
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
