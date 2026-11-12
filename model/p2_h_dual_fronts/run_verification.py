#!/usr/bin/env python3
"""Run the analytical sweep and conservative dual-front verification matrix."""

from __future__ import annotations

import csv
import json
import math
import platform
from pathlib import Path
import sys

import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import PIL
from PIL import Image

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from model.p2_h_dual_fronts.analytical import (  # noqa: E402
    AnalyticalParameters,
    finite_difference_profile,
    modified_velocity_factor,
    normalized_profile,
    prat_equivalent_u,
    prat_velocity_factor,
    sech,
    velocity_ratio,
)
from model.p2_h_dual_fronts.dual_front import (  # noqa: E402
    BoundaryKinematics,
    DualFrontModel,
)


OUT = REPO / "analysis" / "p2_h_dual_moving_fronts_20261112"
FIG = OUT / "figures"
DATA = OUT / "data"
DA_VALUES = (0.1, 1.0, 4.0, 10.0, 20.0, 45.0, 100.0)
U_VALUES = (0.0, 0.25, 0.5, 0.75, 1.0, 1.25)
X_VALUES = (0.5, 1.0, 2.0, 5.0)
LINESTYLES = ("-", "--", "-.", ":")


def save_figure(fig, name: str) -> None:
    """Export an opaque raster preview and an editable vector companion."""

    png = FIG / f"{name}.png"
    fig.savefig(png, dpi=220, facecolor="white", transparent=False)
    with Image.open(png) as rendered:
        if rendered.mode != "RGB":
            rendered.convert("RGB").save(png, dpi=(220, 220))
    fig.savefig(FIG / f"{name}.pdf", facecolor="white")


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"no rows for {path}")
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def analytical_verification() -> dict:
    rows: list[dict] = []
    worst = 0.0
    for da in DA_VALUES:
        for x in X_VALUES:
            for u in U_VALUES:
                xi, numerical = finite_difference_profile(da, x, u, cells=500)
                exact_profile = normalized_profile(xi, da, x, u)
                exact_factor = modified_velocity_factor(da, x, u)
                numerical_factor = x * numerical[-1]
                error = abs(numerical_factor - exact_factor)
                profile_error = float(np.max(np.abs(numerical - exact_profile)))
                worst = max(worst, error, profile_error)
                rows.append(
                    {
                        "Da": da,
                        "X": x,
                        "u": u,
                        "exact_velocity_factor": exact_factor,
                        "fd_velocity_factor": numerical_factor,
                        "velocity_abs_error": error,
                        "profile_max_abs_error": profile_error,
                        "velocity_ratio_exact": velocity_ratio(da, u),
                        "velocity_ratio_from_factors": exact_factor
                        / modified_velocity_factor(da, x, 0.0),
                        "large_Da_ratio": 1.0,
                        "Prat_factor": prat_velocity_factor(da, x),
                        "Prat_equivalent_u": prat_equivalent_u(da, x),
                    }
                )
    write_csv(DATA / "analytical_sweep.csv", rows)
    if worst >= 2.5e-4:
        raise AssertionError(f"finite-difference analytical check failed: {worst}")
    return {"worst_finite_difference_abs_error": worst, "rows": len(rows)}


def analytical_figures() -> None:
    colors = plt.cm.viridis(np.linspace(0.05, 0.95, len(DA_VALUES)))
    u_dense = np.linspace(0.0, 1.25, 200)
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for index, (da, color) in enumerate(zip(DA_VALUES, colors)):
        ax.plot(
            u_dense,
            [velocity_ratio(da, u) for u in u_dense],
            color=color,
            linestyle=LINESTYLES[index % len(LINESTYLES)],
            label=f"Da={da:g}",
        )
    ax.axhline(0.0, color="black", lw=0.8)
    ax.set(xlabel=r"$u=B_{soil}/(C_{eq}j)$", ylabel=r"$v_t(u)/v_t(0)$", ylim=(-0.25, 1.04))
    ax.grid(alpha=0.25)
    ax.legend(ncol=2, fontsize=8)
    fig.tight_layout()
    save_figure(fig, "normalized_top_growth_vs_u")
    plt.close(fig)

    da_dense = np.logspace(-1, 2, 300)
    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for index, u in enumerate(U_VALUES[1:]):
        reduction = [1.0 - velocity_ratio(float(da), u) for da in da_dense]
        ax.loglog(da_dense, reduction, linestyle=LINESTYLES[index % len(LINESTYLES)], label=f"u={u:g}")
    ax.set(xlabel=r"$Da$", ylabel=r"$1-v_t(u)/v_t(0)$")
    ax.grid(which="both", alpha=0.25)
    ax.legend(fontsize=8)
    fig.tight_layout()
    save_figure(fig, "growth_reduction_vs_Da")
    plt.close(fig)

    xi = np.linspace(0.0, 1.0, 300)
    fig, axes = plt.subplots(1, 3, figsize=(12.0, 3.8), sharey=True)
    for ax, da in zip(axes, (1.0, 20.0, 100.0)):
        for index, u in enumerate((0.0, 0.75, 1.25)):
            ax.plot(normalized_profile(xi, da, 2.0, u), xi, linestyle=LINESTYLES[index], label=f"u={u:g}")
        ax.axvline(0.0, color="0.5", lw=0.8)
        ax.set(title=f"Da={da:g}", xlabel=r"normalized $C-C_{eq}$")
        ax.grid(alpha=0.25)
    axes[0].set_ylabel(r"$ξ=(z-z_b)/h$")
    axes[-1].legend(fontsize=8)
    fig.tight_layout()
    save_figure(fig, "concentration_profiles")
    plt.close(fig)

    fig, axes = plt.subplots(1, 3, figsize=(12.0, 3.8), sharey=True)
    for ax, da in zip(axes, (4.0, 20.0, 100.0)):
        for index, u in enumerate((0.0, 0.75, 1.25)):
            reaction = normalized_profile(xi, da, 2.0, u)
            ax.plot(reaction, xi, linestyle=LINESTYLES[index], label=f"u={u:g}")
        ax.axvline(0.0, color="black", lw=0.8)
        ax.set(title=f"Da={da:g}", xlabel="signed reaction (normalized)")
        ax.grid(alpha=0.25)
    axes[0].set_ylabel(r"$ξ$")
    axes[-1].legend(fontsize=8)
    fig.tight_layout()
    save_figure(fig, "signed_internal_reaction_profiles")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    for index, x in enumerate(X_VALUES):
        exact = [modified_velocity_factor(float(da), x, 0.0) for da in da_dense]
        asymptote = x / (x + 1.0)
        ax.semilogx(da_dense, exact, linestyle=LINESTYLES[index], label=f"exact, X={x:g}")
        ax.axhline(asymptote, color=ax.lines[-1].get_color(), ls=":", alpha=0.7)
    ax.set(xlabel=r"$Da$", ylabel="growth-efficiency factor")
    ax.grid(which="both", alpha=0.25)
    ax.legend(ncol=2, fontsize=8)
    fig.tight_layout()
    save_figure(fig, "finite_Da_vs_asymptotic_growth")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(7.2, 4.8))
    ax.loglog(da_dense, 1.0 / np.sqrt(da_dense), color="tab:purple")
    ax.scatter(DA_VALUES, [1.0 / math.sqrt(da) for da in DA_VALUES], color="black", s=18)
    ax.set(xlabel=r"$Da$", ylabel=r"penetration ratio $ℓ_r/h=Da^{-1/2}$")
    ax.grid(which="both", alpha=0.25)
    fig.tight_layout()
    save_figure(fig, "penetration_length_ratio")
    plt.close(fig)


def conservation_error(model: DualFrontModel, initial, final) -> tuple[float, float]:
    a = model.inventory(initial)
    b = model.inventory(final)
    return (
        b["total_salt_with_external_kg_m2"] - a["total_salt_with_external_kg_m2"],
        b["total_water_with_external_kg_m2"] - a["total_water_with_external_kg_m2"],
    )


def run_case(name: str, model: DualFrontModel, initial, motion, duration: float, dt: float) -> dict:
    final, history = model.run(initial, motion, duration, dt)
    salt_error, water_error = conservation_error(model, initial, final)
    write_csv(DATA / f"{name}_history.csv", history)
    before = model.inventory(initial)
    after = model.inventory(final)
    _, eps_final, c_final = model.primitive(final)
    return {
        "test": name,
        "cells": model.cells,
        "dt_s": dt,
        "duration_s": duration,
        "z_bottom_final_m": final.z_bottom_m,
        "z_top_final_m": final.z_top_m,
        "height_initial_m": before["height_m"],
        "height_final_m": after["height_m"],
        "solid_initial_kg_m2": before["solid_kg_m2"],
        "solid_final_kg_m2": after["solid_kg_m2"],
        "average_porosity_initial": before["average_porosity"],
        "average_porosity_final": after["average_porosity"],
        "porosity_min_final": float(np.min(eps_final)),
        "porosity_max_final": float(np.max(eps_final)),
        "concentration_min_final": float(np.min(c_final)),
        "concentration_max_final": float(np.max(c_final)),
        "salt_residual_kg_m2": salt_error,
        "water_residual_kg_m2": water_error,
        "top_solid_kg_m2": final.cumulative_top_solid_kg_m2,
        "base_dissolution_kg_m2": final.cumulative_base_dissolution_kg_m2,
        "bulk_signed_reaction_kg_m2": final.cumulative_bulk_reaction_kg_m2,
    }


def dual_front_tests() -> tuple[list[dict], list[dict]]:
    c0 = AnalyticalParameters().equilibrium_mass_fraction
    h0 = 100.0e-6
    cases: list[dict] = []

    # A: operator-level recovery of the original fixed-bottom kinematics.
    model_a = DualFrontModel(64)
    eps_a = 0.43 + 0.01 * np.sin(2.0 * np.pi * model_a.centers)
    c_a = c0 + 0.002 * np.cos(np.pi * model_a.centers)
    initial_a = model_a.initial_state(h0, c_a, eps_a)
    motion_a = BoundaryKinematics(1.0e-9, 0.0)
    rates_a = model_a.rates(initial_a, motion_a)
    # Independent M0 geometry: w=xi*v_t and H=-w*eps_upwind.
    _, eps, conc = model_a.primitive(initial_a)
    w = model_a.faces * motion_a.top_velocity_m_s
    hflux = np.empty(model_a.cells + 1)
    hflux[0] = 0.0
    hflux[-1] = -motion_a.top_velocity_m_s * model_a.top_porosity
    hflux[1:-1] = -w[1:-1] * eps[1:]
    expected_pore = hflux[:-1] - hflux[1:]
    if not np.allclose(rates_a.pore, expected_pore, rtol=0.0, atol=1.0e-18):
        raise AssertionError("Test A fixed-bottom pore operator mismatch")
    cface = np.empty(model_a.cells + 1)
    cface[0] = conc[0]
    cface[-1] = conc[-1]
    cface[1:-1] = conc[1:]
    gtop = model_a.rho_cr * (1.0 - model_a.top_porosity) * motion_a.top_velocity_m_s
    expected_salt_flux = model_a.rho * cface * hflux + gtop
    expected_salt = expected_salt_flux[:-1] - expected_salt_flux[1:]
    if not np.allclose(rates_a.salt, expected_salt, rtol=0.0, atol=1.0e-18):
        raise AssertionError("Test A fixed-bottom dissolved-salt operator mismatch")
    cases.append(run_case("test_A_baseline_recovery", model_a, initial_a, motion_a, 20.0, 0.25))

    model_b = DualFrontModel(80)
    initial_b = model_b.initial_state(h0, c0)
    cases.append(
        run_case("test_B_pure_top_growth", model_b, initial_b, BoundaryKinematics(2e-9, 0.0), 1000.0, 1.0)
    )
    if max(
        abs(cases[-1]["porosity_min_final"] - 0.43),
        abs(cases[-1]["porosity_max_final"] - 0.43),
        abs(cases[-1]["concentration_min_final"] - c0),
        abs(cases[-1]["concentration_max_final"] - c0),
    ) > 2.0e-12:
        raise AssertionError("Test B failed uniform concentration/porosity recovery")

    model_c = DualFrontModel(80)
    initial_c = model_c.initial_state(h0, c0)
    cases.append(
        run_case("test_C_pure_bottom_retreat", model_c, initial_c, BoundaryKinematics(0.0, 1e-9, 0.4), 1000.0, 1.0)
    )

    # Equal velocities with unequal interface porosities: constant thickness
    # but a non-zero solid-inventory change, as required by the task.
    model_d = DualFrontModel(80, top_porosity=0.43, bottom_porosity=0.50)
    initial_d = model_d.initial_state(h0, c0, porosity=0.465)
    cases.append(
        run_case("test_D_equal_front_velocities", model_d, initial_d, BoundaryKinematics(1e-9, 1e-9, 0.5), 1000.0, 1.0)
    )

    model_e = DualFrontModel(80)
    initial_e = model_e.initial_state(h0, c0)
    cases.append(
        run_case("test_E_thinning_translation", model_e, initial_e, BoundaryKinematics(0.5e-9, 1e-9, 0.5), 1000.0, 1.0)
    )

    def signed_reaction(xi, _eps, _c, _time):
        return np.where(xi > 0.65, 0.02, np.where((xi > 0.20) & (xi < 0.45), -0.01, 0.0))

    model_f = DualFrontModel(80, reaction_law=signed_reaction, bulk_exclusion_fraction=0.15)
    initial_f = model_f.initial_state(h0, c0)
    cases.append(
        run_case("test_F_signed_internal_reaction", model_f, initial_f, BoundaryKinematics(1e-9, 0.5e-9, 0.5), 100.0, 0.5)
    )

    for case in cases:
        if abs(case["salt_residual_kg_m2"]) > 2.0e-12 or abs(case["water_residual_kg_m2"]) > 2.0e-12:
            raise AssertionError(f"conservation failure in {case['test']}: {case}")
    if abs(cases[3]["height_final_m"] - cases[3]["height_initial_m"]) > 1.0e-14:
        raise AssertionError("Test D did not preserve thickness")
    if abs(cases[3]["solid_final_kg_m2"] - cases[3]["solid_initial_kg_m2"]) < 1.0e-8:
        raise AssertionError("Test D failed to expose solid change at constant thickness")
    if not cases[4]["height_final_m"] < cases[4]["height_initial_m"]:
        raise AssertionError("Test E did not thin")

    def convergence_reaction(xi, _eps, c, _time):
        return 8.0 * (c - c0) * (0.25 + xi)

    convergence: list[dict] = []
    for cells in (40, 80, 160):
        for dt in (0.5, 0.25):
            model = DualFrontModel(
                cells,
                reaction_law=convergence_reaction,
                bulk_exclusion_fraction=0.1,
            )
            eps = 0.43 + 0.015 * np.sin(2.0 * np.pi * model.centers)
            conc = c0 + 0.003 * np.cos(np.pi * model.centers)
            initial = model.initial_state(h0, conc, eps)
            result = run_case(
                f"test_G_n{cells}_dt{str(dt).replace('.', 'p')}",
                model,
                initial,
                BoundaryKinematics(1.0e-9, 0.35e-9, 0.35),
                100.0,
                dt,
            )
            convergence.append(result)
    reference = next(row for row in convergence if row["cells"] == 160 and row["dt_s"] == 0.25)
    for row in convergence:
        for key in ("z_bottom_final_m", "z_top_final_m", "height_final_m", "solid_final_kg_m2", "average_porosity_final"):
            row[f"abs_error_{key}"] = abs(row[key] - reference[key])
    write_csv(DATA / "dual_front_test_summary.csv", cases)
    write_csv(DATA / "grid_time_convergence.csv", convergence)
    return cases, convergence


def dual_front_figures(cases: list[dict], convergence: list[dict]) -> None:
    selected = ("test_B_pure_top_growth", "test_C_pure_bottom_retreat", "test_D_equal_front_velocities", "test_E_thinning_translation")
    fig, axes = plt.subplots(2, 2, figsize=(10.5, 7.2))
    for index, name in enumerate(selected):
        rows = list(csv.DictReader((DATA / f"{name}_history.csv").open(encoding="utf-8")))
        t = np.array([float(r["time_s"]) for r in rows])
        label = name.replace("test_", "").replace("_", " ")
        style = LINESTYLES[index]
        axes[0, 0].plot(t, [float(r["z_top_m"]) * 1e6 for r in rows], linestyle=style, label=label)
        axes[0, 1].plot(t, [float(r["z_bottom_m"]) * 1e6 for r in rows], linestyle=style, label=label)
        axes[1, 0].plot(t, [float(r["height_m"]) * 1e6 for r in rows], linestyle=style, label=label)
        axes[1, 1].plot(t, [float(r["solid_kg_m2"]) for r in rows], linestyle=style, label=label)
    axes[0, 0].set(ylabel=r"$z_t$ (μm)")
    axes[0, 1].set(ylabel=r"$z_b$ (μm)")
    axes[1, 0].set(xlabel="time (s)", ylabel=r"$h$ (μm)")
    axes[1, 1].set(xlabel="time (s)", ylabel=r"solid (kg m$^{-2}$)")
    for ax in axes.flat:
        ax.grid(alpha=0.25)
    axes[0, 0].legend(fontsize=7)
    fig.tight_layout()
    save_figure(fig, "dual_front_geometry_and_inventory")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.0))
    for dt, marker in ((0.5, "o"), (0.25, "s")):
        rows = [r for r in convergence if r["dt_s"] == dt]
        axes[0].loglog([r["cells"] for r in rows], [max(r["abs_error_solid_final_kg_m2"], 1e-18) for r in rows], marker=marker, label=f"dt={dt:g} s")
        axes[1].loglog([r["cells"] for r in rows], [max(abs(r["salt_residual_kg_m2"]), 1e-18) for r in rows], marker=marker, label=f"dt={dt:g} s")
    axes[0].set(xlabel="ALE cells", ylabel="solid-inventory error (kg m$^{-2}$)")
    axes[1].set(xlabel="ALE cells", ylabel="salt residual (kg m$^{-2}$)")
    for ax in axes:
        ax.grid(which="both", alpha=0.25)
        ax.legend(fontsize=8)
    fig.tight_layout()
    save_figure(fig, "grid_time_convergence")
    plt.close(fig)


def main() -> None:
    FIG.mkdir(parents=True, exist_ok=True)
    DATA.mkdir(parents=True, exist_ok=True)
    params = AnalyticalParameters()
    (DATA / "actual_3p5um_parameters.json").write_text(
        json.dumps(params.to_record(), indent=2) + "\n", encoding="utf-8"
    )
    analytical = analytical_verification()
    analytical_figures()
    cases, convergence = dual_front_tests()
    dual_front_figures(cases, convergence)
    summary = {
        "status": "PASS_DIAGNOSTIC_NOT_PREDICTIVE",
        "baseline_commit": "fdd4d1f306c09caf3f44cac39f9cc7431aadfaa0",
        "environment": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "matplotlib": matplotlib.__version__,
            "pillow": PIL.__version__,
        },
        "analytical": analytical,
        "actual_parameters": params.to_record(),
        "dual_front_tests": cases,
        "convergence_rows": len(convergence),
        "large_Da_reductions_at_u_1p25": {
            str(da): 1.0 - velocity_ratio(da, 1.25) for da in DA_VALUES
        },
        "finite_Da_zero_flux_minus_Prat": {
            str(da): modified_velocity_factor(da, 2.0, 0.0) - prat_velocity_factor(da, 2.0)
            for da in DA_VALUES
        },
    }
    (OUT / "verification_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
