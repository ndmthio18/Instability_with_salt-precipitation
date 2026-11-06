#!/usr/bin/env python3
"""Compare the single authorized P1-R2 2D run with its frozen R24 baseline."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np


REFERENCE_TIMES = {
    "global_lsa_crossing_s": 31409.1577927,
    "mode28_lsa_crossing_s": 45474.3545908,
    "original_detection_s": 50394.3336331,
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_sigma(path: Path) -> dict[str, np.ndarray]:
    values = np.loadtxt(path, comments="#")
    keys = ("time_s", "surface_mean_x", "surface_sigma_x", "n_top", "fourier_x", "projection_x")
    return {key: values[:, index] for index, key in enumerate(keys)}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def moving_average(values: np.ndarray, width: int = 41) -> np.ndarray:
    width = min(width, len(values) // 2 * 2 - 1)
    if width < 3:
        return values
    kernel = np.ones(width) / width
    padded = np.pad(values, (width // 2, width // 2), mode="edge")
    return np.convolve(padded, kernel, mode="valid")


def growth_rate(series: dict[str, np.ndarray]) -> np.ndarray:
    raw = np.gradient(np.log(np.abs(series["projection_x"])), series["time_s"])
    return moving_average(raw)


def nearest_index(times: np.ndarray, target: float) -> int:
    return int(np.argmin(np.abs(times - target)))


def data_array(parent: ET.Element, name: str) -> np.ndarray:
    element = parent.find(f"DataArray[@Name='{name}']")
    if element is None or element.text is None:
        raise KeyError(name)
    values = np.fromstring(element.text, sep=" ")
    components = int(element.attrib.get("NumberOfComponents", "1"))
    return values.reshape((-1, components)) if components > 1 else values


def read_point_field(path: Path) -> dict[str, np.ndarray]:
    piece = ET.parse(path).getroot().find(".//Piece")
    if piece is None:
        raise ValueError(f"missing Piece in {path}")
    points = data_array(piece.find("Points"), "Coordinates")
    point_data = piece.find("PointData")
    if point_data is None:
        raise ValueError(f"missing PointData in {path}")
    return {
        "x_m": points[:, 0],
        "z_m": points[:, 1],
        "rho_kg_m3": data_array(point_data, "rhoW"),
        "x_salt": data_array(point_data, "moleFraction"),
        "w_salt": data_array(point_data, "massFraction"),
    }


def read_pvd(directory: Path) -> list[tuple[float, Path]]:
    pvd = next(directory.glob("*.pvd"))
    root = ET.parse(pvd).getroot()
    return [
        (float(node.attrib["timestep"]), directory / node.attrib["file"])
        for node in root.findall(".//DataSet")
    ]


def final_quarter_slope(series: dict[str, np.ndarray]) -> float:
    start = int(0.75 * len(series["time_s"]))
    return float(np.polyfit(series["time_s"][start:], np.log(np.abs(series["projection_x"][start:])), 1)[0])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--new-run", type=Path, required=True)
    parser.add_argument("--baseline-run", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    new = read_sigma(next(args.new_run.glob("*_sigmatop.csv")))
    old = read_sigma(next(args.baseline_run.glob("*_sigmatop.csv")))
    if not np.allclose(new["time_s"], old["time_s"], rtol=0.0, atol=1e-7):
        raise ValueError("new and frozen R24 time grids do not match")

    new_norm = new["projection_x"] / new["projection_x"][0]
    old_norm = old["projection_x"] / old["projection_x"][0]
    new_growth = growth_rate(new)
    old_growth = growth_rate(old)

    comparison_rows = []
    for label, reference in {**REFERENCE_TIMES, "final_s": new["time_s"][-1]}.items():
        index = nearest_index(new["time_s"], reference)
        comparison_rows.append(
            {
                "diagnostic_time": label,
                "reference_time_s": reference,
                "actual_time_s": new["time_s"][index],
                "new_projection_x": new["projection_x"][index],
                "old_projection_x": old["projection_x"][index],
                "new_normalized_amplitude": new_norm[index],
                "old_normalized_amplitude": old_norm[index],
                "new_growth_rate_s-1": new_growth[index],
                "old_growth_rate_s-1": old_growth[index],
                "new_surface_mean_x": new["surface_mean_x"][index],
                "old_surface_mean_x": old["surface_mean_x"][index],
                "new_surface_sigma_x": new["surface_sigma_x"][index],
                "old_surface_sigma_x": old["surface_sigma_x"][index],
            }
        )
    with (args.output / "instability_comparison.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=comparison_rows[0].keys())
        writer.writeheader()
        writer.writerows(comparison_rows)

    budget = read_csv(args.new_run / "component_budget.csv")
    residual_keys = ("total_residual_kg", "salt_residual_kg", "water_residual_kg")
    max_residuals = {
        key: max(abs(float(row[key])) for row in budget) for key in residual_keys
    }
    accepted = read_csv(args.new_run / "accepted_steps.csv")

    cuts = read_csv(args.new_run / "cut_fluxes.csv")
    cut_names = sorted({row["cut"] for row in cuts if row["cut"].startswith("requested_")})
    with (args.output / "downward_salt_flux.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=cuts[0].keys())
        writer.writeheader()
        writer.writerows(row for row in cuts if row["cut"] in cut_names)

    fig, axes = plt.subplots(2, 2, figsize=(12, 8), constrained_layout=True)
    time_h = new["time_s"] / 3600.0
    axes[0, 0].plot(time_h, old_norm, label="Frozen normalized R24")
    axes[0, 0].plot(time_h, new_norm, label="Conservative Mode A")
    axes[0, 0].set(ylabel="A(t) / A(0)")
    axes[0, 0].legend()
    axes[0, 1].plot(time_h, old_growth * 3600.0, label="Frozen normalized R24")
    axes[0, 1].plot(time_h, new_growth * 3600.0, label="Conservative Mode A")
    axes[0, 1].set(ylabel="Smoothed d ln(A) / dt (h$^{-1}$)")
    axes[1, 0].plot(time_h, old["surface_mean_x"], label="Frozen normalized R24")
    axes[1, 0].plot(time_h, new["surface_mean_x"], label="Conservative Mode A")
    axes[1, 0].set(xlabel="Time (h)", ylabel="Mean surface NaCl mole fraction")
    axes[1, 1].plot(time_h, old["surface_sigma_x"], label="Frozen normalized R24")
    axes[1, 1].plot(time_h, new["surface_sigma_x"], label="Conservative Mode A")
    axes[1, 1].set(xlabel="Time (h)", ylabel="Surface horizontal standard deviation")
    for axis in axes.flat:
        axis.grid(alpha=0.25)
        for reference in REFERENCE_TIMES.values():
            axis.axvline(reference / 3600.0, color="0.6", lw=0.7, ls="--")
    fig.savefig(args.output / "INSTABILITY_COMPARISON.png", dpi=200)
    plt.close(fig)

    pvd = read_pvd(args.new_run)
    field_targets = {
        "initial": 0.0,
        "global_lsa": REFERENCE_TIMES["global_lsa_crossing_s"],
        "mode28_lsa": REFERENCE_TIMES["mode28_lsa_crossing_s"],
        "final": new["time_s"][-1],
    }
    selected_fields = []
    for label, target in field_targets.items():
        actual_time, path = min(pvd, key=lambda record: abs(record[0] - target))
        selected_fields.append((label, actual_time, read_point_field(path)))

    fig, axes = plt.subplots(1, 4, figsize=(16, 4.1), sharex=True, sharey=True, constrained_layout=True)
    all_density = np.concatenate([item[2]["rho_kg_m3"] for item in selected_fields])
    levels = np.linspace(float(all_density.min()), float(all_density.max()), 25)
    last = None
    for axis, (label, actual_time, field) in zip(axes, selected_fields):
        mask = field["z_m"] >= 0.18
        triangulation = mtri.Triangulation(field["x_m"][mask] * 1000.0, (0.2 - field["z_m"][mask]) * 1000.0)
        last = axis.tricontourf(triangulation, field["rho_kg_m3"][mask], levels=levels, cmap="viridis", extend="both")
        axis.set_title(f"{label}\n{actual_time / 3600.0:.2f} h")
        axis.set_xlabel("x (mm)")
        axis.invert_yaxis()
    axes[0].set_ylabel("Depth below surface (mm)")
    fig.colorbar(last, ax=axes, label="Liquid density (kg m$^{-3}$)")
    fig.savefig(args.output / "DENSITY_FIELDS.png", dpi=200)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 5), constrained_layout=True)
    for cut in cut_names:
        rows = [row for row in cuts if row["cut"] == cut]
        depth_mm = (0.2 - float(rows[0]["cut_z_m"])) * 1000.0
        ax.plot(
            [float(row["time_s"]) / 3600.0 for row in rows],
            [float(row["salt_down_kg_s"]) for row in rows],
            label=f"{depth_mm:.2f} mm",
        )
    ax.set(xlabel="Time (h)", ylabel="Conservative downward NaCl flux (kg s$^{-1}$)")
    ax.grid(alpha=0.25)
    ax.legend(title="Depth below surface", ncol=2)
    fig.savefig(args.output / "DOWNWARD_SALT_FLUX.png", dpi=200)
    plt.close(fig)

    summary = {
        "schema_version": 1,
        "run_scope": "single bounded one-wavelength Case 36877 comparison; no production campaign",
        "grid": {"cells_x": 24, "cells_z": 200, "vertical_grading": -0.98},
        "baseline": "frozen normalized R24-M1-full fixed-grid reference; not continuum truth",
        "accepted_steps": len(accepted),
        "retried_steps": sum(int(row["retried"]) > 0 for row in accepted),
        "final_time_s": float(new["time_s"][-1]),
        "max_component_budget_residuals_kg": max_residuals,
        "final_normalized_amplitude": {"conservative_mode_A": float(new_norm[-1]), "frozen_normalized": float(old_norm[-1])},
        "surface_sigma_minimum_time_s": {
            "conservative_mode_A": float(new["time_s"][int(np.argmin(new["surface_sigma_x"]))]),
            "frozen_normalized": float(old["time_s"][int(np.argmin(old["surface_sigma_x"]))]),
        },
        "final_quarter_log_growth_rate_s-1": {"conservative_mode_A": final_quarter_slope(new), "frozen_normalized": final_quarter_slope(old)},
        "final_surface_mean_x": {"conservative_mode_A": float(new["surface_mean_x"][-1]), "frozen_normalized": float(old["surface_mean_x"][-1])},
        "final_relative_differences": {
            "normalized_amplitude": float(new_norm[-1] / old_norm[-1] - 1.0),
            "final_quarter_log_growth_rate": float(final_quarter_slope(new) / final_quarter_slope(old) - 1.0),
            "surface_mean_x": float(new["surface_mean_x"][-1] / old["surface_mean_x"][-1] - 1.0),
        },
        "input_sha256": sha256(args.new_run / "V8_case36877_R24_mode28_mass_A.input"),
        "eigenprofile_sha256": sha256(args.new_run / "eigenprofile.dat"),
        "claims": {
            "component_conservation": "supported at this fixed grid and accepted time-step sequence",
            "numerical_convergence": "not established by this single R24 run",
            "physical_validation": "not established",
            "instability_growth": "bounded fixed-grid comparison only; reference times are diagnostics, not guaranteed predictions",
        },
    }
    (args.output / "instability_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
