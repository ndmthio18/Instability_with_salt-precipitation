#!/usr/bin/env python3
"""Extract and compare matched 1D P1-R2 background runs."""

from __future__ import annotations

import argparse
import csv
import json
import math
import xml.etree.ElementTree as ET
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


M_WATER = 0.01801528
M_SALT = 0.05844
POROSITY = 0.4057687912195745
REFERENCE_TIMES = {
    "global_lsa_crossing_s": 31409.1577927,
    "mode28_lsa_crossing_s": 45474.3545908,
    "original_detection_s": 50394.3336331,
}


def mass_fraction(x_salt: np.ndarray) -> np.ndarray:
    mixture = (1.0 - x_salt) * M_WATER + x_salt * M_SALT
    return x_salt * M_SALT / mixture


def read_pvd(directory: Path) -> list[tuple[float, Path]]:
    pvd = next(directory.glob("*.pvd"))
    root = ET.parse(pvd).getroot()
    records = []
    for dataset in root.findall(".//DataSet"):
        records.append(
            (
                float(dataset.attrib["timestep"]),
                directory / dataset.attrib["file"],
            )
        )
    return records


def data_array(parent: ET.Element, name: str) -> np.ndarray:
    element = parent.find(f"DataArray[@Name='{name}']")
    if element is None or element.text is None:
        raise KeyError(name)
    values = np.fromstring(element.text, sep=" ")
    components = int(element.attrib.get("NumberOfComponents", "1"))
    return values.reshape((-1, components)) if components > 1 else values


def grouped_mean(z: np.ndarray, values: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    rounded = np.round(z, 14)
    levels = np.unique(rounded)
    means = np.array([np.mean(values[rounded == level], axis=0) for level in levels])
    return levels, means


def read_vtu(path: Path) -> dict[str, np.ndarray]:
    piece = ET.parse(path).getroot().find(".//Piece")
    if piece is None:
        raise ValueError(f"missing Piece in {path}")
    points = data_array(piece.find("Points"), "Coordinates")
    point_data = piece.find("PointData")
    cell_data = piece.find("CellData")
    cells = piece.find("Cells")
    if point_data is None or cell_data is None or cells is None:
        raise ValueError(f"incomplete VTK arrays in {path}")

    z, pw = grouped_mean(points[:, 1], data_array(point_data, "pw"))
    _, sw = grouped_mean(points[:, 1], data_array(point_data, "sw"))
    _, x_salt = grouped_mean(
        points[:, 1], data_array(point_data, "moleFraction")
    )
    try:
        _, w_salt = grouped_mean(
            points[:, 1], data_array(point_data, "massFraction")
        )
    except KeyError:
        w_salt = mass_fraction(x_salt)
    _, rho = grouped_mean(points[:, 1], data_array(point_data, "rhoW"))

    connectivity = data_array(cells, "connectivity").astype(int)
    offsets = data_array(cells, "offsets").astype(int)
    starts = np.r_[0, offsets[:-1]]
    centers = np.array(
        [points[connectivity[start:stop]].mean(axis=0) for start, stop in zip(starts, offsets)]
    )
    velocity = data_array(cell_data, "velocity")
    z_cell, velocity_mean = grouped_mean(centers[:, 1], velocity)
    return {
        "z": z,
        "pw": pw,
        "sw": sw,
        "x": x_salt,
        "w": w_salt,
        "rho": rho,
        "z_cell": z_cell,
        "qx": velocity_mean[:, 0],
        "qz": velocity_mean[:, 1],
    }


def vertical_weights(z: np.ndarray) -> np.ndarray:
    weights = np.empty_like(z)
    weights[0] = 0.5 * (z[1] - z[0])
    weights[-1] = 0.5 * (z[-1] - z[-2])
    weights[1:-1] = 0.5 * (z[2:] - z[:-2])
    return weights


def inventory(profile: dict[str, np.ndarray], width: float) -> dict[str, float]:
    weights = vertical_weights(profile["z"]) * width
    total_density = POROSITY * profile["sw"] * profile["rho"]
    salt_density = total_density * profile["w"]
    return {
        "total_kg": float(np.sum(weights * total_density)),
        "salt_kg": float(np.sum(weights * salt_density)),
        "water_kg": float(np.sum(weights * (total_density - salt_density))),
    }


def relative_l2(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.linalg.norm(a - b) / max(np.linalg.norm(b), 1e-300))


def comparison(
    left: dict[str, np.ndarray],
    right: dict[str, np.ndarray],
) -> dict[str, float]:
    result: dict[str, float] = {}
    for name in ("pw", "sw", "x", "w", "rho"):
        result[f"{name}_relative_l2"] = relative_l2(left[name], right[name])
        result[f"{name}_max_abs"] = float(np.max(np.abs(left[name] - right[name])))
    result["qz_relative_l2"] = relative_l2(left["qz"], right["qz"])
    result["qz_max_abs_m_s"] = float(np.max(np.abs(left["qz"] - right["qz"])))
    return result


def read_budget_max(path: Path) -> dict[str, float] | None:
    budget = path / "component_budget.csv"
    if not budget.exists():
        return None
    with budget.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    return {
        key: max(abs(float(row[key])) for row in rows)
        for key in (
            "total_residual_kg",
            "salt_residual_kg",
            "water_residual_kg",
        )
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--workspace", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    directories = {
        "original_normalized": args.workspace / "run_V5_original_normalized_1d",
        "mass_mode_A": args.workspace / "V5_mass_mode_A_1d",
        "mass_mode_B": args.workspace / "V5_mass_mode_B_1d",
        "mass_mode_A_dt_half": args.workspace / "V7_mass_mode_A_dt_half_1d",
    }
    widths = {name: 0.021428571428571432 for name in directories}
    records: dict[str, dict[float, dict[str, np.ndarray]]] = {}
    inventories: dict[str, dict[float, dict[str, float]]] = {}

    with (args.output / "profiles.csv").open(
        "w", newline="", encoding="utf-8"
    ) as profile_handle, (args.output / "velocity_profiles.csv").open(
        "w", newline="", encoding="utf-8"
    ) as velocity_handle:
        profile_writer = csv.writer(profile_handle)
        velocity_writer = csv.writer(velocity_handle)
        profile_writer.writerow(["run", "time_s", "z_m", "pw_Pa", "Sw", "x_NaCl", "w_NaCl", "rho_kg_m3"])
        velocity_writer.writerow(["run", "time_s", "z_cell_m", "qx_m_s", "qz_m_s"])
        for name, directory in directories.items():
            records[name] = {}
            inventories[name] = {}
            for time, path in read_pvd(directory):
                profile = read_vtu(path)
                records[name][time] = profile
                inventories[name][time] = inventory(profile, widths[name])
                for values in zip(
                    profile["z"],
                    profile["pw"],
                    profile["sw"],
                    profile["x"],
                    profile["w"],
                    profile["rho"],
                ):
                    profile_writer.writerow((name, time, *values))
                for values in zip(
                    profile["z_cell"], profile["qx"], profile["qz"]
                ):
                    velocity_writer.writerow((name, time, *values))

    base_times = sorted(records["mass_mode_A"])
    common = sorted(
        set(base_times)
        & set(records["original_normalized"])
        & set(records["mass_mode_B"])
        & set(records["mass_mode_A_dt_half"])
    )
    comparisons = []
    pairs = {
        "mass_A_vs_original": ("mass_mode_A", "original_normalized"),
        "mass_B_vs_mass_A": ("mass_mode_B", "mass_mode_A"),
        "dt_half_vs_dt": ("mass_mode_A_dt_half", "mass_mode_A"),
    }
    for time in common:
        for label, (left, right) in pairs.items():
            comparisons.append(
                {
                    "pair": label,
                    "time_s": time,
                    **comparison(records[left][time], records[right][time]),
                }
            )
    with (args.output / "comparisons.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=comparisons[0].keys())
        writer.writeheader()
        writer.writerows(comparisons)

    with (args.output / "inventories.csv").open(
        "w", newline="", encoding="utf-8"
    ) as handle:
        writer = csv.writer(handle)
        writer.writerow(["run", "time_s", "total_kg", "salt_kg", "water_kg"])
        for name, by_time in inventories.items():
            for time, values in sorted(by_time.items()):
                writer.writerow(
                    (
                        name,
                        time,
                        values["total_kg"],
                        values["salt_kg"],
                        values["water_kg"],
                    )
                )

    nearest_reference_times = {
        label: min(common, key=lambda value: abs(value - reference))
        for label, reference in REFERENCE_TIMES.items()
    }
    final_time = common[-1]
    final_comparisons = {
        row["pair"]: row
        for row in comparisons
        if row["time_s"] == final_time
    }
    budget_residuals = {
        name: read_budget_max(directory)
        for name, directory in directories.items()
    }
    retry_path = args.workspace / "V7_forced_retry_attempt1" / "accepted_steps.csv"
    with retry_path.open(newline="", encoding="utf-8") as handle:
        retry_rows = list(csv.DictReader(handle))
    retry_evidence = {
        "accepted_steps": len(retry_rows),
        "retry_count": sum(int(row["retried"]) for row in retry_rows),
        "first_requested_dt_s": float(retry_rows[0]["requested_dt_s"]),
        "first_accepted_dt_s": float(retry_rows[0]["accepted_dt_s"]),
    }

    summary = {
        "schema_version": 1,
        "actual_common_output_times_s": common,
        "reference_times_s": REFERENCE_TIMES,
        "nearest_actual_times_s": nearest_reference_times,
        "final_time_s": final_time,
        "final_comparisons": final_comparisons,
        "max_conservative_budget_residuals_kg": budget_residuals,
        "forced_retry": retry_evidence,
        "claims": {
            "numerical_component_conservation": "tested separately from profile agreement",
            "old_solver": "normalized baseline, not physical component conservation",
            "physical_validation": "not established by these numerical comparisons",
        },
    }
    (args.output / "background_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    colors = {
        "original_normalized": "#555555",
        "mass_mode_A": "#0072B2",
        "mass_mode_B": "#D55E00",
        "mass_mode_A_dt_half": "#009E73",
    }
    labels = {
        "original_normalized": "Original normalized",
        "mass_mode_A": "Mass Mode A",
        "mass_mode_B": "Mass Mode B",
        "mass_mode_A_dt_half": "Mass Mode A dt/2",
    }
    fig, axes = plt.subplots(1, 3, figsize=(10.5, 3.4), constrained_layout=True)
    for name in directories:
        profile = records[name][final_time]
        depth_mm = 1000.0 * (0.2 - profile["z"])
        axes[0].plot(profile["x"], depth_mm, color=colors[name], label=labels[name])
        axes[1].plot(profile["sw"], depth_mm, color=colors[name])
        axes[2].plot(profile["rho"], depth_mm, color=colors[name])
    axes[0].set_xlabel("NaCl mole fraction")
    axes[1].set_xlabel("Liquid saturation")
    axes[2].set_xlabel(r"Density (kg m$^{-3}$)")
    for axis in axes:
        axis.set_ylabel("Depth below surface (mm)")
        axis.invert_yaxis()
        axis.grid(alpha=0.25)
    axes[0].legend(frameon=False, fontsize=8)
    fig.savefig(args.output / "BACKGROUND_FINAL_PROFILES.png", dpi=220)
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(9, 3.4), constrained_layout=True)
    for name in directories:
        times = sorted(records[name])
        top_x = [records[name][time]["x"][-1] for time in times]
        total_mass = [inventories[name][time]["total_kg"] for time in times]
        axes[0].plot(np.array(times) / 3600.0, top_x, color=colors[name], label=labels[name])
        axes[1].plot(np.array(times) / 3600.0, total_mass, color=colors[name])
    axes[0].set_xlabel("Time (h)")
    axes[0].set_ylabel("Surface NaCl mole fraction")
    axes[1].set_xlabel("Time (h)")
    axes[1].set_ylabel("Independent total liquid mass (kg)")
    for axis in axes:
        axis.grid(alpha=0.25)
    axes[0].legend(frameon=False, fontsize=8)
    fig.savefig(args.output / "BACKGROUND_TRAJECTORIES.png", dpi=220)
    plt.close(fig)

    fig, axis = plt.subplots(figsize=(6.2, 3.5), constrained_layout=True)
    pair_colors = {
        "mass_A_vs_original": "#7B3294",
        "mass_B_vs_mass_A": "#D55E00",
        "dt_half_vs_dt": "#009E73",
    }
    for label in pairs:
        rows = [row for row in comparisons if row["pair"] == label]
        axis.semilogy(
            [row["time_s"] / 3600.0 for row in rows],
            [max(row["x_relative_l2"], 1e-16) for row in rows],
            color=pair_colors[label],
            label=label.replace("_", " "),
        )
    axis.set_xlabel("Time (h)")
    axis.set_ylabel("NaCl profile relative L2 difference")
    axis.grid(alpha=0.25)
    axis.legend(frameon=False, fontsize=8)
    fig.savefig(args.output / "BACKGROUND_DIFFERENCES.png", dpi=220)
    plt.close(fig)

    print(json.dumps(summary, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
