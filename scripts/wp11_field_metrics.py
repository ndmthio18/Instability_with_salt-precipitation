"""Connectivity-aware field and modal reduction for WP11 trajectories.

VTK point concentrations are mapped through the declared cell connectivity before
they are paired with cell-centre velocity.  The velocity products remain qualitative:
the conservative transport evidence comes only from transport_budget.csv.
"""
import argparse
import csv
import hashlib
import json
import xml.etree.ElementTree as ET
from pathlib import Path

import numpy as np

from single_wavelength_post import widths
from wp6_4t_metrics import neumann_spectrum

T_S = 545692.2550893936 / 4.0
TRACKED_MODES = (1, 2, 4, 6, 10)


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle))


def parse_vtu(path):
    root = ET.fromstring(Path(path).read_bytes())
    arrays = {}
    metadata = {}
    for element in root.iter("DataArray"):
        name = element.attrib.get("Name")
        if name:
            arrays[name] = np.fromstring(element.text or "", sep=" ")
            metadata[name] = element.attrib
    points = arrays["Coordinates"].reshape(-1, 3)
    offsets = arrays["offsets"].astype(int)
    connectivity_flat = arrays["connectivity"].astype(int)
    starts = np.r_[0, offsets[:-1]]
    cell_sizes = offsets - starts
    if not np.all(cell_sizes == 4):
        raise AssertionError(f"non-quadrilateral cells in {path}: {np.unique(cell_sizes)}")
    connectivity = np.vstack([connectivity_flat[a:b] for a, b in zip(starts, offsets)])
    if connectivity.min() < 0 or connectivity.max() >= len(points):
        raise AssertionError(f"invalid connectivity index in {path}")
    velocity_name = next(name for name in arrays if "velocity" in name.lower())
    components = int(metadata[velocity_name].get("NumberOfComponents", 3))
    velocity = arrays[velocity_name].reshape(-1, components)
    if len(velocity) != len(connectivity):
        raise AssertionError(f"cell velocity/connectivity length mismatch in {path}")
    point_x = arrays["moleFraction"]
    cell_x = point_x[connectivity].mean(axis=1)
    centers = points[connectivity].mean(axis=1)
    return points, point_x, connectivity, centers, cell_x, velocity, velocity_name


def point_grid(points, values):
    xs, zs = np.unique(points[:, 0]), np.unique(points[:, 1])
    ix, iz = np.searchsorted(xs, points[:, 0]), np.searchsorted(zs, points[:, 1])
    linear = iz * len(xs) + ix
    if len(np.unique(linear)) != len(points) or len(points) != len(xs) * len(zs):
        raise AssertionError("point coordinates do not form a bijective structured grid")
    grid = np.empty((len(zs), len(xs)))
    grid[iz, ix] = values
    return xs, zs, grid


def cell_grid(centers, values):
    xs, zs = np.unique(centers[:, 0]), np.unique(centers[:, 1])
    ix, iz = np.searchsorted(xs, centers[:, 0]), np.searchsorted(zs, centers[:, 1])
    linear = iz * len(xs) + ix
    if len(np.unique(linear)) != len(centers) or len(centers) != len(xs) * len(zs):
        raise AssertionError("cell centres do not form a bijective structured grid")
    grid = np.empty((len(zs), len(xs)))
    grid[iz, ix] = values
    return xs, zs, grid


def write_csv(path, records):
    if not records:
        return
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def first_reappearance_after_absence(records):
    seen_frame_without_downflow = False
    for row in records:
        if row["peak_downward_velocity_m_s"] <= 0:
            seen_frame_without_downflow = True
        elif seen_frame_without_downflow:
            return row
    return None


def main(run, out, label):
    run, out = Path(run), Path(out)
    out.mkdir(parents=True, exist_ok=True)
    fields = sorted(run.glob("*.vtu"))
    diagnostics = read_csv(run / "post" / "weighted_diagnostics.csv")
    manifest = json.loads((run / "post" / "field_hashes.json").read_text())
    if not (fields and len(fields) == len(diagnostics) == len(manifest)):
        raise AssertionError("field/diagnostic/hash-manifest count mismatch")

    target_T = [2, 3, 4, 5, 6, 10, 15, 20, float(diagnostics[-1]["time_s"]) / T_S]
    target_T = sorted(set(q for q in target_T if q <= float(diagnostics[-1]["time_s"]) / T_S + 1e-9))
    frame_times = np.array([float(r["time_s"]) / T_S for r in diagnostics])
    selected = {int(np.argmin(np.abs(frame_times - q))): q for q in target_T}

    mode_rows, velocity_rows, panel_arrays = [], [], {}
    connectivity_reference = None
    connectivity_consistent = True
    geometry_consistent = True
    reference_points = None
    velocity_name = None
    for frame, path in enumerate(fields):
        expected_path = Path(manifest[frame]["path"]).name
        if path.name != expected_path or sha256(path) != manifest[frame]["sha256"]:
            raise AssertionError(f"field hash/path mismatch at frame {frame}")
        points, point_x, connectivity, centers, cell_x, velocity, this_velocity_name = parse_vtu(path)
        velocity_name = velocity_name or this_velocity_name
        if this_velocity_name != velocity_name:
            raise AssertionError("velocity array name changed across frames")
        if connectivity_reference is None:
            connectivity_reference = connectivity.copy()
            reference_points = points.copy()
        else:
            connectivity_consistent &= np.array_equal(connectivity, connectivity_reference)
            geometry_consistent &= np.array_equal(points, reference_points)
        xs, zs, x_grid = point_grid(points, point_x)
        x_weights = widths(xs)
        anomaly = x_grid - (x_grid @ x_weights / x_weights.sum())[:, None]
        _, energy = neumann_spectrum(xs, zs, anomaly)
        total_energy = float(energy.sum())
        dominant = int(np.argmax(energy) + 1)
        for mode, value in enumerate(energy, 1):
            mode_rows.append({
                "grid": label,
                "frame": frame,
                "step": int(diagnostics[frame]["step"]),
                "time_s": float(diagnostics[frame]["time_s"]),
                "time_over_T": frame_times[frame],
                "mode": mode,
                "physical_wavenumber_m_inv": mode * np.pi / x_weights.sum(),
                "wavelength_m": 2.0 * x_weights.sum() / mode,
                "cells_per_wavelength": 2.0 * (len(xs) - 1) / mode,
                "energy": float(value),
                "energy_fraction": float(value / total_energy) if total_energy > 0 else 0.0,
                "dominant_mode": dominant,
                "tracked_mode": int(mode in TRACKED_MODES),
            })

        cxs, czs, x_cell_grid = cell_grid(centers, cell_x)
        _, _, down_grid = cell_grid(centers, -velocity[:, 1])
        row_mean = x_cell_grid.mean(axis=1, keepdims=True)
        x_anomaly = x_cell_grid - row_mean
        down_anomaly = down_grid - down_grid.mean(axis=1, keepdims=True)
        positive_anomaly = x_anomaly[x_anomaly > 0]
        anomaly_cut = float(np.quantile(positive_anomaly, 0.9)) if positive_anomaly.size else np.inf
        high_finger = x_anomaly >= anomaly_cut
        overlap = high_finger & (down_grid > 0)
        penetration_mm = float((zs[-1] - czs[overlap.any(axis=1)].min()) * 1000.0) if overlap.any() else None
        flat_x, flat_down = x_cell_grid.ravel(), down_grid.ravel()
        high_global = flat_x >= np.quantile(flat_x, 0.9)
        velocity_rows.append({
            "grid": label,
            "nominal_time_over_T": selected.get(frame, ""),
            "frame": frame,
            "step": int(diagnostics[frame]["step"]),
            "time_s": float(diagnostics[frame]["time_s"]),
            "time_over_T": frame_times[frame],
            "downward_cell_fraction": float(np.mean(down_grid > 0)),
            "peak_downward_velocity_m_s": float(down_grid.max()),
            "peak_upward_velocity_m_s": float(-down_grid.min()),
            "salt_down_correlation": float(np.corrcoef(flat_x, flat_down)[0, 1]),
            "horizontal_anomaly_down_correlation": float(np.corrcoef(x_anomaly.ravel(), down_anomaly.ravel())[0, 1]),
            "high_concentration_downward_fraction": float(np.mean(flat_down[high_global] > 0)),
            "high_concentration_mean_downward_velocity_m_s": float(flat_down[high_global].mean()),
            "finger_penetration_depth_mm": penetration_mm,
            "finger_definition": "deepest cell with downward velocity and horizontal concentration anomaly in the top decile of positive anomalies",
            "source_sha256": manifest[frame]["sha256"],
        })
        if frame in selected:
            prefix = f"frame_{frame}"
            panel_arrays[f"{prefix}_point_x_m"] = xs
            panel_arrays[f"{prefix}_point_z_m"] = zs
            panel_arrays[f"{prefix}_concentration"] = x_grid
            panel_arrays[f"{prefix}_cell_x_m"] = cxs
            panel_arrays[f"{prefix}_cell_z_m"] = czs
            panel_arrays[f"{prefix}_downward_velocity_m_s"] = down_grid
            panel_arrays[f"{prefix}_actual_time_over_T"] = np.array(frame_times[frame])
            panel_arrays[f"{prefix}_nominal_time_over_T"] = np.array(selected[frame])

    tracked = [row for row in mode_rows if row["tracked_mode"]]
    write_csv(out / "mode_spectrum.csv", mode_rows)
    write_csv(out / "tracked_mode_series.csv", tracked)
    write_csv(out / "velocity_checkpoints.csv", velocity_rows)
    np.savez_compressed(out / "field_panels.npz", **panel_arrays)
    result = {
        "schema_version": 1,
        "grid": label,
        "frames": len(fields),
        "first_step": int(diagnostics[0]["step"]),
        "last_step": int(diagnostics[-1]["step"]),
        "last_time_over_T": float(frame_times[-1]),
        "point_count": int(len(reference_points)),
        "cell_count": int(len(connectivity_reference)),
        "cell_connectivity": {
            "quad_cell_size_verified": True,
            "indices_in_range_verified": True,
            "connectivity_identical_across_frames": bool(connectivity_consistent),
            "geometry_identical_across_frames": bool(geometry_consistent),
            "point_to_cell_concentration_mapping": "arithmetic mean of the four point values declared by each VTK connectivity row",
            "velocity_mapping": "cell-data row paired with the same VTK connectivity row, then bijectively mapped by geometric cell centre",
        },
        "velocity_array_name": velocity_name,
        "velocity_scope": "qualitative cell-centre Darcy reconstruction only; not a conservative face flux",
        "mode_basis": "Neumann cosine modes; physical wavenumber k_m=m*pi/W",
        "tracked_modes": list(TRACKED_MODES),
        "first_observed_downward_flow_sign_only": next((row for row in velocity_rows if row["peak_downward_velocity_m_s"] > 0), None),
        "first_reappearing_downward_flow_after_absence": first_reappearance_after_absence(velocity_rows),
        "downward_onset_semantics": "sign-only includes the imposed initial eigenprofile; mechanism timing uses the first positive saved frame after at least one saved frame with no downward cell",
        "selected_velocity_checkpoints": [row for row in velocity_rows if row["nominal_time_over_T"] != ""],
        "hashes": {
            "field_manifest": sha256(run / "post" / "field_hashes.json"),
            "weighted_diagnostics": sha256(run / "post" / "weighted_diagnostics.csv"),
        },
    }
    (out / "field_metrics.json").write_text(json.dumps(result, indent=2))
    print(json.dumps({"grid": label, "frames": len(fields), "last_T": frame_times[-1], "selected": len(selected)}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("run")
    parser.add_argument("out")
    parser.add_argument("label")
    args = parser.parse_args()
    main(args.run, args.out, args.label)
