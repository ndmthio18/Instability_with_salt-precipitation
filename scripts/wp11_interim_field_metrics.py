"""Reduce the immutable sparse VTU set from a WP11 interim snapshot.

This intentionally does not invent unsaved frames or a terminal solver event.  Modal
and reconstructed-velocity histories are sampled only at the fields listed in the
snapshot manifest; conservative transport remains sourced from transport_budget.csv.
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np

from single_wavelength_post import widths
from wp6_4t_metrics import neumann_spectrum
from wp11_field_metrics import TRACKED_MODES, T_S, cell_grid, parse_vtu, point_grid, sha256


def write_csv(path, records):
    if not records:
        return
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def main(run, out, label):
    run, out = Path(run), Path(out)
    out.mkdir(parents=True, exist_ok=True)
    manifest_path = run / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    declared = {row["path"]: row for row in manifest["files"]}

    mode_rows, velocity_rows, percentile_rows, panel_arrays = [], [], [], {}
    reference_connectivity = None
    reference_points = None
    velocity_name = None
    for index, selected in enumerate(manifest["selected_fields"]):
        path = run / selected["filename"]
        if not path.is_file() or sha256(path) != selected["sha256"]:
            raise AssertionError(f"missing or hash-mismatched field: {path}")
        if declared[path.name]["sha256"] != selected["sha256"]:
            raise AssertionError(f"manifest disagreement for {path.name}")
        points, point_x, connectivity, centers, cell_x, velocity, this_velocity_name = parse_vtu(path)
        if reference_connectivity is None:
            reference_connectivity = connectivity.copy()
            reference_points = points.copy()
            velocity_name = this_velocity_name
        else:
            if not np.array_equal(connectivity, reference_connectivity):
                raise AssertionError("connectivity changed across sparse snapshot fields")
            if not np.array_equal(points, reference_points):
                raise AssertionError("geometry changed across sparse snapshot fields")
            if this_velocity_name != velocity_name:
                raise AssertionError("velocity array name changed")

        time_T = float(selected["actual_time_over_T"])
        frame = int(selected["frame"])
        xs, zs, x_grid = point_grid(points, point_x)
        x_weights = widths(xs)
        anomaly = x_grid - (x_grid @ x_weights / x_weights.sum())[:, None]
        _, energy = neumann_spectrum(xs, zs, anomaly)
        total_energy = float(energy.sum())
        dominant = int(np.argmax(energy) + 1)
        odd_energy = float(energy[0::2].sum())
        for mode, value in enumerate(energy, 1):
            mode_rows.append({
                "grid": label, "frame": frame, "step": int(selected["step"]),
                "time_s": float(selected["actual_time_s"]), "time_over_T": time_T,
                "mode": mode, "physical_wavenumber_m_inv": mode * np.pi / x_weights.sum(),
                "wavelength_m": 2.0 * x_weights.sum() / mode,
                "cells_per_wavelength": 2.0 * (len(xs) - 1) / mode,
                "energy": float(value),
                "energy_fraction": float(value / total_energy) if total_energy else 0.0,
                "dominant_mode": dominant, "tracked_mode": int(mode in TRACKED_MODES),
                "odd_mode_energy_fraction": odd_energy / total_energy if total_energy else 0.0,
                "sampling_scope": "sparse immutable interim VTU snapshot",
            })

        cxs, czs, x_cell_grid = cell_grid(centers, cell_x)
        _, _, down_grid = cell_grid(centers, -velocity[:, 1])
        row_mean = x_cell_grid.mean(axis=1, keepdims=True)
        x_anomaly = x_cell_grid - row_mean
        down_anomaly = down_grid - down_grid.mean(axis=1, keepdims=True)
        positive_anomaly = x_anomaly[x_anomaly > 0]
        anomaly_cut = float(np.quantile(positive_anomaly, 0.9)) if positive_anomaly.size else np.inf
        overlap = (x_anomaly >= anomaly_cut) & (down_grid > 0)
        penetration_mm = float((zs[-1] - czs[overlap.any(axis=1)].min()) * 1000.0) if overlap.any() else None
        flat_x, flat_down = x_cell_grid.ravel(), down_grid.ravel()
        high_global = flat_x >= np.quantile(flat_x, 0.9)
        velocity_rows.append({
            "grid": label, "nominal_time_over_T": selected["requested_time_over_T"],
            "frame": frame, "step": int(selected["step"]),
            "time_s": float(selected["actual_time_s"]), "time_over_T": time_T,
            "downward_cell_fraction": float(np.mean(down_grid > 0)),
            "peak_downward_velocity_m_s": float(down_grid.max()),
            "peak_upward_velocity_m_s": float(-down_grid.min()),
            "salt_down_correlation": float(np.corrcoef(flat_x, flat_down)[0, 1]),
            "horizontal_anomaly_down_correlation": float(np.corrcoef(x_anomaly.ravel(), down_anomaly.ravel())[0, 1]),
            "high_concentration_downward_fraction": float(np.mean(flat_down[high_global] > 0)),
            "high_concentration_mean_downward_velocity_m_s": float(flat_down[high_global].mean()),
            "finger_penetration_depth_mm": penetration_mm,
            "finger_definition": "deepest cell with downward velocity and horizontal concentration anomaly in the top decile of positive anomalies",
            "source_sha256": selected["sha256"],
            "sampling_scope": "sparse immutable interim VTU snapshot",
        })
        top = x_grid[-1]
        percentile_rows.append({
            "grid": label, "frame": frame, "step": int(selected["step"]),
            "time_over_T": time_T, "p05": float(np.quantile(top, 0.05)),
            "p25": float(np.quantile(top, 0.25)), "p50": float(np.quantile(top, 0.50)),
            "p75": float(np.quantile(top, 0.75)), "p95": float(np.quantile(top, 0.95)),
            "minimum": float(top.min()), "maximum": float(top.max()),
        })
        prefix = f"frame_{frame}"
        panel_arrays[f"{prefix}_point_x_m"] = xs
        panel_arrays[f"{prefix}_point_z_m"] = zs
        panel_arrays[f"{prefix}_concentration"] = x_grid
        panel_arrays[f"{prefix}_cell_x_m"] = cxs
        panel_arrays[f"{prefix}_cell_z_m"] = czs
        panel_arrays[f"{prefix}_downward_velocity_m_s"] = down_grid
        panel_arrays[f"{prefix}_actual_time_over_T"] = np.array(time_T)
        panel_arrays[f"{prefix}_nominal_time_over_T"] = np.array(selected["requested_time_over_T"] if selected["requested_time_over_T"] != "latest_stable" else time_T)

    write_csv(out / "mode_spectrum.csv", mode_rows)
    write_csv(out / "tracked_mode_series.csv", [r for r in mode_rows if r["tracked_mode"]])
    write_csv(out / "velocity_checkpoints.csv", velocity_rows)
    write_csv(out / "surface_percentiles.csv", percentile_rows)
    np.savez_compressed(out / "field_panels.npz", **panel_arrays)
    result = {
        "schema_version": 1, "grid": label, "frames": len(manifest["selected_fields"]),
        "first_step": int(manifest["selected_fields"][0]["step"]),
        "last_step": int(manifest["selected_fields"][-1]["step"]),
        "last_time_over_T": float(manifest["selected_fields"][-1]["actual_time_over_T"]),
        "snapshot_boundary_time_over_T": float(manifest["snapshot_end_time_over_T"]),
        "point_count": int(len(reference_points)), "cell_count": int(len(reference_connectivity)),
        "cell_connectivity": {"quad_cell_size_verified": True, "indices_in_range_verified": True,
            "connectivity_identical_across_frames": True, "geometry_identical_across_frames": True,
            "point_to_cell_concentration_mapping": "arithmetic mean of the four point values declared by each VTK connectivity row",
            "velocity_mapping": "cell-data row paired with the same VTK connectivity row, then bijectively mapped by geometric cell centre"},
        "velocity_array_name": velocity_name,
        "velocity_scope": "qualitative cell-centre Darcy reconstruction only; not a conservative face flux",
        "mode_basis": "Neumann cosine modes; physical wavenumber k_m=m*pi/W",
        "tracked_modes": list(TRACKED_MODES),
        "sampling_scope": "15 complete immutable fields; no interpolation is used to create missing field states",
        "selected_velocity_checkpoints": velocity_rows,
        "hashes": {"snapshot_manifest": sha256(manifest_path)},
    }
    (out / "field_metrics.json").write_text(json.dumps(result, indent=2))
    print(json.dumps({"grid": label, "frames": result["frames"], "last_field_T": result["last_time_over_T"]}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("run")
    parser.add_argument("out")
    parser.add_argument("label")
    args = parser.parse_args()
    main(args.run, args.out, args.label)
