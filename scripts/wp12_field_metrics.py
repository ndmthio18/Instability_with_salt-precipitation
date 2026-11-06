"""Reduce WP12 fields to modal, symmetry and resolved downflow time series."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

from single_wavelength_post import widths
from wp11_field_metrics import parse_vtu, point_grid, cell_grid, sha256
from wp6_4t_metrics import neumann_spectrum


T_S = 545692.2550893936 / 4.0
TRACKED = (2, 4, 6, 10)


def read_csv(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def weighted_norm2(values: np.ndarray, xs: np.ndarray, zs: np.ndarray) -> float:
    weights = np.outer(widths(zs), widths(xs))
    return float(np.sum(values * values * weights))


def first_two_consecutive(rows: list[dict], key: str):
    for left, right in zip(rows, rows[1:]):
        if left[key] and right[key]:
            return left
    return None


def main(run: Path, out: Path, label: str, max_T: float | None = None, post_dir: Path | None = None) -> None:
    out.mkdir(parents=True, exist_ok=True)
    fields = sorted(run.glob("*.vtu"))
    post_dir = post_dir or run / "post"
    diagnostics = read_csv(post_dir / "weighted_diagnostics.csv")
    manifest = json.loads((post_dir / "field_hashes.json").read_text())
    if not fields or len(fields) != len(diagnostics) or len(fields) != len(manifest):
        raise AssertionError(f"field mapping mismatch: {len(fields)}, {len(diagnostics)}, {len(manifest)}")
    if max_T is not None:
        selected = [index for index, row in enumerate(diagnostics) if float(row["time_s"]) / T_S <= max_T + 1e-8]
        fields = [fields[index] for index in selected]
        diagnostics = [diagnostics[index] for index in selected]
        manifest = [manifest[index] for index in selected]

    rows = []
    for frame, (path, diagnostic, recorded) in enumerate(zip(fields, diagnostics, manifest)):
        if path.name != Path(recorded["path"]).name or sha256(path) != recorded["sha256"]:
            raise AssertionError(f"field hash/path mismatch at frame {frame}")
        points, point_x, _, centers, _, velocity, _ = parse_vtu(path)
        xs, zs, grid = point_grid(points, point_x)
        horizontal_mean = grid @ widths(xs) / widths(xs).sum()
        anomaly = grid - horizontal_mean[:, None]
        _, energy = neumann_spectrum(xs, zs, anomaly)
        total = float(energy.sum())
        odd = float(sum(value for mode, value in enumerate(energy, 1) if mode % 2 == 1))
        mirror_num = weighted_norm2(grid - grid[:, ::-1], xs, zs)
        mirror_den = weighted_norm2(anomaly, xs, zs)
        _, _, down = cell_grid(centers, -velocity[:, 1])
        peak_down = float(down.max())
        resolved_fraction = float(np.mean(down > 1e-10))
        row = {
            "grid": label,
            "frame": frame,
            "step": int(diagnostic["step"]),
            "time_s": float(diagnostic["time_s"]),
            "time_over_T": float(diagnostic["time_s"]) / T_S,
            "dominant_mode": int(np.argmax(energy) + 1) if total > 0 else 0,
            "total_nonzero_mode_energy": total,
            "odd_mode_energy_fraction": odd / total if total > 0 else 0.0,
            "mirror_symmetry_error": np.sqrt(mirror_num / mirror_den) if mirror_den > 0 else 0.0,
            "peak_downward_velocity_m_s": peak_down,
            "resolved_downward_cell_fraction": resolved_fraction,
            "resolved_downflow_frame": int(peak_down > 1e-10 and resolved_fraction > 0.001),
            "field_sha256": recorded["sha256"],
        }
        for mode in TRACKED:
            row[f"mode_{mode}_energy"] = float(energy[mode - 1]) if len(energy) >= mode else 0.0
            row[f"mode_{mode}_energy_fraction"] = float(energy[mode - 1] / total) if total > 0 and len(energy) >= mode else 0.0
        rows.append(row)

    write_csv(out / "field_state_series.csv", rows)
    onset = first_two_consecutive(rows, "resolved_downflow_frame")
    summary = {
        "schema_version": 1,
        "grid": label,
        "frames": len(rows),
        "first_step": rows[0]["step"],
        "last_step": rows[-1]["step"],
        "last_time_over_T": rows[-1]["time_over_T"],
        "resolved_downflow_definition": "first of two consecutive saved fields with peak downward velocity >1e-10 m/s and downward cell fraction >0.001",
        "resolved_downflow_onset": onset,
        "mirror_error_definition": "area-weighted L2[x(x,z)-x(W-x,z)] / area-weighted L2[x-horizontal_mean(z)]",
        "odd_mode_definition": "odd nonzero Neumann-cosine energy divided by all nonzero Neumann-cosine energy",
        "tracked_modes": list(TRACKED),
        "maximum_mirror_symmetry_error": max(row["mirror_symmetry_error"] for row in rows),
        "maximum_odd_mode_energy_fraction": max(row["odd_mode_energy_fraction"] for row in rows),
    }
    (out / "wp12_field_metrics.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps({"grid": label, "frames": len(rows), "resolved_downflow_T": None if onset is None else onset["time_over_T"]}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("run", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("label")
    parser.add_argument("--max-T", type=float)
    parser.add_argument("--post-dir", type=Path)
    args = parser.parse_args()
    main(args.run, args.out, args.label, args.max_T, args.post_dir)
