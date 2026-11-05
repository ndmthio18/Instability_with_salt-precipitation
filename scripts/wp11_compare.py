"""Assemble auditable Stage-A/Stage-B comparison tables from WP11 reductions."""
import argparse
import csv
import json
from collections import Counter
from pathlib import Path


def read_csv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path, records):
    if not records:
        return
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]))
        writer.writeheader()
        writer.writerows(records)


def nearest(records, target, key="time_over_T"):
    return min(records, key=lambda row: abs(float(row[key]) - target))


def load(label, path):
    path = Path(path)
    return {
        "label": label,
        "path": path,
        "scalar": json.loads((path / "scalar_metrics.json").read_text()),
        "field": json.loads((path / "field_metrics.json").read_text()),
        "surface": read_csv(path / "surface_series.csv"),
        "flux": read_csv(path / "flux_series.csv"),
        "modes": read_csv(path / "mode_spectrum.csv"),
        "velocity": read_csv(path / "velocity_checkpoints.csv"),
    }


def main(out, specifications):
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    grids = [load(label, path) for label, path in specifications]
    shared_T = min(grid["scalar"]["stop"]["time_over_T"] for grid in grids)
    surface_summary, flux_summary, velocity_summary, mode_windows, common_time = [], [], [], [], []
    for grid in grids:
        label = grid["label"]
        surface = [row for row in grid["surface"] if float(row["time_over_T"]) <= shared_T + 1e-9]
        peak = max(surface, key=lambda row: float(row["surface_mean"]))
        terminal = surface[-1]
        stop = grid["scalar"]["stop"]
        surface_summary.append({
            "grid": label,
            "shared_end_T": shared_T,
            "surface_mean_peak": float(peak["surface_mean"]),
            "surface_mean_peak_T": float(peak["time_over_T"]),
            "surface_mean_at_shared_end": float(terminal["surface_mean"]),
            "decline_from_peak": float(peak["surface_mean"]) - float(terminal["surface_mean"]),
            "decline_fraction_from_peak": (float(peak["surface_mean"]) - float(terminal["surface_mean"])) / float(peak["surface_mean"]),
            "surface_max_at_shared_end": float(terminal["surface_max"]),
            "domain_max_at_shared_end": float(terminal["domain_max"]),
            "stop_reason": stop["reason"],
            "stop_T": stop["time_over_T"],
            "threshold_low_T": None if stop.get("global_threshold_low_s") is None else stop["global_threshold_low_s"] / (545692.2550893936 / 4.0),
            "threshold_high_T": None if stop.get("global_threshold_high_s") is None else stop["global_threshold_high_s"] / (545692.2550893936 / 4.0),
        })
        for depth, plane in grid["scalar"]["planes"].items():
            for component in ("advective", "total"):
                summary = plane[component]["above_1e-15"]
                longest = summary["longest_continuous_positive_interval"]
                flux_summary.append({
                    "grid": label,
                    "depth_mm": float(depth),
                    "component": component,
                    "threshold_m_s": 1e-15,
                    "first_positive_step": summary["first_positive_step"],
                    "first_positive_T": summary.get("first_positive_time_over_T"),
                    "last_positive_step": summary["last_positive_step"],
                    "last_positive_T": summary.get("last_positive_time_over_T"),
                    "positive_step_count": summary["positive_step_count"],
                    "longest_continuous_first_T": None if longest is None else longest["first_time_over_T"],
                    "longest_continuous_last_T": None if longest is None else longest["last_time_over_T"],
                    "longest_continuous_step_count": 0 if longest is None else longest["accepted_step_count"],
                    "longest_continuous_duration_T": 0.0 if longest is None else longest["accepted_step_support_duration_over_T"],
                    "peak_positive_flux_m_s": summary["peak_positive_flux_m_s"],
                    "positive_step_integral_m": summary["positive_step_integral_m"],
                })
        velocity = grid["velocity"]
        peak_down = max(velocity, key=lambda row: float(row["peak_downward_velocity_m_s"]))
        max_penetration_rows = [row for row in velocity if row["finger_penetration_depth_mm"] not in ("", None)]
        max_penetration = max(max_penetration_rows, key=lambda row: float(row["finger_penetration_depth_mm"])) if max_penetration_rows else None
        reappearance = grid["field"]["first_reappearing_downward_flow_after_absence"]
        velocity_summary.append({
            "grid": label,
            "sign_only_first_T": grid["field"]["first_observed_downward_flow_sign_only"]["time_over_T"],
            "reappearance_after_absence_T": None if reappearance is None else reappearance["time_over_T"],
            "reappearance_step": None if reappearance is None else reappearance["step"],
            "reappearance_peak_down_m_s": None if reappearance is None else reappearance["peak_downward_velocity_m_s"],
            "overall_peak_down_m_s": float(peak_down["peak_downward_velocity_m_s"]),
            "overall_peak_down_T": float(peak_down["time_over_T"]),
            "area_fraction_at_peak": float(peak_down["downward_cell_fraction"]),
            "salt_down_correlation_at_peak": float(peak_down["salt_down_correlation"]),
            "maximum_finger_penetration_mm": None if max_penetration is None else float(max_penetration["finger_penetration_depth_mm"]),
            "maximum_finger_penetration_T": None if max_penetration is None else float(max_penetration["time_over_T"]),
        })
        one_mode_per_frame = [row for row in grid["modes"] if int(row["mode"]) == 1]
        bounds = ((0, 4), (4, 6), (6, 10), (10, shared_T + 1e-9))
        for start, end in bounds:
            values = [int(row["dominant_mode"]) for row in one_mode_per_frame if start <= float(row["time_over_T"]) < end]
            counts = Counter(values)
            mode_windows.append({
                "grid": label,
                "start_T": start,
                "end_T": min(end, shared_T),
                "frame_count": len(values),
                "dominant_mode_counts": json.dumps(dict(sorted(counts.items()))),
                "most_frequent_dominant_mode": None if not values else counts.most_common(1)[0][0],
            })
        for target in (4.0, 5.0, 6.0, shared_T):
            if target > shared_T + 1e-9:
                continue
            srow = nearest(grid["surface"], target)
            vrow = nearest(grid["velocity"], target)
            common_time.append({
                "grid": label,
                "nominal_T": target,
                "surface_actual_T": float(srow["time_over_T"]),
                "surface_mean": float(srow["surface_mean"]),
                "surface_max": float(srow["surface_max"]),
                "velocity_actual_T": float(vrow["time_over_T"]),
                "downward_area_fraction": float(vrow["downward_cell_fraction"]),
                "peak_downward_velocity_m_s": float(vrow["peak_downward_velocity_m_s"]),
                "salt_down_correlation": float(vrow["salt_down_correlation"]),
                "finger_penetration_depth_mm": None if vrow["finger_penetration_depth_mm"] == "" else float(vrow["finger_penetration_depth_mm"]),
            })
    write_csv(out / "surface_summary.csv", surface_summary)
    write_csv(out / "positive_flux_summary.csv", flux_summary)
    write_csv(out / "velocity_summary.csv", velocity_summary)
    write_csv(out / "mode_window_counts.csv", mode_windows)
    write_csv(out / "common_time_checkpoints.csv", common_time)
    result = {
        "schema_version": 1,
        "grids": [grid["label"] for grid in grids],
        "shared_comparison_end_T": shared_T,
        "surface": surface_summary,
        "velocity": velocity_summary,
        "positive_flux_definition": "accepted steps with signed conservative plane-integrated flux >1e-15 m/s; first-last span is not treated as continuous",
        "velocity_onset_definition": "first positive saved cell-centre downward velocity after at least one saved frame with no downward cell; t=0 sign-only seed is retained separately",
        "scope": "comparison tables only; outcome A/B/C requires evidence-bounded scientific judgment",
    }
    (out / "comparison_summary.json").write_text(json.dumps(result, indent=2))
    print(json.dumps({"grids": result["grids"], "shared_T": shared_T}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("out")
    parser.add_argument("grids", nargs="+")
    args = parser.parse_args()
    specs = [specification.split("=", 1) for specification in args.grids]
    main(args.out, specs)
