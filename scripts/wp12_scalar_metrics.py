"""Compute preregistered WP12 scalar and conservative-flux metrics."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


T_S = 545692.2550893936 / 4.0
DEPTHS = {
    "upper_0.995000": 1.0,
    "upper_0.990000": 2.0,
    "upper_0.975000": 5.0,
    "upper_0.950000": 10.0,
}


def read_csv(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def sustained(records, predicate, length=50):
    run = []
    previous = None
    for row in records:
        step = int(row["step"])
        if predicate(row):
            if previous is None or step == previous + 1:
                run.append(row)
            else:
                run = [row]
            if len(run) >= length:
                return run[0]
            previous = step
        else:
            run = []
            previous = None
    return None


def flux_metrics(records, field: str):
    accepted = [row for row in records if int(row["step"]) > 0]
    values = [float(row[field]) for row in accepted]
    dts = [float(row["dt_s"]) for row in accepted]
    first = sustained(accepted, lambda row: float(row[field]) > 1e-15)
    peak_index = max(range(len(values)), key=values.__getitem__)
    return {
        "first_sustained_above_1e-15_step": None if first is None else int(first["step"]),
        "first_sustained_above_1e-15_time_s": None if first is None else float(first["time_s"]),
        "first_sustained_above_1e-15_time_over_T": None if first is None else float(first["time_s"]) / T_S,
        "peak_positive_m_s": max(0.0, values[peak_index]),
        "peak_step": int(accepted[peak_index]["step"]),
        "peak_time_over_T": float(accepted[peak_index]["time_s"]) / T_S,
        "integrated_positive_m": sum(max(value, 0.0) * dt for value, dt in zip(values, dts)),
        "integrated_resolved_positive_m": sum(value * dt for value, dt in zip(values, dts) if value > 1e-15),
        "integrated_signed_m": sum(value * dt for value, dt in zip(values, dts)),
        "accepted_steps": len(accepted),
    }


def main(run: Path, out: Path, label: str, amplitude: float, max_T: float | None = None, post_dir: Path | None = None) -> None:
    out.mkdir(parents=True, exist_ok=True)
    surface = read_csv(run / "surface_events.csv")
    online = read_csv(run / "weighted_online.csv")
    budget = read_csv(run / "transport_budget.csv")
    source_stop = json.loads((run / "stop_record.json").read_text())
    if max_T is not None:
        surface = [row for row in surface if float(row["time_s"]) / T_S <= max_T + 1e-8]
        online = [row for row in online if float(row["time_s"]) / T_S <= max_T + 1e-8]
        budget = [row for row in budget if float(row["time_s"]) / T_S <= max_T + 1e-8]
        stop = {
            "reason": "registered_analysis_boundary",
            "last_accepted_time_s": float(surface[-1]["time_s"]),
            "accepted_steps": int(surface[-1]["step"]),
            "requested_analysis_end_T": max_T,
            "source_stop": source_stop,
        }
    else:
        stop = source_stop
    if len(surface) != len(online):
        raise AssertionError("surface/online row mismatch")

    harmonic_rows = []
    for row in online:
        h2, h3 = float(row["harmonic2_top"]), float(row["harmonic3_top"])
        harmonic_rows.append(dict(row, H23=math.hypot(h2, h3)))
    milestone = sustained(harmonic_rows[1:], lambda row: float(row["H23"]) >= 1e-5)

    peak = max(surface, key=lambda row: float(row["x_top_area_mean"]))
    terminal = surface[-1]
    planes = {}
    for region, depth in DEPTHS.items():
        records = [row for row in budget if row["region"] == region]
        if len(records) != len(surface):
            raise AssertionError(f"{region}: budget/surface row mismatch")
        advective = flux_metrics(records, "cut_salt_adv_down_m_s")
        total = flux_metrics(records, "cut_salt_total_down_m_s")
        final = records[-1]
        advective["reported_cumulative_m"] = float(final["cumulative_cut_adv_m"])
        advective["reintegration_error_m"] = advective["integrated_signed_m"] - advective["reported_cumulative_m"]
        total["reported_cumulative_m"] = float(final["cumulative_cut_salt_m"])
        total["reintegration_error_m"] = total["integrated_signed_m"] - total["reported_cumulative_m"]
        planes[str(depth)] = {
            "region": region,
            "actual_depth_mm": (1.0 - float(final["actual_z_over_H"])) * 200.0,
            "advective": advective,
            "total": total,
            "maximum_abs_salt_residual_m": max(abs(float(row["salt_residual_m"])) for row in records),
            "maximum_abs_water_residual_m": max(abs(float(row["water_residual_m"])) for row in records),
            "maximum_abs_whole_salt_residual_m": max(abs(float(row["whole_salt_residual_m"])) for row in records),
            "maximum_abs_whole_water_residual_m": max(abs(float(row["whole_water_residual_m"])) for row in records),
        }

    post_dir = post_dir or run / "post"
    independent = json.loads((post_dir / "independent_checks.json").read_text())
    result = {
        "schema_version": 1,
        "grid": label,
        "amplitude": amplitude,
        "accepted_steps": len(surface) - 1,
        "stop": dict(stop, time_over_T=float(stop["last_accepted_time_s"]) / T_S),
        "nonlinear_milestone": {
            "definition": "first accepted step of 50 consecutive accepted steps with H23=sqrt(h2^2+h3^2) >= 1e-5",
            "step": None if milestone is None else int(milestone["step"]),
            "time_s": None if milestone is None else float(milestone["time_s"]),
            "time_over_T": None if milestone is None else float(milestone["time_s"]) / T_S,
            "H23": None if milestone is None else float(milestone["H23"]),
        },
        "surface": {
            "mean_peak": float(peak["x_top_area_mean"]),
            "mean_peak_step": int(peak["step"]),
            "mean_peak_time_over_T": float(peak["time_s"]) / T_S,
            "mean_terminal": float(terminal["x_top_area_mean"]),
            "decline_after_peak": float(peak["x_top_area_mean"]) - float(terminal["x_top_area_mean"]),
            "decline_fraction": (float(peak["x_top_area_mean"]) - float(terminal["x_top_area_mean"])) / float(peak["x_top_area_mean"]),
            "max_terminal": float(terminal["x_top_max"]),
            "domain_max_terminal": float(terminal["x_domain_max"]),
        },
        "planes": planes,
        "independent_conservation_passed": bool(independent["passed"]),
        "independent_checks": independent,
        "sign_convention": "positive is downward; all plane fluxes are conservative SCVF integrations",
    }
    (out / "wp12_scalar_metrics.json").write_text(json.dumps(result, indent=2))
    with (out / "harmonic_signal.csv").open("w", newline="") as handle:
        fields = ["step", "time_s", "time_over_T", "weighted_fourier", "harmonic2_top", "harmonic3_top", "H23"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in harmonic_rows:
            writer.writerow({
                "step": row["step"], "time_s": row["time_s"], "time_over_T": float(row["time_s"]) / T_S,
                "weighted_fourier": row["weighted_fourier"], "harmonic2_top": row["harmonic2_top"],
                "harmonic3_top": row["harmonic3_top"], "H23": row["H23"],
            })
    with (out / "surface_series.csv").open("w", newline="") as handle:
        fields = ["grid", "amplitude", "step", "time_s", "time_over_T", "surface_mean", "surface_max", "domain_max", "surface_sigma"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in surface:
            writer.writerow({
                "grid": label, "amplitude": amplitude, "step": row["step"], "time_s": row["time_s"],
                "time_over_T": float(row["time_s"]) / T_S, "surface_mean": row["x_top_area_mean"],
                "surface_max": row["x_top_max"], "domain_max": row["x_domain_max"], "surface_sigma": row["sigma_top_node"],
            })
    with (out / "flux_series.csv").open("w", newline="") as handle:
        fields = ["grid", "amplitude", "depth_mm", "step", "time_s", "time_over_T", "advective_down_m_s", "total_down_m_s", "cumulative_advective_m", "cumulative_total_m"]
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for region, depth in DEPTHS.items():
            for row in budget:
                if row["region"] != region:
                    continue
                writer.writerow({
                    "grid": label, "amplitude": amplitude, "depth_mm": depth, "step": row["step"],
                    "time_s": row["time_s"], "time_over_T": float(row["time_s"]) / T_S,
                    "advective_down_m_s": row["cut_salt_adv_down_m_s"], "total_down_m_s": row["cut_salt_total_down_m_s"],
                    "cumulative_advective_m": row["cumulative_cut_adv_m"], "cumulative_total_m": row["cumulative_cut_salt_m"],
                })
    print(json.dumps({"grid": label, "amplitude": amplitude, "steps": result["accepted_steps"], "milestone_T": result["nonlinear_milestone"]["time_over_T"]}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("run", type=Path)
    parser.add_argument("out", type=Path)
    parser.add_argument("label")
    parser.add_argument("amplitude", type=float)
    parser.add_argument("--max-T", type=float)
    parser.add_argument("--post-dir", type=Path)
    args = parser.parse_args()
    main(args.run, args.out, args.label, args.amplitude, args.max_T, args.post_dir)
