"""Reduce one WP11-compatible trajectory to auditable scalar/flux products.

The positive-flux summaries deliberately operate on accepted steps, not on the
span between the first and last positive sample.  Two definitions are retained:
strict sign (>0) and the previously registered 1e-15 m/s resolvability screen.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path

T_S = 545692.2550893936 / 4.0
DEPTH_REGIONS = {
    "upper_0.995000": 1.0,
    "upper_0.990000": 2.0,
    "upper_0.975000": 5.0,
    "upper_0.950000": 10.0,
    "upper_0.900000": 20.0,
}


def rows(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle))


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def positive_summary(records, field, threshold):
    accepted = [r for r in records if int(r["step"]) > 0]
    hits = [(r, float(r[field])) for r in accepted if float(r[field]) > threshold]
    if not hits:
        return {
            "threshold_m_s": threshold,
            "first_positive_step": None,
            "last_positive_step": None,
            "positive_step_count": 0,
            "longest_continuous_positive_interval": None,
            "peak_positive_flux_m_s": None,
            "peak_positive_step": None,
            "peak_positive_time_s": None,
            "peak_positive_time_over_T": None,
            "positive_step_support_duration_s": 0.0,
            "positive_step_support_duration_over_T": 0.0,
            "positive_step_integral_m": 0.0,
            "integration_semantics": "sum of F_i*dt_i over accepted steps with F_i above threshold",
        }

    runs, current = [], []
    previous_step = None
    hit_steps = {int(r["step"]): (r, value) for r, value in hits}
    for record in accepted:
        step = int(record["step"])
        if step in hit_steps:
            if current and previous_step is not None and step != previous_step + 1:
                runs.append(current)
                current = []
            current.append(hit_steps[step])
            previous_step = step
        elif current:
            runs.append(current)
            current = []
            previous_step = None
    if current:
        runs.append(current)
    longest = max(runs, key=lambda run: (sum(float(r["dt_s"]) for r, _ in run), len(run)))
    first, last = hits[0][0], hits[-1][0]
    longest_first, longest_last = longest[0][0], longest[-1][0]
    return {
        "threshold_m_s": threshold,
        "first_positive_step": int(first["step"]),
        "first_positive_time_s": float(first["time_s"]),
        "first_positive_time_over_T": float(first["time_s"]) / T_S,
        "last_positive_step": int(last["step"]),
        "last_positive_time_s": float(last["time_s"]),
        "last_positive_time_over_T": float(last["time_s"]) / T_S,
        "positive_step_count": len(hits),
        "positive_step_fraction": len(hits) / len(accepted),
        "longest_continuous_positive_interval": {
            "first_step": int(longest_first["step"]),
            "last_step": int(longest_last["step"]),
            "first_time_s": float(longest_first["time_s"]),
            "last_time_s": float(longest_last["time_s"]),
            "first_time_over_T": float(longest_first["time_s"]) / T_S,
            "last_time_over_T": float(longest_last["time_s"]) / T_S,
            "accepted_step_count": len(longest),
            "accepted_step_support_duration_s": sum(float(r["dt_s"]) for r, _ in longest),
            "accepted_step_support_duration_over_T": sum(float(r["dt_s"]) for r, _ in longest) / T_S,
        },
        "peak_positive_flux_m_s": max(value for _, value in hits),
        "peak_positive_step": int(max(hits, key=lambda item: item[1])[0]["step"]),
        "peak_positive_time_s": float(max(hits, key=lambda item: item[1])[0]["time_s"]),
        "peak_positive_time_over_T": float(max(hits, key=lambda item: item[1])[0]["time_s"]) / T_S,
        "positive_step_support_duration_s": sum(float(r["dt_s"]) for r, _ in hits),
        "positive_step_support_duration_over_T": sum(float(r["dt_s"]) for r, _ in hits) / T_S,
        "positive_step_integral_m": sum(value * float(r["dt_s"]) for r, value in hits),
        "integration_semantics": "sum of F_i*dt_i over accepted steps with F_i above threshold",
    }


def write_csv(path, fieldnames, records):
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)


def main(run, out, label, snapshot_manifest=None):
    run, out = Path(run), Path(out)
    out.mkdir(parents=True, exist_ok=True)
    surface = rows(run / "surface_events.csv")
    online = rows(run / "weighted_online.csv")
    budget = rows(run / "transport_budget.csv")
    if snapshot_manifest is None:
        stop = json.loads((run / "stop_record.json").read_text())
        boundary_file = run / "stop_record.json"
    else:
        boundary_file = Path(snapshot_manifest)
        snapshot = json.loads(boundary_file.read_text())
        stop = {
            "reason": "interim_snapshot_boundary",
            "last_accepted_step": int(snapshot["snapshot_end_step"]),
            "last_accepted_time_s": float(snapshot["snapshot_end_time_s"]),
            "global_threshold_low_s": None,
            "global_threshold_high_s": None,
            "terminal_state": False,
            "production_job_id": snapshot.get("production_job_id"),
            "scope": snapshot.get("scope", "read-only interim snapshot; not a terminal simulation state"),
        }
    online_by_step = {int(r["step"]): r for r in online}

    surface_series = []
    for record in surface:
        step = int(record["step"])
        mode = online_by_step[step]
        surface_series.append({
            "grid": label,
            "step": step,
            "time_s": float(record["time_s"]),
            "time_over_T": float(record["time_s"]) / T_S,
            "surface_mean": float(record["x_top_area_mean"]),
            "surface_max": float(record["x_top_max"]),
            "domain_max": float(record["x_domain_max"]),
            "surface_sigma": float(record["sigma_top_node"]),
            "weighted_fourier_m1": float(mode["weighted_fourier"]),
            "weighted_harmonic_m2": float(mode["harmonic2_top"]),
            "weighted_harmonic_m3": float(mode["harmonic3_top"]),
        })

    flux_series, plane_metrics = [], {}
    by_region = {}
    for record in budget:
        by_region.setdefault(record["region"], []).append(record)
    for region, depth_mm in DEPTH_REGIONS.items():
        records = by_region[region]
        for record in records:
            flux_series.append({
                "grid": label,
                "depth_mm": depth_mm,
                "step": int(record["step"]),
                "time_s": float(record["time_s"]),
                "time_over_T": float(record["time_s"]) / T_S,
                "advective_down_m_s": float(record["cut_salt_adv_down_m_s"]),
                "diffusive_down_m_s": float(record["cut_salt_diff_down_m_s"]),
                "total_down_m_s": float(record["cut_salt_total_down_m_s"]),
                "cumulative_advective_m": float(record["cumulative_cut_adv_m"]),
                "cumulative_diffusive_m": float(record["cumulative_cut_diff_m"]),
                "inventory_salt_m": float(record["inventory_salt_m"]),
                "salt_residual_m": float(record["salt_residual_m"]),
                "water_residual_m": float(record["water_residual_m"]),
            })
        plane_metrics[str(depth_mm)] = {
            "region": region,
            "requested_z_over_H": float(records[0]["requested_z_over_H"]),
            "actual_z_over_H": float(records[0]["actual_z_over_H"]),
            "actual_depth_mm": (1.0 - float(records[0]["actual_z_over_H"])) * 200.0,
            "advective": {
                "strict_positive": positive_summary(records, "cut_salt_adv_down_m_s", 0.0),
                "above_1e-15": positive_summary(records, "cut_salt_adv_down_m_s", 1e-15),
            },
            "total": {
                "strict_positive": positive_summary(records, "cut_salt_total_down_m_s", 0.0),
                "above_1e-15": positive_summary(records, "cut_salt_total_down_m_s", 1e-15),
            },
            "final": flux_series[-1],
            "maximum_abs_salt_residual_m": max(abs(float(r["salt_residual_m"])) for r in records),
            "maximum_abs_water_residual_m": max(abs(float(r["water_residual_m"])) for r in records),
        }

    result = {
        "schema_version": 1,
        "grid": label,
        "stop": dict(stop, time_over_T=stop["last_accepted_time_s"] / T_S),
        "accepted_steps": len(surface) - 1,
        "saved_scalar_rows": len(surface),
        "sign_convention": "positive is downward; values are signed plane-integrated conservative SCVF fluxes",
        "local_positive_negative_scvf_split_available": False,
        "positive_interval_note": "first/last are bounds only; longest continuous interval and actual positive-step count are computed independently",
        "planes": plane_metrics,
        "hashes": {name: sha256(run / name) for name in (
            "params.input", "eigenprofile.dat", "surface_events.csv", "weighted_online.csv",
            "transport_budget.csv")},
    }
    result["hashes"][boundary_file.name] = sha256(boundary_file)
    write_csv(out / "surface_series.csv", list(surface_series[0]), surface_series)
    write_csv(out / "flux_series.csv", list(flux_series[0]), flux_series)
    (out / "scalar_metrics.json").write_text(json.dumps(result, indent=2))
    print(json.dumps({"grid": label, "accepted_steps": result["accepted_steps"], "stop": result["stop"]}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("run")
    parser.add_argument("out")
    parser.add_argument("label")
    parser.add_argument("--snapshot-manifest")
    args = parser.parse_args()
    main(args.run, args.out, args.label, args.snapshot_manifest)
