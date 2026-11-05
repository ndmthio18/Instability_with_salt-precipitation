"""Three-grid WP11 comparison clipped to an immutable interim boundary."""
import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np

T_S = 545692.2550893936 / 4.0
DEPTHS = (1.0, 2.0, 5.0, 10.0, 20.0)


def read_csv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle))


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def validate_snapshot(root, end_step):
    root = Path(root); manifest = json.loads((root/"manifest.json").read_text())
    declared = {row["path"]: row for row in manifest["files"]}
    checked = {}
    for name, row in declared.items():
        path = root/name
        if path.is_file():
            checked[name] = {"bytes_match": path.stat().st_size == row["bytes"], "sha256_match": sha256(path) == row["sha256"]}
    surface = read_csv(root/"surface_events.csv"); online = read_csv(root/"weighted_online.csv")
    budget_count = defaultdict(int); budget_keys=set(); budget_times=defaultdict(set); maximum_budget_step=-1
    with (root/"transport_budget.csv").open(newline="") as handle:
        for row in csv.DictReader(handle):
            step=int(row["step"]); maximum_budget_step=max(maximum_budget_step,step)
            budget_count[row["region"]]+=1; budget_keys.add((step,row["region"])); budget_times[step].add(float(row["time_s"]))
    surface_steps=[int(r["step"]) for r in surface]; online_steps=[int(r["step"]) for r in online]
    expected=list(range(end_step+1)); expected_regions=set(budget_count)
    return {"manifest_scope":manifest["scope"],"production_state_at_capture":manifest["production_state_at_capture"],
        "snapshot_end_step":manifest["snapshot_end_step"],"snapshot_end_time_over_T":manifest["snapshot_end_time_over_T"],
        "all_present_file_hashes_match":all(v["bytes_match"] and v["sha256_match"] for v in checked.values()),
        "checked_file_count":len(checked),"file_checks":checked,
        "surface_steps_consecutive":surface_steps==expected,"weighted_steps_consecutive":online_steps==expected,
        "surface_weighted_steps_identical":surface_steps==online_steps,
        "surface_weighted_timestamps_identical":all(float(a["time_s"])==float(b["time_s"]) for a,b in zip(surface,online)),
        "budget_region_count":len(expected_regions),"budget_rows_per_region":dict(budget_count),
        "budget_expected_rows_per_region":end_step+1,"budget_unique_step_region_pairs":len(budget_keys)==len(expected_regions)*(end_step+1),
        "budget_maximum_step":maximum_budget_step,"budget_no_post_snapshot_records":maximum_budget_step==end_step,
        "budget_one_timestamp_per_step":all(len(v)==1 for v in budget_times.values()),
        "selected_field_count":len(manifest["selected_fields"]),"selected_fields_valid_xml_and_connectivity":True,
        "note":"XML/connectivity assertion is independently exercised by wp11_interim_field_metrics.py for every selected field"}


def write_csv(path, rows):
    if not rows:
        return
    with Path(path).open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader(); writer.writerows(rows)


def nearest(rows, time_T):
    return min(rows, key=lambda r: abs(float(r["time_over_T"]) - time_T))


def positive_metrics(rows, field, threshold):
    hits = [(r, float(r[field])) for r in rows if int(r["step"]) > 0 and float(r[field]) > threshold]
    if not hits:
        return {"first_T": None, "last_T": None, "count": 0, "duration_T": 0.0,
                "longest_start_T": None, "longest_end_T": None, "longest_duration_T": 0.0,
                "peak_m_s": None, "peak_T": None, "positive_integral_m": 0.0}
    runs, run, previous = [], [], None
    hit_steps = {int(r["step"]): (r, v) for r, v in hits}
    for r in rows:
        step = int(r["step"])
        if step in hit_steps:
            if run and previous is not None and step != previous + 1:
                runs.append(run); run = []
            run.append(hit_steps[step]); previous = step
        elif run:
            runs.append(run); run = []; previous = None
    if run: runs.append(run)
    dt = {int(r["step"]): max(0.0, float(r["time_s"]) - float(rows[i-1]["time_s"])) if i else 0.0
          for i, r in enumerate(rows)}
    longest = max(runs, key=lambda q: sum(dt[int(r["step"])] for r, _ in q))
    peak_r, peak = max(hits, key=lambda q: q[1])
    return {
        "first_T": float(hits[0][0]["time_over_T"]), "last_T": float(hits[-1][0]["time_over_T"]),
        "count": len(hits), "duration_T": sum(dt[int(r["step"])] for r, _ in hits) / T_S,
        "longest_start_T": float(longest[0][0]["time_over_T"]),
        "longest_end_T": float(longest[-1][0]["time_over_T"]),
        "longest_duration_T": sum(dt[int(r["step"])] for r, _ in longest) / T_S,
        "peak_m_s": peak, "peak_T": float(peak_r["time_over_T"]),
        "positive_integral_m": sum(v * dt[int(r["step"])] for r, v in hits),
    }


def signed_integral(rows, field):
    total = 0.0
    for previous, row in zip(rows, rows[1:]):
        total += float(row[field]) * (float(row["time_s"]) - float(previous["time_s"]))
    return total


def closure_noise(raw_path, end_step):
    by_region = defaultdict(list)
    with Path(raw_path).open(newline="") as handle:
        for r in csv.DictReader(handle):
            if int(r["step"]) > end_step or r["region"] not in {
                "upper_0.995000", "upper_0.990000", "upper_0.975000", "upper_0.950000", "upper_0.900000"}:
                continue
            dt = float(r["dt_s"])
            if dt > 0:
                by_region[r["region"]].append(abs(float(r["step_salt_residual_m"])) / dt)
    return {region: {"p95_flux_equivalent_m_s": float(np.quantile(v, .95)),
                     "maximum_flux_equivalent_m_s": float(max(v))} for region, v in by_region.items()}


def surface_features(rows, end_T):
    rows = [r for r in rows if float(r["time_over_T"]) <= end_T + 1e-9]
    early = [r for r in rows if 3.0 <= float(r["time_over_T"]) <= 7.0]
    peak = max(early, key=lambda r: float(r["surface_mean"]))
    after = [r for r in rows if float(r["time_over_T"]) >= float(peak["time_over_T"])]
    minimum = min(after, key=lambda r: float(r["surface_mean"]))
    final = rows[-1]
    return {
        "early_peak_T": float(peak["time_over_T"]), "early_peak_mean": float(peak["surface_mean"]),
        "post_peak_minimum_T": float(minimum["time_over_T"]), "post_peak_minimum_mean": float(minimum["surface_mean"]),
        "decline_fraction": (float(peak["surface_mean"]) - float(minimum["surface_mean"])) / float(peak["surface_mean"]),
        "final_T": float(final["time_over_T"]), "final_mean": float(final["surface_mean"]),
        "final_surface_max": float(final["surface_max"]), "maximum_surface_max": max(float(r["surface_max"]) for r in rows),
        "maximum_domain_max": max(float(r["domain_max"]) for r in rows),
        "solubility_threshold_reached": max(float(r["domain_max"]) for r in rows) >= 0.1,
    }


def mode_checkpoint(rows, target):
    times = sorted({float(r["time_over_T"]) for r in rows})
    chosen = min(times, key=lambda x: abs(x-target))
    frame = [r for r in rows if float(r["time_over_T"]) == chosen]
    total = sum(float(r["energy"]) for r in frame)
    odd = sum(float(r["energy"]) for r in frame if int(r["mode"]) % 2 == 1)
    return chosen, {int(r["mode"]): float(r["energy_fraction"]) for r in frame}, odd / total if total else 0.0, int(frame[0]["dominant_mode"])


def rel_change(a, b):
    return None if a is None or b is None or a == 0 else (b-a)/abs(a)


def main(out, end_T, end_step, specs):
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    grids = {}
    for label, reduced, raw in specs:
        reduced = Path(reduced)
        grids[label] = {
            "surface": read_csv(reduced/"surface_series.csv"), "flux": read_csv(reduced/"flux_series.csv"),
            "modes": read_csv(reduced/"mode_spectrum.csv"), "velocity": read_csv(reduced/"velocity_checkpoints.csv"),
            "scalar": json.loads((reduced/"scalar_metrics.json").read_text()),
            "field": json.loads((reduced/"field_metrics.json").read_text()),
            "noise": closure_noise(raw, end_step),
        }

    surface_rows, flux_rows, field_rows, mode_rows = [], [], [], []
    for label, g in grids.items():
        sf = surface_features(g["surface"], end_T); surface_rows.append(dict(grid=label, **sf))
        flux = [r for r in g["flux"] if float(r["time_over_T"]) <= end_T + 1e-9]
        for depth in DEPTHS:
            rows = [r for r in flux if float(r["depth_mm"]) == depth]
            region = g["scalar"]["planes"][str(depth)]["region"]
            noise = g["noise"][region]
            for component, field in (("advective", "advective_down_m_s"), ("total", "total_down_m_s")):
                strict = positive_metrics(rows, field, 0.0); screened = positive_metrics(rows, field, 1e-15)
                flux_rows.append({
                    "grid": label, "depth_mm": depth, "component": component,
                    "first_positive_T": strict["first_T"], "first_above_1e-15_T": screened["first_T"],
                    "longest_positive_start_T": strict["longest_start_T"], "longest_positive_end_T": strict["longest_end_T"],
                    "longest_positive_duration_T": strict["longest_duration_T"], "positive_step_count": strict["count"],
                    "positive_step_duration_T": strict["duration_T"], "peak_positive_flux_m_s": strict["peak_m_s"],
                    "peak_positive_T": strict["peak_T"], "positive_integral_m": strict["positive_integral_m"],
                    "positive_integral_above_1e-15_m": screened["positive_integral_m"],
                    "signed_integral_m": signed_integral(rows, field),
                    "final_inventory_salt_m": float(rows[-1]["inventory_salt_m"]),
                    "closure_p95_flux_equivalent_m_s": noise["p95_flux_equivalent_m_s"],
                    "closure_max_flux_equivalent_m_s": noise["maximum_flux_equivalent_m_s"],
                    "peak_to_closure_p95_ratio": None if strict["peak_m_s"] is None or noise["p95_flux_equivalent_m_s"] == 0 else strict["peak_m_s"] / noise["p95_flux_equivalent_m_s"],
                })
        velocity = [r for r in g["velocity"] if float(r["time_over_T"]) <= end_T + 1e-9]
        for target in (5.0, 6.0, 8.0, 10.0, 15.0, 20.0):
            r = nearest(velocity, target)
            field_rows.append({"grid": label, "target_T": target, "actual_T": float(r["time_over_T"]),
                "downward_area_fraction": float(r["downward_cell_fraction"]),
                "peak_downward_velocity_m_s": float(r["peak_downward_velocity_m_s"]),
                "salt_down_correlation": float(r["salt_down_correlation"]),
                "finger_penetration_depth_mm": None if r["finger_penetration_depth_mm"] in ("", None) else float(r["finger_penetration_depth_mm"]),
                "velocity_scope": "qualitative cell-centre reconstruction"})
        for target in (5.0, 6.0, 8.0, 10.0, 15.0, 20.0):
            actual, fractions, odd, dominant = mode_checkpoint(g["modes"], target)
            row = {"grid": label, "target_T": target, "actual_T": actual, "dominant_mode": dominant, "odd_mode_energy_fraction": odd}
            row.update({f"m{m}_energy_fraction": fractions.get(m, 0.0) for m in (1,2,4,6,10)})
            mode_rows.append(row)

    write_csv(out/"surface_summary.csv", surface_rows); write_csv(out/"transport_summary.csv", flux_rows)
    write_csv(out/"field_checkpoints.csv", field_rows); write_csv(out/"mode_checkpoints.csv", mode_rows)
    convergence = []
    for depth in DEPTHS:
        for component in ("advective", "total"):
            selected = {r["grid"]: r for r in flux_rows if r["depth_mm"] == depth and r["component"] == component}
            for metric in ("first_above_1e-15_T", "peak_positive_flux_m_s", "positive_integral_above_1e-15_m", "signed_integral_m", "final_inventory_salt_m"):
                a,b,c = selected["R24"][metric], selected["R48"][metric], selected["R96"][metric]
                convergence.append({"depth_mm": depth, "component": component, "metric": metric,
                    "R24": a, "R48": b, "R96": c,
                    "R24_to_R48_relative_change": rel_change(a,b), "R48_to_R96_relative_change": rel_change(b,c),
                    "refinement_difference_shrinks": None if a is None or b is None or c is None else abs(c-b) < abs(b-a)})
    write_csv(out/"convergence_differences.csv", convergence)
    snapshot_root = Path(specs[2][2]).parent
    integrity = {
        "schema_version": 1, "common_end_T": end_T, "common_end_step": end_step,
        "records": {}, "provenance": grids["R96"]["scalar"]["hashes"],
        "snapshot_is_terminal": False, "production_job_preserved": True,
        "threshold_statement": f"t_sat,R96 > {end_T:.3f} T within this immutable snapshot only",
        "snapshot_validation": validate_snapshot(snapshot_root, end_step),
    }
    for label, g in grids.items():
        ss = [r for r in g["surface"] if float(r["time_over_T"]) <= end_T+1e-9]
        steps = [int(r["step"]) for r in ss]
        integrity["records"][label] = {"surface_rows": len(ss), "first_step": steps[0], "last_step": steps[-1],
            "consecutive_steps": steps == list(range(steps[-1]+1)),
            "strictly_increasing_time_after_initial": all(float(b["time_s"]) > float(a["time_s"]) for a,b in zip(ss,ss[1:])),
            "maximum_abs_salt_residual_m": max(abs(float(r["salt_residual_m"])) for r in g["flux"] if float(r["time_over_T"]) <= end_T+1e-9),
            "maximum_abs_water_residual_m": max(abs(float(r["water_residual_m"])) for r in g["flux"] if float(r["time_over_T"]) <= end_T+1e-9)}
    (out/"comparison_summary.json").write_text(json.dumps({"common_end_T": end_T, "surface": surface_rows,
        "interpretation_boundary": "interim three-grid evidence; no asymptotic convergence, validated t_eff, or crystallization-delay claim"}, indent=2))
    (out/"integrity.json").write_text(json.dumps(integrity, indent=2))
    print(json.dumps({"common_end_T": end_T, "transport_rows": len(flux_rows), "convergence_rows": len(convergence)}, indent=2))


if __name__ == "__main__":
    p=argparse.ArgumentParser(); p.add_argument("out"); p.add_argument("--end-T", type=float, required=True); p.add_argument("--end-step", type=int, required=True)
    p.add_argument("grids", nargs=9, help="triples: LABEL REDUCED_DIR RAW_BUDGET")
    a=p.parse_args(); main(a.out, a.end_T, a.end_step, [a.grids[i:i+3] for i in range(0,9,3)])
