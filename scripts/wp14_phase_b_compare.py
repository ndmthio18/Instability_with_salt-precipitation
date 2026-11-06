"""Compare the single WP14 R48 dt/2 run with frozen R48/R96 dt references."""
from __future__ import annotations

import csv
import io
import json
import math
import tarfile
from pathlib import Path

import numpy as np


PAPER2 = Path(__file__).resolve().parents[1]
BASE = PAPER2 / "analysis" / "wp14_nonlinear_timestep_20261106"
PHASE_A = BASE / "phase_a"
PHASE_B = BASE / "phase_b"
OUT = PHASE_B / "comparison"
WP12 = PAPER2 / "analysis" / "wp12_amplitude_convergence_20261106"
REF_R48 = WP12 / "enhanced_metrics" / "R48"
REF_R96 = WP12 / "enhanced_metrics" / "R96"
NEW_R48 = PHASE_B / "metrics"
RAW = PHASE_B / "raw"
ARCHIVES = {
    "R48_dt": RAW / "reference_R48_dt_products.tar.gz",
    "R96_dt": RAW / "reference_R96_dt_products.tar.gz",
    "R48_dt2": RAW / "refined_R48_dt2_products.tar.gz",
}
T_S = 545692.2550893936 / 4.0
DT = 113.68588647695695
FLUX_FLOOR = 1e-15
SUSTAIN_DURATION_S = 50 * DT
RATIO_THRESHOLD = 0.1
RATIO_DURATION_T = 0.25
DEPTH_REGION = {1.0: "upper_0.995000", 2.0: "upper_0.990000", 5.0: "upper_0.975000", 10.0: "upper_0.950000"}


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise AssertionError(f"empty output: {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def tar_member(handle: tarfile.TarFile, suffix: str):
    matches = [m for m in handle.getmembers() if m.name == suffix or m.name.endswith("/" + suffix)]
    if len(matches) != 1:
        raise RuntimeError(f"expected one {suffix}, got {[m.name for m in matches]}")
    stream = handle.extractfile(matches[0])
    if stream is None:
        raise RuntimeError(suffix)
    return stream


def load_budget(archive: Path) -> list[dict[str, str]]:
    with tarfile.open(archive, "r:gz") as handle:
        stream = tar_member(handle, "transport_budget.csv")
        text = io.TextIOWrapper(stream, encoding="utf-8", newline="")
        return [row for row in csv.DictReader(text) if row["region"] in DEPTH_REGION.values()]


def load_npz(archive: Path, suffix: str) -> dict[str, np.ndarray]:
    with tarfile.open(archive, "r:gz") as handle:
        data = tar_member(handle, suffix).read()
    with np.load(io.BytesIO(data)) as bundle:
        return {key: bundle[key].copy() for key in bundle.files}


def load_case(root: Path, label: str, budget: list[dict[str, str]]) -> dict:
    return {
        "label": label,
        "root": root,
        "scalar": json.loads((root / "wp12_scalar_metrics.json").read_text()),
        "field": json.loads((root / "wp12_field_metrics.json").read_text()),
        "field_series": read_csv(root / "field_state_series.csv"),
        "surface": read_csv(root / "surface_series.csv"),
        "harmonic": read_csv(root / "harmonic_signal.csv"),
        "flux": read_csv(root / "flux_series.csv"),
        "budget": budget,
    }


def normalized_fields(case: dict) -> list[dict]:
    initial = float(case["field_series"][0]["mode_2_energy"])
    floor = max(1e-30, 1e-8 * initial)
    rows = []
    for raw in case["field_series"]:
        seed = float(raw["mode_2_energy"])
        total = float(raw["total_nonzero_mode_energy"])
        nonseed = max(0.0, total - seed)
        rows.append({
            "step": int(raw["step"]), "time_s": float(raw["time_s"]),
            "time_over_T": float(raw["time_over_T"]),
            "R_harm": math.sqrt(nonseed / seed) if seed >= floor and seed > 0 else None,
            "dominant_mode": int(raw["dominant_mode"]),
            "total_nonzero_energy": total,
            "m2": float(raw["mode_2_energy_fraction"]),
            "m4": float(raw["mode_4_energy_fraction"]),
            "m6": float(raw["mode_6_energy_fraction"]),
            "m10": float(raw["mode_10_energy_fraction"]),
            "odd": float(raw["odd_mode_energy_fraction"]),
            "mirror": float(raw["mirror_symmetry_error"]),
            "peak_downward_velocity_m_s": float(raw["peak_downward_velocity_m_s"]),
        })
    return rows


def ratio_crossing(rows: list[dict]) -> dict | None:
    for index, row in enumerate(rows):
        if row["R_harm"] is None or row["R_harm"] < RATIO_THRESHOLD:
            continue
        target = row["time_over_T"] + RATIO_DURATION_T
        for later in rows[index + 1:]:
            if later["R_harm"] is None or later["R_harm"] < RATIO_THRESHOLD:
                break
            if later["time_over_T"] + 1e-12 >= target:
                return row
    return None


def nearest(rows: list[dict], target: float) -> dict:
    return min(rows, key=lambda row: abs(float(row["time_over_T"]) - target))


def first_sustained_time(times: list[float], values: list[float]) -> float | None:
    for start in range(1, len(values)):
        if values[start] <= FLUX_FLOOR:
            continue
        target = times[start] + SUSTAIN_DURATION_S
        for index in range(start + 1, len(values)):
            if values[index] <= FLUX_FLOOR:
                break
            if times[index] + 1e-8 >= target:
                return times[start] / T_S
    return None


def pulse_metrics(case: dict, depth: float) -> dict:
    raw = sorted((r for r in case["flux"] if float(r["depth_mm"]) == depth), key=lambda r: int(r["step"]))
    times = [float(r["time_s"]) for r in raw]
    adv = [float(r["advective_down_m_s"]) for r in raw]
    total = [float(r["total_down_m_s"]) for r in raw]
    dts = [0.0] + [right - left for left, right in zip(times, times[1:])]
    pulses, active = [], []
    for index in range(1, len(adv)):
        if adv[index] > FLUX_FLOOR:
            active.append(index)
        elif active:
            pulses.append(active)
            active = []
    if active:
        pulses.append(active)
    summaries = []
    for pulse in pulses:
        peak = max(pulse, key=lambda i: adv[i])
        summaries.append({
            "area": sum(adv[i] * dts[i] for i in pulse),
            "start_T": (times[pulse[0]] - dts[pulse[0]]) / T_S,
            "end_T": times[pulse[-1]] / T_S,
            "peak_T": times[peak] / T_S,
            "peak": adv[peak],
        })
    principal = max(summaries, key=lambda item: item["area"]) if summaries else None
    positive = [max(value, 0.0) * dt for value, dt in zip(adv, dts)]
    cumulative, running = [], 0.0
    for increment in positive:
        running += increment
        cumulative.append(running)
    fractions = {}
    resolved_integral = running >= FLUX_FLOOR * SUSTAIN_DURATION_S
    for fraction in (0.05, 0.10, 0.50, 0.90):
        value = None
        if resolved_integral:
            target = fraction * running
            value = next(times[i] / T_S for i, amount in enumerate(cumulative) if amount + 1e-30 >= target)
        fractions[f"t{int(100*fraction):02d}_T"] = value
    peak_index = max(range(1, len(adv)), key=lambda i: adv[i])
    return {
        "first_resolved_T": next((times[i] / T_S for i in range(1, len(adv)) if adv[i] > FLUX_FLOOR), None),
        "first_sustained_physical_duration_T": first_sustained_time(times, adv),
        "principal_start_T": None if principal is None else principal["start_T"],
        "principal_end_T": None if principal is None else principal["end_T"],
        "principal_peak_T": None if principal is None else principal["peak_T"],
        "principal_integral_m": None if principal is None else principal["area"],
        "peak_positive_m_s": max(0.0, adv[peak_index]),
        "peak_positive_T": times[peak_index] / T_S if adv[peak_index] > 0 else None,
        "positive_integral_m": running,
        "signed_advective_integral_m": sum(v * dt for v, dt in zip(adv, dts)),
        "signed_total_integral_m": sum(v * dt for v, dt in zip(total, dts)),
        **fractions,
    }


def final_inventory(case: dict, depth: float) -> float:
    region = DEPTH_REGION[depth]
    rows = [r for r in case["budget"] if r["region"] == region]
    return float(max(rows, key=lambda r: int(r["step"]))["inventory_salt_m"])


def case_qois(case: dict) -> tuple[dict[str, float | None], dict[float, dict]]:
    normalized = normalized_fields(case)
    crossing = ratio_crossing(normalized)
    state6 = nearest(normalized, 6.0)
    surface = case["scalar"]["surface"]
    down = case["field"]["resolved_downflow_onset"]
    q = {
        "normalized_R0.1_T": None if crossing is None else crossing["time_over_T"],
        "resolved_downflow_T": None if down is None else float(down["time_over_T"]),
        "surface_mean_peak": float(surface["mean_peak"]),
        "surface_peak_T": float(surface["mean_peak_time_over_T"]),
        "surface_decline": float(surface["decline_after_peak"]),
        "surface_terminal_mean": float(surface["mean_terminal"]),
        "surface_terminal_max": float(surface["max_terminal"]),
        "maximum_mirror_error": float(case["field"]["maximum_mirror_symmetry_error"]),
        "m2_fraction_6T": state6["m2"], "m4_fraction_6T": state6["m4"],
        "m6_fraction_6T": state6["m6"], "m10_fraction_6T": state6["m10"],
        "odd_fraction_6T": state6["odd"],
    }
    pulses = {}
    for depth in (1.0, 2.0, 5.0, 10.0):
        p = pulse_metrics(case, depth)
        pulses[depth] = p
        prefix = f"{int(depth)}mm"
        for key in ("principal_peak_T", "principal_integral_m", "peak_positive_m_s", "positive_integral_m",
                    "signed_advective_integral_m", "signed_total_integral_m", "t50_T", "t90_T"):
            q[f"{prefix}_{key}"] = p[key]
        q[f"{prefix}_final_salt_inventory_m"] = final_inventory(case, depth)
    return q, pulses


def paired_rows(reference: list[dict], refined: list[dict], refined_step_factor: int) -> tuple[list[tuple[dict, dict]], float]:
    lookup = {int(row["step"]): row for row in refined}
    pairs, max_time_error = [], 0.0
    for row in reference:
        other = lookup[int(row["step"]) * refined_step_factor]
        error = abs(float(row["time_s"]) - float(other["time_s"]))
        max_time_error = max(max_time_error, error)
        pairs.append((row, other))
    return pairs, max_time_error


def budget_by_depth(rows: list[dict], depth: float) -> list[dict]:
    return sorted((r for r in rows if r["region"] == DEPTH_REGION[depth]), key=lambda r: int(r["step"]))


def widths(values: np.ndarray) -> np.ndarray:
    return np.diff(np.r_[values[0], (values[1:] + values[:-1]) / 2, values[-1]])


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for path in ARCHIVES.values():
        if not path.is_file():
            raise FileNotFoundError(path)
    budgets = {label: load_budget(path) for label, path in ARCHIVES.items()}
    cases = {
        "R48_dt": load_case(REF_R48, "R48_dt", budgets["R48_dt"]),
        "R96_dt": load_case(REF_R96, "R96_dt", budgets["R96_dt"]),
        "R48_dt2": load_case(NEW_R48, "R48_dt2", budgets["R48_dt2"]),
    }
    qois, pulses = {}, {}
    for label, case in cases.items():
        qois[label], pulses[label] = case_qois(case)

    surface_pairs, surface_time_error = paired_rows(cases["R48_dt"]["surface"], cases["R48_dt2"]["surface"], 2)
    harmonic_pairs, harmonic_time_error = paired_rows(cases["R48_dt"]["harmonic"], cases["R48_dt2"]["harmonic"], 2)
    surface_rows = []
    for (ref_s, new_s), (_, new_h) in zip(surface_pairs, harmonic_pairs):
        surface_rows.append({
            "reference_step": int(ref_s["step"]), "refined_step": int(new_s["step"]),
            "time_s": float(ref_s["time_s"]), "time_over_T": float(ref_s["time_over_T"]),
            "surface_mean_dt": float(ref_s["surface_mean"]), "surface_mean_dt2": float(new_s["surface_mean"]),
            "surface_mean_difference": float(new_s["surface_mean"]) - float(ref_s["surface_mean"]),
            "surface_max_dt": float(ref_s["surface_max"]), "surface_max_dt2": float(new_s["surface_max"]),
            "H23_dt": float(cases["R48_dt"]["harmonic"][int(ref_s["step"])]["H23"]),
            "H23_dt2": float(new_h["H23"]),
        })
    write_csv(OUT / "matched_surface_and_H23.csv", surface_rows)

    field_ref, field_new = normalized_fields(cases["R48_dt"]), normalized_fields(cases["R48_dt2"])
    field_pairs, field_time_error = paired_rows(field_ref, field_new, 2)
    field_rows = []
    modal_tvs = []
    dominant_disagreements = 0
    for ref, new in field_pairs:
        ref_modes = [ref[key] for key in ("m2", "m4", "m6", "m10")]
        new_modes = [new[key] for key in ("m2", "m4", "m6", "m10")]
        ref_modes.append(max(0.0, 1 - sum(ref_modes)))
        new_modes.append(max(0.0, 1 - sum(new_modes)))
        tv = 0.5 * sum(abs(a - b) for a, b in zip(ref_modes, new_modes))
        modal_tvs.append(tv)
        dominant_disagreements += int(ref["dominant_mode"] != new["dominant_mode"])
        field_rows.append({
            "reference_step": ref["step"], "refined_step": new["step"],
            "time_s": ref["time_s"], "time_over_T": ref["time_over_T"],
            "R_harm_dt": ref["R_harm"], "R_harm_dt2": new["R_harm"],
            "dominant_mode_dt": ref["dominant_mode"], "dominant_mode_dt2": new["dominant_mode"],
            "modal_total_variation": tv,
            "m2_dt": ref["m2"], "m2_dt2": new["m2"],
            "m4_dt": ref["m4"], "m4_dt2": new["m4"],
            "m6_dt": ref["m6"], "m6_dt2": new["m6"],
            "m10_dt": ref["m10"], "m10_dt2": new["m10"],
            "odd_dt": ref["odd"], "odd_dt2": new["odd"],
            "mirror_dt": ref["mirror"], "mirror_dt2": new["mirror"],
            "peak_downward_velocity_dt": ref["peak_downward_velocity_m_s"],
            "peak_downward_velocity_dt2": new["peak_downward_velocity_m_s"],
        })
    write_csv(OUT / "matched_field_modal_series.csv", field_rows)

    flux_rows, inventory_rows = [], []
    flux_time_errors, inventory_time_errors = [], []
    for depth in (1.0, 2.0, 5.0, 10.0):
        ref_flux = [r for r in cases["R48_dt"]["flux"] if float(r["depth_mm"]) == depth]
        new_flux = [r for r in cases["R48_dt2"]["flux"] if float(r["depth_mm"]) == depth]
        pairs, error = paired_rows(ref_flux, new_flux, 2)
        flux_time_errors.append(error)
        for ref, new in pairs:
            flux_rows.append({
                "depth_mm": depth, "reference_step": int(ref["step"]), "refined_step": int(new["step"]),
                "time_s": float(ref["time_s"]), "time_over_T": float(ref["time_over_T"]),
                "advective_dt": float(ref["advective_down_m_s"]), "advective_dt2": float(new["advective_down_m_s"]),
                "advective_difference": float(new["advective_down_m_s"]) - float(ref["advective_down_m_s"]),
                "total_dt": float(ref["total_down_m_s"]), "total_dt2": float(new["total_down_m_s"]),
                "cumulative_advective_dt": float(ref["cumulative_advective_m"]),
                "cumulative_advective_dt2": float(new["cumulative_advective_m"]),
                "cumulative_total_dt": float(ref["cumulative_total_m"]),
                "cumulative_total_dt2": float(new["cumulative_total_m"]),
            })
        ref_budget = budget_by_depth(cases["R48_dt"]["budget"], depth)
        new_budget = budget_by_depth(cases["R48_dt2"]["budget"], depth)
        pairs, error = paired_rows(ref_budget, new_budget, 2)
        inventory_time_errors.append(error)
        for ref, new in pairs:
            inventory_rows.append({
                "depth_mm": depth, "reference_step": int(ref["step"]), "refined_step": int(new["step"]),
                "time_s": float(ref["time_s"]), "time_over_T": float(ref["time_s"]) / T_S,
                "salt_inventory_dt": float(ref["inventory_salt_m"]),
                "salt_inventory_dt2": float(new["inventory_salt_m"]),
                "salt_inventory_difference": float(new["inventory_salt_m"]) - float(ref["inventory_salt_m"]),
            })
    write_csv(OUT / "matched_conservative_flux_series.csv", flux_rows)
    write_csv(OUT / "matched_salt_inventory_series.csv", inventory_rows)

    profiles = {label: load_npz(archive, "post/profiles.npz") for label, archive in ARCHIVES.items()}
    profile_rows = []
    for target in (3.0, 4.0, 5.0, 6.0):
        reference = profiles["R48_dt"]
        refined = profiles["R48_dt2"]
        i_ref = int(np.argmin(abs(reference["time_s"] / T_S - target)))
        i_new = int(np.argmin(abs(refined["time_s"] / T_S - target)))
        if abs(reference["time_s"][i_ref] - refined["time_s"][i_new]) > 1e-6:
            raise AssertionError("profile time mismatch")
        for index, z in enumerate(reference["z_m"]):
            profile_rows.append({
                "target_T": target, "actual_T": reference["time_s"][i_ref] / T_S, "z_m": float(z),
                "x_mean_dt": float(reference["x_mean"][i_ref, index]),
                "x_mean_dt2": float(refined["x_mean"][i_new, index]),
                "x_sigma_dt": float(reference["x_sigma"][i_ref, index]),
                "x_sigma_dt2": float(refined["x_sigma"][i_new, index]),
            })
    write_csv(OUT / "matched_vertical_profiles.csv", profile_rows)

    representative = {label: load_npz(archive, "post/representative_fields.npz") for label, archive in ARCHIVES.items()}
    full_field_rows = []
    ref_rep, new_rep = representative["R48_dt"], representative["R48_dt2"]
    common = sorted(set(key for key in ref_rep if key.startswith("x_")) & set(key for key in new_rep if key.startswith("x_")))
    wx, wz = widths(ref_rep["xs"]), widths(ref_rep["z_m"])
    weights = np.outer(wz, wx)
    for key in common:
        difference = new_rep[key] - ref_rep[key]
        reference_norm = math.sqrt(float(np.sum(ref_rep[key] ** 2 * weights)))
        full_field_rows.append({
            "field_key": key,
            "max_abs_difference": float(np.max(abs(difference))),
            "weighted_L2_difference": math.sqrt(float(np.sum(difference ** 2 * weights))),
            "relative_weighted_L2": math.sqrt(float(np.sum(difference ** 2 * weights))) / reference_norm,
        })
    write_csv(OUT / "matched_representative_full_fields.csv", full_field_rows)

    qoi_rows = []
    all_metrics = sorted(qois["R48_dt"])
    for metric in all_metrics:
        r48, r96, dt2 = qois["R48_dt"][metric], qois["R96_dt"][metric], qois["R48_dt2"][metric]
        if r48 is None or r96 is None or dt2 is None:
            dx = dt_error = ratio = None
        else:
            dx, dt_error = abs(r96 - r48), abs(dt2 - r48)
            ratio = dt_error / dx if dx > 0 else None
        if ratio is None:
            interpretation = "undefined"
        elif ratio < 0.25 - 1e-8:
            interpretation = "temporal_difference_much_smaller"
        elif ratio <= 4 + 1e-12:
            interpretation = "comparable_finite_resolution_differences"
        else:
            interpretation = "temporal_difference_larger"
        qoi_rows.append({
            "metric": metric, "R48_dt": r48, "R96_dt": r96, "R48_dt2": dt2,
            "delta_x": dx, "delta_t": dt_error, "delta_t_over_delta_x": ratio,
            "interpretation": interpretation,
        })
    write_csv(OUT / "finite_resolution_error_budget.csv", qoi_rows)

    pulse_rows = []
    for label in ("R48_dt", "R48_dt2", "R96_dt"):
        for depth, values in pulses[label].items():
            pulse_rows.append({"case": label, "depth_mm": depth, **values})
    write_csv(OUT / "transport_pulse_dt_comparison.csv", pulse_rows)

    summary = {
        "schema_version": 1,
        "matching": "R48 dt step k versus R48 dt/2 step 2k; field outputs 50 versus 100 steps",
        "maximum_matching_time_error_s": {
            "surface": surface_time_error, "harmonic": harmonic_time_error,
            "field": field_time_error, "flux": max(flux_time_errors), "inventory": max(inventory_time_errors),
        },
        "normalized_milestone_T": {label: qois[label]["normalized_R0.1_T"] for label in qois},
        "resolved_downflow_T": {label: qois[label]["resolved_downflow_T"] for label in qois},
        "dominant_mode_6T": {label: int(nearest(normalized_fields(cases[label]), 6.0)["dominant_mode"])
                              for label in cases},
        "R48_dt_vs_dt2_modal": {
            "dominant_mode_disagreement_frames": dominant_disagreements,
            "frames": len(field_pairs),
            "maximum_modal_total_variation": max(modal_tvs),
            "median_modal_total_variation": float(np.median(modal_tvs)),
        },
        "conservation": {
            label: {
                "passed": case["scalar"]["independent_conservation_passed"],
                "max_salt_reconstruction_m": max(
                    plane["maximum_abs_salt_residual_m"] for plane in case["scalar"]["planes"].values()),
                "max_water_reconstruction_m": max(
                    plane["maximum_abs_water_residual_m"] for plane in case["scalar"]["planes"].values()),
            } for label, case in cases.items()
        },
        "interpretation_thresholds": {
            "delta_t_much_smaller_than_delta_x": "delta_t/delta_x < 0.25",
            "comparable": "0.25 <= delta_t/delta_x <= 4",
            "temporal_larger": "delta_t/delta_x > 4",
        },
        "limitations": [
            "finite-resolution differences are not truncation-error bounds",
            "spatial and temporal errors are not assumed additive or independent",
            "no convergence order is inferred if modal branches differ",
            "pulse timing is not a validated physical t_eff",
        ],
    }
    (OUT / "phase_b_comparison_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps({
        "matched_surface_rows": len(surface_rows), "matched_field_rows": len(field_rows),
        "matched_flux_rows": len(flux_rows), "matched_inventory_rows": len(inventory_rows),
        "qois": len(qoi_rows), "summary": summary["R48_dt_vs_dt2_modal"],
    }, indent=2))


if __name__ == "__main__":
    main()
