"""WP14 Phase A: read-only normalized nonlinear and transport-pulse audit.

This script consumes only the frozen WP12/WP13 reduced products.  It does not
open or modify any DuMuX production run and never fits a time shift.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path


PAPER2 = Path(__file__).resolve().parents[1]
WP12 = PAPER2 / "analysis" / "wp12_amplitude_convergence_20261106"
WP13 = PAPER2 / "analysis" / "wp13_amplitude_1e5_20261106"
OUT = PAPER2 / "analysis" / "wp14_nonlinear_timestep_20261106" / "phase_a"
T_S = 545692.2550893936 / 4.0
GRIDS = ("R24", "R48", "R96")
AMPLITUDES = (1e-6, 1e-5, 1e-4)
RATIO_THRESHOLDS = (0.01, 0.1, 1.0)
PRIMARY_RATIO_THRESHOLD = 0.1
REGISTERED_DURATION_T = 0.25
ADVECTIVE_FLUX_FLOOR = 1e-15
SUSTAINED_STEPS = 50
NOMINAL_DT_S = 113.68588647695695
CUMULATIVE_RESOLUTION_FLOOR_M = ADVECTIVE_FLUX_FLOOR * SUSTAINED_STEPS * NOMINAL_DT_S


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise AssertionError(f"refusing to write empty table: {path}")
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def case_dir(grid: str, amplitude: float) -> Path:
    if amplitude == 1e-6:
        return WP12 / "reference_metrics" / grid
    if amplitude == 1e-5:
        return WP13 / "metrics" / grid
    if amplitude == 1e-4:
        return WP12 / "enhanced_metrics" / grid
    raise ValueError(amplitude)


def load_case(grid: str, amplitude: float) -> dict:
    root = case_dir(grid, amplitude)
    scalar = json.loads((root / "wp12_scalar_metrics.json").read_text())
    field = json.loads((root / "wp12_field_metrics.json").read_text())
    fields = read_csv(root / "field_state_series.csv")
    harmonics = read_csv(root / "harmonic_signal.csv")
    flux = read_csv(root / "flux_series.csv")
    if len(fields) < 2 or len(harmonics) != 7201:
        raise AssertionError(f"unexpected frozen trajectory length for {grid}, {amplitude:g}")
    return {
        "grid": grid,
        "amplitude": amplitude,
        "root": root,
        "scalar": scalar,
        "field": field,
        "fields": fields,
        "harmonics": harmonics,
        "flux": flux,
    }


def nearest(rows: list[dict], target_t: float) -> dict:
    return min(rows, key=lambda row: abs(float(row["time_over_T"]) - target_t))


def normalized_rows(case: dict) -> list[dict]:
    initial_seed = float(case["fields"][0]["mode_2_energy"])
    seed_floor = max(1e-30, initial_seed * 1e-8)
    rows = []
    for raw in case["fields"]:
        seed = float(raw["mode_2_energy"])
        total = float(raw["total_nonzero_mode_energy"])
        nonseed = max(0.0, total - seed)
        valid = seed >= seed_floor
        ratio = math.sqrt(nonseed / seed) if valid and seed > 0 else None
        odd = float(raw["odd_mode_energy_fraction"])
        rows.append({
            "grid": case["grid"],
            "amplitude": case["amplitude"],
            "frame": int(raw["frame"]),
            "step": int(raw["step"]),
            "time_s": float(raw["time_s"]),
            "time_over_T": float(raw["time_over_T"]),
            "seed_energy_E2": seed,
            "nonseed_energy": nonseed,
            "seed_energy_floor": seed_floor,
            "ratio_valid": int(valid),
            "R_harm": ratio,
            "total_nonzero_energy": total,
            "even_mode_fraction": max(0.0, 1.0 - odd),
            "odd_mode_fraction": odd,
            "mirror_symmetry_error": float(raw["mirror_symmetry_error"]),
            "dominant_mode": int(raw["dominant_mode"]),
            "m2_fraction": float(raw["mode_2_energy_fraction"]),
            "m4_fraction": float(raw["mode_4_energy_fraction"]),
            "m6_fraction": float(raw["mode_6_energy_fraction"]),
            "m10_fraction": float(raw["mode_10_energy_fraction"]),
        })
    return rows


def sustained_ratio_crossing(rows: list[dict], threshold: float) -> dict | None:
    """First saved state remaining above threshold for 0.25T without gaps."""
    for start, row in enumerate(rows):
        if not row["ratio_valid"] or row["R_harm"] is None or row["R_harm"] < threshold:
            continue
        target = row["time_over_T"] + REGISTERED_DURATION_T
        previous_step = row["step"]
        reached_duration = False
        for later in rows[start + 1:]:
            if later["step"] <= previous_step:
                raise AssertionError("non-monotone field steps")
            previous_step = later["step"]
            if not later["ratio_valid"] or later["R_harm"] is None or later["R_harm"] < threshold:
                break
            if later["time_over_T"] + 1e-12 >= target:
                reached_duration = True
                break
        if reached_duration:
            return row
    return None


def state_row(case: dict, normalized: list[dict], coordinate: str, target: float | None) -> dict:
    if target is None:
        return {
            "grid": case["grid"], "amplitude": case["amplitude"],
            "coordinate": coordinate, "target": None, "actual_T": None,
            "dominant_mode": None, "R_harm": None, "total_nonzero_energy": None,
            "m2_fraction": None, "m4_fraction": None, "m6_fraction": None,
            "m10_fraction": None, "odd_mode_fraction": None,
            "mirror_symmetry_error": None,
        }
    row = nearest(normalized, target)
    return {
        "grid": case["grid"], "amplitude": case["amplitude"],
        "coordinate": coordinate, "target": target,
        "actual_T": row["time_over_T"], "dominant_mode": row["dominant_mode"],
        "R_harm": row["R_harm"], "total_nonzero_energy": row["total_nonzero_energy"],
        "m2_fraction": row["m2_fraction"], "m4_fraction": row["m4_fraction"],
        "m6_fraction": row["m6_fraction"], "m10_fraction": row["m10_fraction"],
        "odd_mode_fraction": row["odd_mode_fraction"],
        "mirror_symmetry_error": row["mirror_symmetry_error"],
    }


def first_sustained_index(values: list[float], steps: list[int]) -> int | None:
    run_start = None
    run_length = 0
    previous_step = None
    for index, (value, step) in enumerate(zip(values, steps)):
        if value > ADVECTIVE_FLUX_FLOOR and (previous_step is None or step == previous_step + 1):
            if run_length == 0:
                run_start = index
            run_length += 1
            if run_length >= SUSTAINED_STEPS:
                return run_start
        elif value > ADVECTIVE_FLUX_FLOOR:
            run_start, run_length = index, 1
        else:
            run_start, run_length = None, 0
        previous_step = step
    return None


def threshold_time(times: list[float], cumulative: list[float], target: float) -> float | None:
    for time, value in zip(times, cumulative):
        if value + 1e-30 >= target:
            return time / T_S
    return None


def transport_metrics(case: dict, depth: float) -> dict:
    raw = [row for row in case["flux"] if float(row["depth_mm"]) == depth]
    raw.sort(key=lambda row: int(row["step"]))
    steps = [int(row["step"]) for row in raw]
    times = [float(row["time_s"]) for row in raw]
    adv = [float(row["advective_down_m_s"]) for row in raw]
    total = [float(row["total_down_m_s"]) for row in raw]
    if steps != list(range(steps[-1] + 1)) or steps[-1] != 7200:
        raise AssertionError(f"non-contiguous frozen flux series: {case['grid']} {case['amplitude']} {depth}")
    dts = [0.0] + [right - left for left, right in zip(times, times[1:])]
    positive_increments = [max(value, 0.0) * dt for value, dt in zip(adv, dts)]
    cumulative = []
    running = 0.0
    for increment in positive_increments:
        running += increment
        cumulative.append(running)
    final_positive = cumulative[-1]
    first_resolved = next((i for i in range(1, len(adv)) if adv[i] > ADVECTIVE_FLUX_FLOOR), None)
    sustained = first_sustained_index(adv[1:], steps[1:])
    sustained = None if sustained is None else sustained + 1

    pulses = []
    active = []
    for index in range(1, len(adv)):
        if adv[index] > ADVECTIVE_FLUX_FLOOR:
            active.append(index)
        elif active:
            pulses.append(active)
            active = []
    if active:
        pulses.append(active)
    pulse_summaries = []
    for pulse in pulses:
        peak = max(pulse, key=lambda i: adv[i])
        pulse_summaries.append({
            "indices": pulse,
            "area": sum(adv[i] * dts[i] for i in pulse),
            "duration_s": sum(dts[i] for i in pulse),
            "start_s": times[pulse[0]] - dts[pulse[0]],
            "end_s": times[pulse[-1]],
            "peak": adv[peak],
            "peak_s": times[peak],
        })
    principal = max(pulse_summaries, key=lambda item: item["area"]) if pulse_summaries else None
    peak_index = max(range(1, len(adv)), key=lambda i: adv[i])
    resolved_final = final_positive >= CUMULATIVE_RESOLUTION_FLOOR_M
    fractions = {}
    for fraction in (0.05, 0.10, 0.50, 0.90):
        fractions[f"t_{int(fraction * 100):02d}_T"] = (
            threshold_time(times, cumulative, fraction * final_positive) if resolved_final else None
        )
    return {
        "grid": case["grid"], "amplitude": case["amplitude"], "depth_mm": depth,
        "flux_resolution_floor_m_s": ADVECTIVE_FLUX_FLOOR,
        "cumulative_resolution_floor_m": CUMULATIVE_RESOLUTION_FLOOR_M,
        "first_resolved_positive_T": None if first_resolved is None else times[first_resolved] / T_S,
        "first_sustained_positive_T": None if sustained is None else times[sustained] / T_S,
        "positive_pulse_count": len(pulse_summaries),
        "principal_pulse_start_T": None if principal is None else principal["start_s"] / T_S,
        "principal_pulse_end_T": None if principal is None else principal["end_s"] / T_S,
        "principal_pulse_duration_T": None if principal is None else principal["duration_s"] / T_S,
        "principal_pulse_integral_m": None if principal is None else principal["area"],
        "principal_pulse_peak_m_s": None if principal is None else principal["peak"],
        "principal_pulse_peak_T": None if principal is None else principal["peak_s"] / T_S,
        "peak_positive_m_s": max(0.0, adv[peak_index]),
        "peak_positive_T": times[peak_index] / T_S if adv[peak_index] > 0 else None,
        "positive_advective_integral_m": final_positive,
        "resolved_positive_advective_integral_m": sum(
            value * dt for value, dt in zip(adv, dts) if value > ADVECTIVE_FLUX_FLOOR
        ),
        "signed_advective_integral_m": sum(value * dt for value, dt in zip(adv, dts)),
        "signed_total_integral_m": sum(value * dt for value, dt in zip(total, dts)),
        "cumulative_fraction_times_defined": int(resolved_final),
        **fractions,
    }


def relative_difference(left: float | None, right: float | None, floor: float) -> tuple[float | None, float | None]:
    if left is None or right is None:
        return None, None
    absolute = abs(left - right)
    return absolute, absolute / max(abs(right), floor)


def metric_difference(metric: str, left: float | None, right: float | None, floor: float) -> tuple[float | None, float | None]:
    """Finite difference, with unresolved positive-transport pairs left undefined."""
    positive_signal_metric = (
        "positive_advective" in metric
        or "principal_pulse_integral" in metric
        or "peak_advective" in metric
    )
    if positive_signal_metric and left is not None and right is not None:
        if max(abs(left), abs(right)) < floor:
            return None, None
    return relative_difference(left, right, floor)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    cases = {(g, a): load_case(g, a) for a in AMPLITUDES for g in GRIDS}
    source_files = (
        "wp12_scalar_metrics.json", "wp12_field_metrics.json", "field_state_series.csv",
        "harmonic_signal.csv", "flux_series.csv", "surface_series.csv",
    )
    source_manifest = []
    for amplitude in AMPLITUDES:
        for grid in GRIDS:
            root = cases[(grid, amplitude)]["root"]
            for name in source_files:
                path = root / name
                source_manifest.append({
                    "grid": grid, "amplitude": amplitude,
                    "path_relative_to_paper2": path.relative_to(PAPER2).as_posix(),
                    "bytes": path.stat().st_size, "sha256": sha256(path),
                })
    (OUT / "source_manifest.json").write_text(json.dumps(source_manifest, indent=2))
    normalized = {(g, a): normalized_rows(case) for (g, a), case in cases.items()}
    write_csv(OUT / "normalized_harmonic_timeseries.csv", [
        row for a in AMPLITUDES for g in GRIDS for row in normalized[(g, a)]
    ])

    milestones = []
    crossing_lookup = {}
    mode_rows = []
    for amplitude in AMPLITUDES:
        for grid in GRIDS:
            case = cases[(grid, amplitude)]
            norm = normalized[(grid, amplitude)]
            absolute_t = case["scalar"]["nonlinear_milestone"]["time_over_T"]
            for threshold in RATIO_THRESHOLDS:
                crossing = sustained_ratio_crossing(norm, threshold)
                crossing_lookup[(grid, amplitude, threshold)] = crossing
                milestones.append({
                    "grid": grid, "amplitude": amplitude,
                    "diagnostic": "R_harm", "threshold": threshold,
                    "registered_duration_T": REGISTERED_DURATION_T,
                    "crossing_T": None if crossing is None else crossing["time_over_T"],
                    "crossing_step": None if crossing is None else crossing["step"],
                    "R_harm_at_crossing": None if crossing is None else crossing["R_harm"],
                    "absolute_H23_milestone_T": absolute_t,
                })
            for target in (2.0, 3.0, 4.0, 5.0, 6.0):
                mode_rows.append(state_row(case, norm, f"physical_{target:g}T", target))
            primary = crossing_lookup[(grid, amplitude, PRIMARY_RATIO_THRESHOLD)]
            mode_rows.append(state_row(
                case, norm, "normalized_Rharm_0.1",
                None if primary is None else primary["time_over_T"],
            ))
            mode_rows.append(state_row(
                case, norm, "absolute_H23_plus_1T_auxiliary",
                None if absolute_t is None else absolute_t + 1.0,
            ))
    write_csv(OUT / "normalized_milestones.csv", milestones)
    write_csv(OUT / "mode_checkpoints.csv", mode_rows)

    transport = [
        transport_metrics(cases[(grid, amplitude)], depth)
        for amplitude in AMPLITUDES for grid in GRIDS for depth in (1.0, 2.0, 5.0, 10.0)
    ]
    write_csv(OUT / "transport_pulse_metrics.csv", transport)
    transport_lookup = {(r["grid"], r["amplitude"], r["depth_mm"]): r for r in transport}
    state_lookup = {(r["grid"], r["amplitude"], r["coordinate"]): r for r in mode_rows}

    values = {}
    for amplitude in AMPLITUDES:
        for grid in GRIDS:
            case = cases[(grid, amplitude)]
            crossing = crossing_lookup[(grid, amplitude, PRIMARY_RATIO_THRESHOLD)]
            scalar = case["scalar"]
            state6 = state_lookup[(grid, amplitude, "physical_6T")]
            item = {
                "normalized_milestone_R0.1_T": None if crossing is None else crossing["time_over_T"],
                "surface_mean_peak": scalar["surface"]["mean_peak"],
                "surface_peak_T": scalar["surface"]["mean_peak_time_over_T"],
                "surface_decline": scalar["surface"]["decline_after_peak"],
                "maximum_mirror_error": case["field"]["maximum_mirror_symmetry_error"],
                "m2_fraction_6T": state6["m2_fraction"],
                "m4_fraction_6T": state6["m4_fraction"],
                "m6_fraction_6T": state6["m6_fraction"],
                "m10_fraction_6T": state6["m10_fraction"],
                "odd_fraction_6T": state6["odd_mode_fraction"],
            }
            for depth in (5.0, 10.0):
                pulse = transport_lookup[(grid, amplitude, depth)]
                prefix = f"{int(depth)}mm"
                item[f"{prefix}_principal_pulse_peak_T"] = pulse["principal_pulse_peak_T"]
                item[f"{prefix}_principal_pulse_integral_m"] = pulse["principal_pulse_integral_m"]
                item[f"{prefix}_peak_advective_m_s"] = pulse["peak_positive_m_s"]
                item[f"{prefix}_positive_advective_integral_m"] = pulse["positive_advective_integral_m"]
                item[f"{prefix}_signed_total_integral_m"] = pulse["signed_total_integral_m"]
            values[(grid, amplitude)] = item

    metric_floors = {
        "normalized_milestone_R0.1_T": 1 / 24,
        "surface_mean_peak": 1e-12,
        "surface_peak_T": 1 / 1200,
        "surface_decline": 1e-12,
        "maximum_mirror_error": 1e-12,
        "m2_fraction_6T": 1e-8,
        "m4_fraction_6T": 1e-8,
        "m6_fraction_6T": 1e-8,
        "m10_fraction_6T": 1e-8,
        "odd_fraction_6T": 1e-12,
    }
    for depth in (5, 10):
        metric_floors[f"{depth}mm_principal_pulse_peak_T"] = 1 / 1200
        metric_floors[f"{depth}mm_principal_pulse_integral_m"] = CUMULATIVE_RESOLUTION_FLOOR_M
        metric_floors[f"{depth}mm_peak_advective_m_s"] = ADVECTIVE_FLUX_FLOOR
        metric_floors[f"{depth}mm_positive_advective_integral_m"] = CUMULATIVE_RESOLUTION_FLOOR_M
        metric_floors[f"{depth}mm_signed_total_integral_m"] = CUMULATIVE_RESOLUTION_FLOOR_M

    grid_rows = []
    for amplitude in AMPLITUDES:
        for metric, floor in metric_floors.items():
            v24, v48, v96 = (values[(grid, amplitude)][metric] for grid in GRIDS)
            e24_48, r24_48 = metric_difference(metric, v24, v48, floor)
            e48_96, r48_96 = metric_difference(metric, v48, v96, floor)
            grid_rows.append({
                "amplitude": amplitude, "metric": metric, "comparison_floor": floor,
                "R24_value": v24, "R48_value": v48, "R96_value": v96,
                "abs_difference_R24_R48": e24_48,
                "abs_difference_R48_R96": e48_96,
                "relative_difference_R24_R48": r24_48,
                "relative_difference_R48_R96": r48_96,
                "adjacent_difference_decreases": (
                    None if e24_48 is None or e48_96 is None else int(e48_96 < e24_48)
                ),
            })
    write_csv(OUT / "grid_difference_by_amplitude.csv", grid_rows)

    amplitude_rows = []
    for grid in GRIDS:
        reference = crossing_lookup[(grid, 1e-6, PRIMARY_RATIO_THRESHOLD)]
        for amplitude in (1e-5, 1e-4):
            enhanced = crossing_lookup[(grid, amplitude, PRIMARY_RATIO_THRESHOLD)]
            ref_t = None if reference is None else reference["time_over_T"]
            new_t = None if enhanced is None else enhanced["time_over_T"]
            amplitude_rows.append({
                "grid": grid, "reference_amplitude": 1e-6,
                "enhanced_amplitude": amplitude,
                "reference_normalized_milestone_T": ref_t,
                "enhanced_normalized_milestone_T": new_t,
                "relative_time_change": None if ref_t is None or new_t is None else (new_t - ref_t) / ref_t,
                "absolute_H23_reference_T": cases[(grid, 1e-6)]["scalar"]["nonlinear_milestone"]["time_over_T"],
                "absolute_H23_enhanced_T": cases[(grid, amplitude)]["scalar"]["nonlinear_milestone"]["time_over_T"],
            })
    write_csv(OUT / "amplitude_normalized_milestone_comparison.csv", amplitude_rows)

    summary = {
        "schema_version": 1,
        "analysis_type": "read-only reanalysis of nine frozen 0-6T trajectories",
        "cases": len(cases),
        "grids": list(GRIDS),
        "amplitudes": list(AMPLITUDES),
        "fourier_basis": "Neumann cosine modes cos(m*pi*x/W), m>=1, applied to x-horizontal_mean(z)",
        "vertical_and_horizontal_weights": "control-volume widths from point coordinates; modal energy is width-weighted over x and z",
        "normalized_ratio": "R_harm=sqrt((sum_m>=1 E_m-E_2)/E_2)",
        "seed_mode": "m=2, the physical A28 mode in the reduced box",
        "seed_floor": "per trajectory max(1e-30, 1e-8 times initial E2); ratio undefined below floor",
        "ratio_thresholds": list(RATIO_THRESHOLDS),
        "primary_ratio_threshold": PRIMARY_RATIO_THRESHOLD,
        "registered_duration_T": REGISTERED_DURATION_T,
        "positive_flux_floor_m_s": ADVECTIVE_FLUX_FLOOR,
        "sustained_positive_steps": SUSTAINED_STEPS,
        "cumulative_positive_fraction_denominator": "frozen C_adv+(z,6T), computed over the same 0-6T window",
        "cumulative_resolution_floor_m": CUMULATIVE_RESOLUTION_FLOOR_M,
        "principal_pulse": "contiguous J_adv>1e-15 interval with largest resolved positive integral",
        "warnings": [
            "normalized thresholds are predefined and not fitted to improve agreement",
            "saved-field cadence is 50 accepted steps; normalized milestone timing is quantized at 1/24 T",
            "odd-mode suppression and grid agreement are associated observations, not a causal identification",
            "adjacent finite-resolution differences are not Richardson convergence or truncation-error bounds",
            "transport-pulse timings are diagnostics, not a validated physical t_eff",
        ],
        "source_commit": "acee5fc6d384b7c6719331a077dd87bbc77afebb",
    }
    (OUT / "phase_a_protocol_and_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps({
        "cases": len(cases),
        "normalized_rows": sum(len(v) for v in normalized.values()),
        "milestones": len(milestones),
        "mode_rows": len(mode_rows),
        "transport_rows": len(transport),
        "grid_rows": len(grid_rows),
    }, indent=2))


if __name__ == "__main__":
    main()
