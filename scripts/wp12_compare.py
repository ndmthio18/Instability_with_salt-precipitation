"""Assemble separated amplitude- and grid-sensitivity tables for WP12."""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path


def read_csv(path: Path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def load(tag: str, path: Path):
    scalar = json.loads((path / "wp12_scalar_metrics.json").read_text())
    field = json.loads((path / "wp12_field_metrics.json").read_text())
    return {
        "tag": tag,
        "path": path,
        "grid": scalar["grid"],
        "amplitude": scalar["amplitude"],
        "scalar": scalar,
        "field": field,
        "field_series": read_csv(path / "field_state_series.csv"),
        "harmonic": read_csv(path / "harmonic_signal.csv"),
    }


def delta(new, reference):
    if new is None or reference is None:
        return None, None
    difference = new - reference
    relative = difference / abs(reference) if reference != 0 else None
    return difference, relative


def nearest(rows, target):
    return min(rows, key=lambda row: abs(float(row["time_over_T"]) - target))


def metric_map(case):
    scalar, field = case["scalar"], case["field"]
    values = {
        "nonlinear_milestone_T": scalar["nonlinear_milestone"]["time_over_T"],
        "resolved_downflow_T": None if field["resolved_downflow_onset"] is None else field["resolved_downflow_onset"]["time_over_T"],
        "surface_mean_peak": scalar["surface"]["mean_peak"],
        "surface_peak_T": scalar["surface"]["mean_peak_time_over_T"],
        "surface_decline": scalar["surface"]["decline_after_peak"],
        "maximum_mirror_error": field["maximum_mirror_symmetry_error"],
        "maximum_odd_mode_fraction": field["maximum_odd_mode_energy_fraction"],
    }
    for depth, plane in scalar["planes"].items():
        for component in ("advective", "total"):
            summary = plane[component]
            prefix = f"{component}_{depth}mm"
            values[f"{prefix}_sustained_positive_T"] = summary["first_sustained_above_1e-15_time_over_T"]
            values[f"{prefix}_peak_m_s"] = summary["peak_positive_m_s"]
            values[f"{prefix}_integrated_positive_m"] = summary["integrated_positive_m"]
            values[f"{prefix}_integrated_signed_m"] = summary["integrated_signed_m"]
    return values


def main(out: Path, specs: list[tuple[str, Path]]) -> None:
    out.mkdir(parents=True, exist_ok=True)
    cases = [load(tag, path) for tag, path in specs]
    lookup = {(case["grid"], case["amplitude"]): case for case in cases}

    case_rows, plane_rows, mode_rows, aligned_rows = [], [], [], []
    for case in cases:
        scalar, field = case["scalar"], case["field"]
        onset = field["resolved_downflow_onset"]
        case_rows.append({
            "grid": case["grid"], "amplitude": case["amplitude"],
            "nonlinear_milestone_T": scalar["nonlinear_milestone"]["time_over_T"],
            "resolved_downflow_T": None if onset is None else onset["time_over_T"],
            "surface_mean_peak": scalar["surface"]["mean_peak"],
            "surface_peak_T": scalar["surface"]["mean_peak_time_over_T"],
            "surface_decline": scalar["surface"]["decline_after_peak"],
            "surface_decline_fraction": scalar["surface"]["decline_fraction"],
            "maximum_mirror_error": field["maximum_mirror_symmetry_error"],
            "maximum_odd_mode_fraction": field["maximum_odd_mode_energy_fraction"],
            "conservation_passed": scalar["independent_conservation_passed"],
        })
        for depth, plane in scalar["planes"].items():
            plane_rows.append({
                "grid": case["grid"], "amplitude": case["amplitude"], "depth_mm": depth,
                "advective_sustained_positive_T": plane["advective"]["first_sustained_above_1e-15_time_over_T"],
                "advective_peak_m_s": plane["advective"]["peak_positive_m_s"],
                "advective_positive_integral_m": plane["advective"]["integrated_positive_m"],
                "advective_signed_integral_m": plane["advective"]["integrated_signed_m"],
                "total_positive_integral_m": plane["total"]["integrated_positive_m"],
                "total_signed_integral_m": plane["total"]["integrated_signed_m"],
                "max_abs_salt_residual_m": plane["maximum_abs_salt_residual_m"],
                "max_abs_water_residual_m": plane["maximum_abs_water_residual_m"],
            })
        for target in (2.0, 3.0, 4.0, 5.0, 6.0):
            row = nearest(case["field_series"], target)
            mode_rows.append({
                "grid": case["grid"], "amplitude": case["amplitude"], "nominal_T": target,
                "actual_T": float(row["time_over_T"]), "dominant_mode": int(row["dominant_mode"]),
                "m2_fraction": float(row["mode_2_energy_fraction"]), "m4_fraction": float(row["mode_4_energy_fraction"]),
                "m6_fraction": float(row["mode_6_energy_fraction"]), "m10_fraction": float(row["mode_10_energy_fraction"]),
                "odd_fraction": float(row["odd_mode_energy_fraction"]), "mirror_error": float(row["mirror_symmetry_error"]),
            })
        milestone = scalar["nonlinear_milestone"]["time_over_T"]
        for row in case["harmonic"]:
            aligned_rows.append({
                "grid": case["grid"], "amplitude": case["amplitude"], "step": row["step"],
                "physical_T": row["time_over_T"],
                "time_minus_nonlinear_milestone_T": None if milestone is None else float(row["time_over_T"]) - milestone,
                "H23": row["H23"], "weighted_fourier": row["weighted_fourier"],
            })

    available_amplitudes = sorted({case["amplitude"] for case in cases})
    enhanced_amplitudes = [amplitude for amplitude in available_amplitudes if amplitude != 1e-6]
    amplitude_rows = []
    for grid in ("R24", "R48", "R96"):
        reference = lookup[(grid, 1e-6)]
        ref_values = metric_map(reference)
        for enhanced_amplitude in enhanced_amplitudes:
            new_values = metric_map(lookup[(grid, enhanced_amplitude)])
            for metric in ref_values:
                difference, relative = delta(new_values[metric], ref_values[metric])
                amplitude_rows.append({
                    "comparison": "amplitude_sensitivity_same_grid", "grid": grid, "metric": metric,
                    "reference_amplitude": 1e-6, "enhanced_amplitude": enhanced_amplitude,
                    "reference_value": ref_values[metric], "enhanced_value": new_values[metric],
                    "difference": difference, "relative_difference": relative,
                })

    grid_rows = []
    for amplitude in enhanced_amplitudes:
        for coarse, fine in (("R24", "R48"), ("R48", "R96"), ("R24", "R96")):
            coarse_values, fine_values = metric_map(lookup[(coarse, amplitude)]), metric_map(lookup[(fine, amplitude)])
            for metric in coarse_values:
                difference, relative = delta(coarse_values[metric], fine_values[metric])
                grid_rows.append({
                    "comparison": "grid_sensitivity_identical_amplitude", "amplitude": amplitude,
                    "coarse_grid": coarse, "fine_grid": fine, "metric": metric,
                    "coarse_value": coarse_values[metric], "fine_value": fine_values[metric],
                    "coarse_minus_fine": difference, "relative_to_fine": relative,
                })

    modal_sensitivity = []
    modes = (2, 4, 6, 10)
    for grid in ("R24", "R48", "R96"):
        for enhanced_amplitude in enhanced_amplitudes:
            for coordinate in ("physical_6T", "one_T_after_nonlinear_milestone"):
                states = {}
                for amplitude in (1e-6, enhanced_amplitude):
                    case = lookup[(grid, amplitude)]
                    target = 6.0 if coordinate == "physical_6T" else case["scalar"]["nonlinear_milestone"]["time_over_T"] + 1.0
                    states[amplitude] = nearest(case["field_series"], target)
                ref = [float(states[1e-6][f"mode_{mode}_energy_fraction"]) for mode in modes]
                new = [float(states[enhanced_amplitude][f"mode_{mode}_energy_fraction"]) for mode in modes]
                ref.append(max(0.0, 1.0 - sum(ref)))
                new.append(max(0.0, 1.0 - sum(new)))
                modal_sensitivity.append({
                    "comparison": "amplitude_sensitivity_same_grid", "grid": grid, "coordinate": coordinate,
                    "enhanced_amplitude": enhanced_amplitude,
                    "reference_actual_T": float(states[1e-6]["time_over_T"]), "enhanced_actual_T": float(states[enhanced_amplitude]["time_over_T"]),
                    "reference_dominant_mode": int(states[1e-6]["dominant_mode"]), "enhanced_dominant_mode": int(states[enhanced_amplitude]["dominant_mode"]),
                    "aggregated_modal_total_variation": 0.5 * sum(abs(left - right) for left, right in zip(ref, new)),
                    **{f"reference_m{mode}_fraction": ref[index] for index, mode in enumerate(modes)},
                    **{f"enhanced_m{mode}_fraction": new[index] for index, mode in enumerate(modes)},
                    "reference_other_fraction": ref[-1], "enhanced_other_fraction": new[-1],
                })

    write_csv(out / "case_summary.csv", case_rows)
    write_csv(out / "plane_transport_summary.csv", plane_rows)
    write_csv(out / "modal_checkpoints.csv", mode_rows)
    write_csv(out / "amplitude_sensitivity.csv", amplitude_rows)
    write_csv(out / "grid_sensitivity_by_amplitude.csv", grid_rows)
    write_csv(out / "modal_amplitude_sensitivity.csv", modal_sensitivity)
    write_csv(out / "milestone_aligned_harmonics.csv", aligned_rows)
    summary = {
        "schema_version": 1,
        "cases": [case["tag"] for case in cases],
        "amplitude_comparison": f"same grid only: 1e-6 versus {enhanced_amplitudes}",
        "grid_comparison": f"R24/R48/R96 separately at amplitudes {enhanced_amplitudes}",
        "alignment": "physical time and preregistered t-t_nl; no fitted or arbitrary shift",
        "case_summary": case_rows,
    }
    (out / "comparison_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps({"cases": len(cases), "amplitude_rows": len(amplitude_rows), "grid_rows": len(grid_rows)}, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("out", type=Path)
    parser.add_argument("cases", nargs="+")
    args = parser.parse_args()
    main(args.out, [(item.split("=", 1)[0], Path(item.split("=", 1)[1])) for item in args.cases])
