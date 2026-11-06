"""Verify the nine-case amplitude benchmark delivery without re-running DuMuX."""
from __future__ import annotations

import csv
import json
import py_compile
from pathlib import Path


ROOT = Path(__file__).resolve().parent.parent
WP12 = ROOT / "analysis" / "wp12_amplitude_convergence_20261106"
WP13 = ROOT / "analysis" / "wp13_amplitude_1e5_20261106"
COMPARISON = WP12 / "comparison"
FIGURES = WP12 / "figures"
GRIDS = ("R24", "R48", "R96")
AMPLITUDES = (1e-6, 1e-5, 1e-4)


def rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def case_path(grid: str, amplitude: float) -> Path:
    if amplitude == 1e-6:
        return WP12 / "reference_metrics" / grid
    if amplitude == 1e-5:
        return WP13 / "metrics" / grid
    return WP12 / "enhanced_metrics" / grid


def main() -> None:
    checks: dict[str, bool] = {}
    for grid in GRIDS:
        for amplitude in AMPLITUDES:
            label = f"{grid}/A{amplitude:.0e}"
            path = case_path(grid, amplitude)
            scalar = json.loads((path / "wp12_scalar_metrics.json").read_text())
            field = json.loads((path / "wp12_field_metrics.json").read_text())
            checks[f"{label}/identity"] = scalar["grid"] == grid and scalar["amplitude"] == amplitude
            checks[f"{label}/conservation"] = scalar["independent_conservation_passed"] is True
            checks[f"{label}/surface_rows"] = len(rows(path / "surface_series.csv")) == 7201
            checks[f"{label}/harmonic_rows"] = len(rows(path / "harmonic_signal.csv")) == 7201
            checks[f"{label}/field_rows"] = len(rows(path / "field_state_series.csv")) == 145
            checks[f"{label}/field_metrics"] = (
                field["grid"] == grid
                and field["frames"] == 145
                and "resolved_downflow_onset" in field
                and "maximum_mirror_symmetry_error" in field
                and "maximum_odd_mode_energy_fraction" in field
            )
            if amplitude != 1e-6:
                acceptance = json.loads((path / "archive_acceptance.json").read_text())
                checks[f"{label}/archive_acceptance"] = acceptance["passed"] and not acceptance["failed"]
                checks[f"{label}/requested_horizon"] = (
                    acceptance["stop"]["reason"] == "requested_horizon"
                    and acceptance["stop"]["accepted_steps"] == 7200
                )
                checks[f"{label}/terminal_clamp"] = acceptance["timestep_audit"]["terminal_horizon_clamp"]

    expected_rows = {
        "case_summary.csv": 9,
        "plane_transport_summary.csv": 36,
        "modal_checkpoints.csv": 45,
        "amplitude_sensitivity.csv": 234,
        "grid_sensitivity_by_amplitude.csv": 234,
        "modal_amplitude_sensitivity.csv": 12,
        "milestone_aligned_harmonics.csv": 64809,
    }
    for name, expected in expected_rows.items():
        checks[f"comparison/{name}"] = len(rows(COMPARISON / name)) == expected

    for stem in (
        "Figure_1_timing_and_registered_alignment",
        "Figure_2_A1e-4_grid_transport",
        "Figure_3_modal_energy_physical_and_aligned",
        "Figure_4_symmetry_and_timing",
        "Figure_5_A1e-5_grid_transport",
    ):
        for suffix in ("pdf", "png", "svg", "alignment.json", "pdf_text_audit.json", "collision_audit.json"):
            checks[f"figure/{stem}.{suffix}"] = (FIGURES / f"{stem}.{suffix}").exists()

    for name in ("REPORT_ZH.md", "NUMERICAL_AND_COST_CHECKS.md", "REPRODUCE.md"):
        checks[f"document/{name}"] = (WP12 / name).exists()
    checks["document/FIGURE_DESCRIPTIONS.md"] = (FIGURES / "FIGURE_DESCRIPTIONS.md").exists()

    for name in (
        "wp12_prepare_amplitude.py", "wp13_prepare_amplitude_1e5.py",
        "wp12_extract_reference.py", "wp12_reduce_archive.py",
        "wp12_archive_acceptance.py", "wp12_scalar_metrics.py",
        "wp12_field_metrics.py", "wp12_compare.py", "wp12_make_figures.py",
        "verify_wp12_delivery.py",
    ):
        py_compile.compile(str(ROOT / "scripts" / name), doraise=True)
    checks["scripts/python_compile"] = True

    result = {
        "schema_version": 1,
        "delivery_passed": all(checks.values()),
        "checks_count": len(checks),
        "failures": [name for name, passed in checks.items() if not passed],
        "checks": checks,
        "claim_boundary": "numerical and delivery acceptance; not physical validation of solubility-crossing time",
    }
    (WP12 / "delivery_verification.json").write_text(json.dumps(result, indent=2))
    print(json.dumps({key: result[key] for key in ("delivery_passed", "checks_count", "failures")}, indent=2))
    if not result["delivery_passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
