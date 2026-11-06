"""Independent delivery acceptance for WP14 Phase B."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


PAPER2 = Path(__file__).resolve().parents[1]
BASE = PAPER2 / "analysis" / "wp14_nonlinear_timestep_20261106" / "phase_b"
INPUTS = BASE / "inputs"
PROVENANCE = BASE / "provenance"
METRICS = BASE / "metrics"
COMPARISON = BASE / "comparison"
FIGURES = BASE / "figures"
RAW = BASE / "raw"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def main() -> None:
    checks = []

    def check(name: str, condition: bool, detail=None) -> None:
        checks.append({"name": name, "passed": bool(condition), "detail": detail})

    plan = json.loads((INPUTS / "plan.json").read_text())
    diff = json.loads((INPUTS / "input_diff_audit.json").read_text())
    check("single authorized simulation", plan["production_simulation_limit"] == 1)
    check("automatic retry disabled", plan["automatic_retry_authorized"] is False)
    check("authorized R48 grid", plan["grid"] == "R48")
    check("authorized half timestep text", plan["nominal_dt_text"] == "56.842943238478475")
    check("authorized 14400 steps", plan["authorized_steps"] == 14400)
    check("authorized 6T horizon", plan["authorized_horizon_T"] == 6.0)
    check("input hash", sha256(INPUTS / "params.input") == plan["input_sha256"])
    check("only four approved input changes", set(diff["changed_keys"]) == {
        "Problem.Name", "TimeLoop.DtInitial", "TimeLoop.MaxTimeStepSize", "TimeLoop.OutputInterval"})
    check("all other input bytes preserved", diff["all_other_parameters_byte_preserved"] is True)

    reference_hashes = {
        "reference_R48_dt_products.tar.gz": "c5b409f6a319b70edccf2cfc82180fbbc25003847c6d751e250553a1977a41b7",
        "reference_R96_dt_products.tar.gz": "1db299d21fb780d99d26c6c1d0164dad9e349b7ec47929c81bd6c779ce3815bf",
    }
    for name, expected in reference_hashes.items():
        path = RAW / name
        check(f"reference archive exists {name}", path.is_file())
        if path.is_file():
            check(f"reference archive hash {name}", sha256(path) == expected)

    required_provenance = (
        "compute_preflight.json", "production_started.json", "finished.json", "archive_acceptance.json",
        "metrics_manifest.json", "sacct.txt", "submission.txt", "log_manifest.json",
        "preflight.stdout", "acceptance.stdout", "reduction.stdout", "remote_plan.json", "remote_params.input",
    )
    for name in required_provenance:
        check(f"provenance file {name}", (PROVENANCE / name).is_file())

    log_manifest = json.loads((PROVENANCE / "log_manifest.json").read_text())
    check("empty log hash registered", log_manifest["empty_sha256"] == hashlib.sha256(b"").hexdigest())
    check("all registered empty logs are zero bytes", all(value == 0 for value in log_manifest["logs"].values()))

    acceptance = json.loads((PROVENANCE / "archive_acceptance.json").read_text())
    finished = json.loads((PROVENANCE / "finished.json").read_text())
    manifest = json.loads((PROVENANCE / "metrics_manifest.json").read_text())
    check("independent archive acceptance", acceptance["passed"] is True, acceptance.get("failed"))
    for name, passed in acceptance["tests"].items():
        check(f"archive gate {name}", passed)
    check("solver completed without retry", finished["solver"]["returncode"] == 0 and not finished["solver"]["timed_out"])
    check("requested horizon stop", finished["solver"]["stop"]["reason"] == "requested_horizon")
    check("accepted 14400 steps", finished["solver"]["stop"]["accepted_steps"] == 14400)
    check("145 archived fields", acceptance["independent_checks"]["frames"] == 145)
    check("reduction acceptance", manifest["numerical_acceptance_passed"] is True)
    refined_archive = RAW / "refined_R48_dt2_products.tar.gz"
    check("refined products archive exists", refined_archive.is_file())
    if refined_archive.is_file():
        check("refined products archive hash", sha256(refined_archive) == finished["products"]["sha256"])

    expected_metric_rows = {
        "surface_series.csv": 14401,
        "harmonic_signal.csv": 14401,
        "field_state_series.csv": 145,
        "flux_series.csv": 57604,
    }
    for name, expected in expected_metric_rows.items():
        rows = csv_rows(METRICS / name)
        check(f"metric rows {name}", len(rows) == expected, len(rows))
    scalar = json.loads((METRICS / "wp12_scalar_metrics.json").read_text())
    field = json.loads((METRICS / "wp12_field_metrics.json").read_text())
    check("reduced scalar conservation", scalar["independent_conservation_passed"] is True)
    check("reduced field grid", field["grid"] == "R48")
    check("reduced field reaches 6T", abs(float(field["last_time_over_T"]) - 6.0) < 1e-9)

    expected_comparison_rows = {
        "matched_surface_and_H23.csv": 7201,
        "matched_field_modal_series.csv": 145,
        "matched_conservative_flux_series.csv": 28804,
        "matched_salt_inventory_series.csv": 28804,
        "matched_vertical_profiles.csv": None,
        "matched_representative_full_fields.csv": None,
        "finite_resolution_error_budget.csv": 49,
        "transport_pulse_dt_comparison.csv": 12,
    }
    for name, expected in expected_comparison_rows.items():
        rows = csv_rows(COMPARISON / name)
        check(f"comparison nonempty {name}", bool(rows), len(rows))
        if expected is not None:
            check(f"comparison rows {name}", len(rows) == expected, len(rows))
    summary = json.loads((COMPARISON / "phase_b_comparison_summary.json").read_text())
    for axis, error in summary["maximum_matching_time_error_s"].items():
        check(f"matching physical time {axis}", float(error) < 3e-7, error)
    check("all three conservation gates", all(item["passed"] for item in summary["conservation"].values()))

    figure_names = (
        "Figure_5_R48_dt_surface_and_modal_evolution",
        "Figure_6_R48_dt_conservative_transport",
        "Figure_7_spatial_vs_temporal_differences",
    )
    for name in figure_names:
        for suffix in (".pdf", ".svg", ".png", ".alignment.json", ".collision_audit.json", ".pdf_text_audit.json"):
            check(f"figure artifact {name}{suffix}", (FIGURES / f"{name}{suffix}").is_file())
        alignment = json.loads((FIGURES / f"{name}.alignment.json").read_text())
        collision = json.loads((FIGURES / f"{name}.collision_audit.json").read_text())
        text_audit = json.loads((FIGURES / f"{name}.pdf_text_audit.json").read_text())
        check(f"alignment pass {name}", alignment.get("verdict") == "PASS", alignment.get("verdict"))
        check(f"collision pass {name}", collision.get("verdict") == "PASS", collision.get("verdict"))
        check(f"font audit pass {name}", bool(text_audit.get("auditable"))
              and int(text_audit.get("below_minimum_count", -1)) == 0,
              {"minimum_found_pt": text_audit.get("minimum_found_pt"),
               "below_minimum_count": text_audit.get("below_minimum_count")})

    for name in ("REPORT_FINAL_ZH.md", "FIGURE_CONTRACT_PHASE_B.md", "QA_NOTES.md", "REPRODUCE.md",
                 "COST_RECORD.json", "PHASE_B_DECISION.json"):
        check(f"delivery file {name}", (BASE / name).is_file())

    failures = [item for item in checks if not item["passed"]]
    result = {
        "phase": "WP14 Phase B",
        "delivery_passed": not failures,
        "checks_count": len(checks),
        "failures": failures,
    }
    (BASE / "phase_b_acceptance.json").write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
