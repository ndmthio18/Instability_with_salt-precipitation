"""Acceptance checks for the P2-H k_r/a_v comparison delivery."""
import csv
import json
import math
import pathlib

root = pathlib.Path.cwd()
analysis = root / "analysis"
results = analysis / "results"
expected = {
    "kr0p1_av1": (1e-4, 1.0),
    "kr1_av1_ref": (1e-3, 1.0),
    "kr10_av1": (1e-2, 1.0),
    "kr1_av0p5": (1e-3, 0.5),
    "kr1_av2": (1e-3, 2.0),
}
checks = []

def check(name, passed, detail):
    checks.append({"check": name, "passed": bool(passed), "detail": str(detail)})
    if not passed:
        raise AssertionError(f"{name}: {detail}")

with (results / "summary.csv").open(newline="") as handle:
    summary = {row["case"]: row for row in csv.DictReader(handle)}
check("five unique runs", set(summary) == set(expected), sorted(summary))
matrix = json.loads((results / "case_matrix.json").read_text())
check("reference reused", matrix["reference_reused"] and len(matrix["unique_runs"]) == 5, matrix)
provenance = json.loads((analysis / "provenance.json").read_text())
check("corrected matrix provenance", provenance["corrected_matrix_job"] == "16408105", provenance)
check("solver checksum", len(provenance["solver_sha256"]) == 64, provenance["solver_sha256"])
with (results / "comparison_spread.csv").open(newline="") as handle:
    spread_rows = list(csv.DictReader(handle))
max_relative_spread = max(float(row["relative_range"]) for row in spread_rows)
check("pre-ALE parameter independence", max_relative_spread <= 1e-12, max_relative_spread)

for name, (kr, av) in expected.items():
    row = summary[name]
    case = root / "production_corrected" / name
    command = json.loads((case / "command.json").read_text())
    options = dict(zip(command[2::2], command[3::2]))
    check(f"{name} kr", math.isclose(float(row["reaction_coefficient_m_s"]), kr), row["reaction_coefficient_m_s"])
    check(f"{name} av", math.isclose(float(row["av_multiplier"]), av), row["av_multiplier"])
    check(f"{name} zero artificial seed", "zero solid" in row["initial_crystal_interpretation"], row["initial_crystal_interpretation"])
    check(f"{name} x0 command", options["-Problem.SaltMoleFraction"] == "0.05", options)
    check(f"{name} x0 accepted field", abs(float(row["initial_salt_mole_fraction"]) - 0.05) <= 1e-12, row["initial_salt_mole_fraction"])
    check(f"{name} activation 0.1 um", float(options["-P2H.ActivationHeight"]) == 1e-7, options)
    check(f"{name} production grid", options["-Grid.Cells1"] == "80" and options["-P2H.CrustCells"] == "64", options)
    check(f"{name} dt", float(options["-TimeLoop.MaxTimeStepSize"]) == 0.5, options)
    check(f"{name} horizon", float(options["-TimeLoop.TEnd"]) == 180000.0, options)
    check(f"{name} salt conservation", float(row["salt_relative_residual"]) <= 1e-8, row["salt_relative_residual"])
    check(f"{name} water conservation", float(row["water_relative_residual"]) <= 1e-8, row["water_relative_residual"])
    transfer = float(row["activation_transfer_relative_residual"])
    activated = math.isfinite(float(row["activation_s"]))
    if activated:
        check(f"{name} activation transfer", math.isfinite(transfer) and transfer <= 1e-12, transfer)
        # A protected stop on the first post-activation trial legitimately leaves
        # only the accepted activation field.  Preserve that failure evidence.
        check(f"{name} crust fields", int(row["crust_field_snapshots"]) >= 1, row["crust_field_snapshots"])
    else:
        check(f"{name} no activation classified", int(row["crust_field_snapshots"]) == 0 and float(row["final_crust_height_m"]) == 0, row)
        expected_status = "FAIL_NO_ALE_ACTIVATION" if int(row["completed_50h"]) else "FAIL_EARLY_STOP"
        check(f"{name} scientific no-go classified", row["scientific_status"] == expected_status, row["scientific_status"])
    check(f"{name} soil fields", int(row["soil_field_snapshots"]) > 2, row["soil_field_snapshots"])
    check(f"{name} field scheduler", float(row["soil_field_max_gap_s"]) <= 100.5, row["soil_field_max_gap_s"])

for stem in (
    "figure1_fixed_av_kr_evolution",
    "figure2_fixed_kr_av_evolution",
    "figure3_evaporation_rate",
    "figure4_height_and_minimum_porosity",
):
    for suffix in ("png", "pdf", "svg"):
        path = analysis / "figures" / f"{stem}.{suffix}"
        check(f"figure {stem}.{suffix}", path.exists() and path.stat().st_size > 1000, path)

output = {"passed": all(item["passed"] for item in checks), "checks": checks}
(analysis / "verification.json").write_text(json.dumps(output, indent=2) + "\n")
print(json.dumps({"passed": output["passed"], "check_count": len(checks)}))
