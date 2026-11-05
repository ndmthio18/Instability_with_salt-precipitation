"""Independent, read-only acceptance gate for one completed WP11 leg."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

EXPECTED_ELF = "112076759fdff0942000307f729c2bd3c7e96963229b43c1b67c1f12bbf9df81"
EXPECTED_PROFILE = "d721e605d126e555f92b43200aa73b999ed4e2b61397edd1e98c0d9fddb36575"
EXPECTED_DT = 113.68588647695695


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_csv(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle))


def prefix_equal(old, new):
    old_lines = Path(old).read_bytes().splitlines()
    new_lines = Path(new).read_bytes().splitlines()
    return {"old_lines": len(old_lines), "new_lines": len(new_lines), "bitwise_prefix": new_lines[:len(old_lines)] == old_lines}


def finite_csv(path):
    rows = read_csv(path)
    for row in rows:
        for value in row.values():
            try:
                number = float(value)
            except (TypeError, ValueError):
                continue
            if not math.isfinite(number):
                return False
    return True


def main(evidence, run, out, label, expected_input, expected_job, prefix_run=None):
    evidence, run, out = Path(evidence), Path(run), Path(out)
    out.mkdir(parents=True, exist_ok=True)
    started = json.loads((evidence / "production_started.json").read_text())
    finished = json.loads((evidence / "finished.json").read_text())
    stop = json.loads((run / "stop_record.json").read_text())
    checks = json.loads((run / "post" / "independent_checks.json").read_text())
    audit = json.loads((run / "post" / "field_audit.json").read_text())
    surface = read_csv(run / "surface_events.csv")
    steps = [int(row["step"]) for row in surface]
    times = [float(row["time_s"]) for row in surface]
    accepted_dt = [float(row["dt_s"]) for row in surface[1:]]
    archive_tests = {
        "fields_archive": sha256(evidence / "fields.tar.gz") == finished["fields"]["sha256"],
        "products_archive": sha256(evidence / "products.tar.gz") == finished["products"]["sha256"],
    }
    prefix = None
    if prefix_run:
        prefix_run = Path(prefix_run)
        prefix = {name: prefix_equal(prefix_run / name, run / name) for name in
                  ("surface_events.csv", "weighted_online.csv", "transport_budget.csv")}
    tests = {
        "registered_job": started["job_id"] == expected_job,
        "registered_grid": started["grid"] == label,
        "four_threads": started["threads"] == 4 and len(started["affinity"]) == 4,
        "executable_hash": started["elf_sha256"] == EXPECTED_ELF,
        "input_hash_started": started["input_sha256"] == expected_input,
        "input_hash_local": sha256(run / "params.input") == expected_input,
        "profile_hash_started": started["profile_sha256"] == EXPECTED_PROFILE,
        "profile_hash_local": sha256(run / "eigenprofile.dat") == EXPECTED_PROFILE,
        "worker_complete": finished["complete"] is True and finished["scratch_removed"] is True,
        "solver_success": finished["solver"]["returncode"] == 0 and not finished["solver"]["timed_out"],
        "registered_stop": stop["reason"] in ("requested_horizon", "surface_solubility_threshold", "interior_solubility_threshold"),
        "accepted_step_count": stop["accepted_steps"] == len(surface) - 1,
        "step_sequence": steps == list(range(len(steps))),
        "strictly_increasing_time": all(b > a for a, b in zip(times, times[1:])),
        "no_unregistered_timestep_reduction": bool(accepted_dt) and all(value == EXPECTED_DT for value in accepted_dt),
        "independent_conservation": checks["passed"],
        "salt_storage_reconstruction": checks["max_inventory_salt_difference_m"] < 1e-14,
        "water_storage_reconstruction": checks["max_inventory_water_difference_m"] < 1e-12,
        "cumulative_flux_reintegration": all(
            item["cumulative_salt_reintegration_error"] == 0.0 and item["cumulative_water_reintegration_error"] == 0.0
            for item in checks["independent_budgets"].values()),
        "terminal_field_preserved": checks["terminal_irregular_stride"] and checks["frames"] == len(list(run.glob("*.vtu"))),
        "field_audit_grid": audit["cells"][0] == int(label[1:]) and audit["cells"][1] == 200,
        "field_audit_profile": audit["profile_sha256"] == EXPECTED_PROFILE,
        "finite_surface": finite_csv(run / "surface_events.csv"),
        "finite_online": finite_csv(run / "weighted_online.csv"),
        "finite_budget": finite_csv(run / "transport_budget.csv"),
        **archive_tests,
    }
    if prefix is not None:
        tests["registered_short_trajectory_prefix"] = all(item["bitwise_prefix"] for item in prefix.values())
    result = {
        "schema_version": 1,
        "grid": label,
        "job_id": expected_job,
        "tests": tests,
        "passed": all(tests.values()),
        "stop": stop,
        "surface_or_interior_event": stop["reason"],
        "prefix_comparison": prefix,
        "independent_checks": checks,
        "scope": "numerical/provenance acceptance; does not establish physical convergence or precipitation kinetics",
    }
    (out / "acceptance.json").write_text(json.dumps(result, indent=2))
    print(json.dumps({"grid": label, "passed": result["passed"], "failed": [key for key, value in tests.items() if not value]}, indent=2))
    if not result["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("evidence")
    parser.add_argument("run")
    parser.add_argument("out")
    parser.add_argument("label")
    parser.add_argument("expected_input")
    parser.add_argument("expected_job")
    parser.add_argument("--prefix-run")
    args = parser.parse_args()
    main(args.evidence, args.run, args.out, args.label, args.expected_input, args.expected_job, args.prefix_run)
