"""Independent archive-level acceptance for the WP14 dt/2 production run."""
from __future__ import annotations

import csv
import hashlib
import io
import json
import math
import os
import sys
import tarfile
from pathlib import Path


ELF = "112076759fdff0942000307f729c2bd3c7e96963229b43c1b67c1f12bbf9df81"
PROFILE = "d721e605d126e555f92b43200aa73b999ed4e2b61397edd1e98c0d9fddb36575"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def bytes_sha(value: bytes) -> str:
    return hashlib.sha256(value).hexdigest()


def member_bytes(handle: tarfile.TarFile, suffix: str) -> bytes:
    matches = [m for m in handle.getmembers() if m.name == suffix or m.name.endswith("/" + suffix)]
    if len(matches) != 1:
        raise RuntimeError(f"expected one {suffix}, got {[m.name for m in matches]}")
    stream = handle.extractfile(matches[0])
    if stream is None:
        raise RuntimeError(f"could not read {suffix}")
    return stream.read()


def parameter_value(params: bytes, key: str) -> str:
    prefix = key + " ="
    matches = []
    for raw_line in params.decode().splitlines():
        line = raw_line.split("#", 1)[0].strip()
        if line.startswith(prefix):
            matches.append(line.split("=", 1)[1].strip())
    if len(matches) != 1:
        raise RuntimeError(f"expected one {key}, got {len(matches)}")
    return matches[0]


def main(evidence: Path, output: Path) -> None:
    started = json.loads((evidence / "production_started.json").read_text())
    finished = json.loads((evidence / "finished.json").read_text())
    plan = json.loads((evidence / "plan.json").read_text())
    expected_job = (evidence / "submission.txt").read_text().strip()
    fields_path, products_path = evidence / "fields.tar.gz", evidence / "products.tar.gz"
    with tarfile.open(products_path, "r:gz") as products:
        params = member_bytes(products, "params.input")
        profile = member_bytes(products, "eigenprofile.dat")
        stop = json.loads(member_bytes(products, "stop_record.json"))
        checks = json.loads(member_bytes(products, "post/independent_checks.json"))
        audit = json.loads(member_bytes(products, "post/field_audit.json"))
        surface = list(csv.DictReader(io.StringIO(member_bytes(products, "surface_events.csv").decode())))
    with tarfile.open(fields_path, "r:gz") as fields:
        members = fields.getmembers()
        vtu_count = sum(member.name.endswith(".vtu") for member in members)
        pvd_count = sum(member.name.endswith(".pvd") for member in members)

    dt = float(plan["nominal_dt_s"])
    steps = [int(row["step"]) for row in surface]
    times = [float(row["time_s"]) for row in surface]
    dts = [float(row["dt_s"]) for row in surface[1:]]
    nonstandard = [index + 1 for index, value in enumerate(dts) if value != dt]
    terminal_clamp = (
        nonstandard == [len(dts)]
        and 0.0 < dts[-1] <= dt
        and math.isclose(times[-1], float(stop["requested_end_s"]), rel_tol=0.0, abs_tol=1e-8)
        and math.isclose(times[-1] - times[-2], dts[-1], rel_tol=0.0, abs_tol=1e-8)
    )
    expected_frames = 1 + stop["accepted_steps"] // int(plan["output_interval_steps"])
    archive_hashes = {"fields": sha256(fields_path), "products": sha256(products_path)}
    base_salt = float(parameter_value(params, "SaltMoleFraction"))
    amplitude = abs(float(parameter_value(params, "Amplitude")))
    tests = {
        "registered_job": started["job_id"] == expected_job,
        "registered_grid": started["grid"] == plan["grid"] == "R48",
        "four_threads": started["threads"] == 4 and len(started["affinity"]) == 4,
        "executable_hash": started["elf_sha256"] == ELF == plan["executable_sha256"],
        "input_hash": started["input_sha256"] == plan["input_sha256"] == bytes_sha(params),
        "profile_hash": started["profile_sha256"] == PROFILE == bytes_sha(profile),
        "archive_hashes": archive_hashes["fields"] == finished["fields"]["sha256"]
                          and archive_hashes["products"] == finished["products"]["sha256"],
        "worker_complete": finished["complete"] is True and finished["scratch_removed"] is True,
        "solver_success": finished["solver"]["returncode"] == 0 and not finished["solver"]["timed_out"],
        "requested_horizon": stop["reason"] == "requested_horizon"
                             and stop["accepted_steps"] == plan["authorized_steps"] == 14400,
        "surface_rows": len(surface) == stop["accepted_steps"] + 1,
        "step_sequence": steps == list(range(stop["accepted_steps"] + 1)),
        "strictly_increasing_time": all(right > left for left, right in zip(times, times[1:])),
        "fixed_timestep_or_terminal_horizon_clamp": not nonstandard or terminal_clamp,
        "input_dt": parameter_value(params, "DtInitial") == plan["nominal_dt_text"]
                    and parameter_value(params, "MaxTimeStepSize") == plan["nominal_dt_text"],
        "output_interval": int(parameter_value(params, "OutputInterval")) == plan["output_interval_steps"] == 100,
        "initial_range": math.isclose(float(surface[0]["x_domain_min"]), base_salt - amplitude, abs_tol=2e-15)
                         and math.isclose(float(surface[0]["x_domain_max"]), base_salt + amplitude, abs_tol=2e-15),
        "independent_conservation": checks["passed"],
        "field_count": vtu_count == checks["frames"] == expected_frames == plan["expected_frames_if_horizon"],
        "single_pvd": pvd_count == 1,
        "field_grid": audit["cells"] == [48, 200],
        "field_profile": audit["profile_sha256"] == PROFILE,
        "salt_storage_reconstruction": checks["max_inventory_salt_difference_m"] < 1e-14,
        "water_storage_reconstruction": checks["max_inventory_water_difference_m"] < 1e-12,
    }
    result = {
        "schema_version": 1,
        "job_id": os.environ.get("SLURM_JOB_ID"),
        "simulation_job_id": expected_job,
        "tests": tests,
        "passed": all(tests.values()),
        "failed": [key for key, value in tests.items() if not value],
        "archive_hashes": archive_hashes,
        "initial_range_audit": {
            "expected_min": base_salt - amplitude,
            "expected_max": base_salt + amplitude,
            "observed_min": float(surface[0]["x_domain_min"]),
            "observed_max": float(surface[0]["x_domain_max"]),
        },
        "timestep_audit": {
            "expected_dt_text": plan["nominal_dt_text"],
            "expected_dt_s": dt,
            "nonstandard_accepted_steps": nonstandard,
            "terminal_dt_s": dts[-1],
            "terminal_horizon_clamp": terminal_clamp,
        },
        "stop": stop,
        "independent_checks": checks,
        "scope": "numerical/provenance acceptance; not physical convergence or solubility-time validation",
    }
    with output.open("w") as handle:
        json.dump(result, handle, indent=2)
        handle.flush()
        os.fsync(handle.fileno())
    print(json.dumps({"passed": result["passed"], "failed": result["failed"]}, indent=2))
    if not result["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main(Path(sys.argv[1]), Path(sys.argv[2]))
