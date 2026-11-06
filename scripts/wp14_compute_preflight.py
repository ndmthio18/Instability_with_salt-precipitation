"""Compute-node resource and frozen-hash preflight for the one WP14 production run."""
from __future__ import annotations

import hashlib
import json
import os
import shutil
import socket
import sys
from datetime import datetime, timezone
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main(root: Path) -> None:
    job = os.environ.get("SLURM_JOB_ID")
    host = socket.gethostname()
    if not job or host.startswith(("ibcln", "test0", "gln0")):
        raise RuntimeError("preflight must run under Slurm on a compute node")
    run = root / "R48"
    plan = json.loads((run / "plan.json").read_text())
    scratch = Path(os.environ.get("TMPDIR", "/tmp"))
    disk = shutil.disk_usage(scratch)
    mem = int(next(line.split()[1] for line in Path("/proc/meminfo").read_text().splitlines()
                   if line.startswith("MemAvailable:"))) * 1024
    hashes = {
        "executable": sha256(root / "bundle" / "bin" / "paper2_transport"),
        "profile": sha256(root / "eigenprofile.dat"),
        "input": sha256(run / "params.input"),
    }
    expected = {
        "executable": plan["executable_sha256"],
        "profile": plan["profile_sha256"],
        "input": plan["input_sha256"],
    }
    checks = {
        "compute_node": True,
        "one_cpu_affinity": len(os.sched_getaffinity(0)) == 1,
        "scratch_free": disk.free >= plan["minimum_scratch_free_bytes"],
        "memory_available": mem >= plan["minimum_available_memory_bytes"],
        "hashes": hashes == expected,
        "new_root_isolated": "paper2_wp14_dt_refinement_20261106" in str(root),
        "production_not_started": (run / "production_started.json").exists()
                                  and (run / "production_started.json").stat().st_size == 0,
    }
    result = {
        "schema_version": 1,
        "utc": datetime.now(timezone.utc).isoformat(),
        "job_id": job,
        "node": host,
        "affinity": sorted(os.sched_getaffinity(0)),
        "scratch_path": str(scratch),
        "scratch_free_bytes": disk.free,
        "available_memory_bytes": mem,
        "hashes": hashes,
        "checks": checks,
        "passed": all(checks.values()),
    }
    target = run / "compute_preflight.json"
    with target.open("w") as handle:
        json.dump(result, handle, indent=2)
        handle.flush()
        os.fsync(handle.fileno())
    print(json.dumps(result, indent=2))
    if not result["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main(Path(sys.argv[1]))
