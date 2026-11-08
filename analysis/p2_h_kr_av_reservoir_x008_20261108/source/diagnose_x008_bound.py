"""Reproduce the x=0.08 pilot bound stop with read-only diagnostics enabled."""
import csv
import json
import os
import pathlib
import shutil
import subprocess
import time


if not os.environ.get("SLURM_JOB_ID"):
    raise RuntimeError("diagnose_x008_bound.py requires a Slurm compute allocation")

root = pathlib.Path.cwd()
source = root / "pilot/pilot_ref_ale64"
target = root / "diagnostic_bound"
target.mkdir(exist_ok=False)
for name in ("params.input", "case_spec.json"):
    shutil.copy2(source / name, target / name)

args = json.loads((source / "command.json").read_text())


def replace_option(name, value):
    index = args.index(name)
    args[index + 1] = str(value).lower() if isinstance(value, bool) else str(value)


replace_option("-TimeLoop.TEnd", 20582.0)
replace_option("-P2H.GridDiagnostics", True)
(target / "command.json").write_text(json.dumps(args, indent=2) + "\n")

started = time.monotonic()
with (target / "solver.out").open("w") as output:
    code = subprocess.run(
        args,
        cwd=target,
        stdout=output,
        stderr=subprocess.STDOUT,
        timeout=3600,
    ).returncode
outcome = {
    "exit_code": code,
    "elapsed_s": time.monotonic() - started,
    "slurm_job_id": os.environ["SLURM_JOB_ID"],
    "hostname": os.uname().nodename,
}
(target / "outcome.json").write_text(json.dumps(outcome, indent=2) + "\n")

with (target / "ale_diagnostics.csv").open(newline="") as handle:
    records = list(csv.DictReader(handle))


def summarize(rows):
    eps = [float(row["epsilon"]) for row in rows]
    concentration = [float(row["c"]) for row in rows]
    return {
        "time_s": float(rows[0]["time_s"]),
        "status": rows[0]["status"],
        "cells": len(rows),
        "epsilon_min": min(eps),
        "epsilon_min_cell": eps.index(min(eps)),
        "epsilon_max": max(eps),
        "concentration_min": min(concentration),
        "concentration_min_cell": concentration.index(min(concentration)),
        "concentration_max": max(concentration),
        "concentration_max_cell": concentration.index(max(concentration)),
        "epsilon_bound_violations": sum(not (0.02 < value < 0.95) for value in eps),
        "concentration_bound_violations": sum(not (0.0 <= value < 0.95) for value in concentration),
    }


groups = []
current = []
current_key = None
for row in records:
    key = (row["time_s"], row["status"])
    if current and key != current_key:
        groups.append(summarize(current))
        current = []
    current.append(row)
    current_key = key
if current:
    groups.append(summarize(current))

accepted = [group for group in groups if group["status"] == "accepted"]
rejected = [group for group in groups if group["status"] == "bound_rejected"]
result = {
    "run": outcome,
    "expected_solver_failure_reproduced": code != 0,
    "last_accepted": accepted[-1] if accepted else None,
    "first_rejected": rejected[0] if rejected else None,
    "rejected_trial_count": len(rejected),
    "diagnostic_passed": code != 0 and bool(accepted) and bool(rejected),
}
(target / "bound_diagnostic.json").write_text(json.dumps(result, indent=2) + "\n")
print(json.dumps(result, indent=2), flush=True)
raise SystemExit(0 if result["diagnostic_passed"] else 1)
