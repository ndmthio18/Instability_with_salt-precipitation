"""Re-run a failed reservoir pilot with full accepted/rejected ALE fields."""
import csv
import json
import os
import pathlib
import shutil
import subprocess
import time


if not os.environ.get("SLURM_JOB_ID"):
    raise RuntimeError("diagnose_bound_from_pilot.py requires a Slurm compute allocation")

root = pathlib.Path.cwd()
source = root / "pilot/pilot_ref_ale64"
target = root / "diagnostic_bound"
target.mkdir(exist_ok=False)
for name in ("params.input", "case_spec.json"):
    shutil.copy2(source / name, target / name)

with (source / "coupled_history.csv").open(newline="") as handle:
    last = None
    for last in csv.DictReader(handle):
        pass
if last is None:
    raise RuntimeError("failed pilot has no accepted history")
diagnostic_end_s = float(last["time_s"]) + 2.0

args = json.loads((source / "command.json").read_text())


def replace_option(name, value):
    index = args.index(name)
    args[index + 1] = str(value).lower() if isinstance(value, bool) else str(value)


replace_option("-TimeLoop.TEnd", diagnostic_end_s)
replace_option("-P2H.GridDiagnostics", True)
(target / "command.json").write_text(json.dumps(args, indent=2) + "\n")

started = time.monotonic()
with (target / "solver.out").open("w") as output:
    code = subprocess.run(
        args, cwd=target, stdout=output, stderr=subprocess.STDOUT, timeout=3600
    ).returncode
outcome = {
    "exit_code": code,
    "elapsed_s": time.monotonic() - started,
    "slurm_job_id": os.environ["SLURM_JOB_ID"],
    "hostname": os.uname().nodename,
    "diagnostic_end_s": diagnostic_end_s,
}
(target / "outcome.json").write_text(json.dumps(outcome, indent=2) + "\n")


def summarize(rows):
    eps = [float(row["epsilon"]) for row in rows]
    concentration = [float(row["c"]) for row in rows]
    min_eps_index = eps.index(min(eps))
    max_eps_index = eps.index(max(eps))
    return {
        "time_s": float(rows[0]["time_s"]),
        "status": rows[0]["status"],
        "cells": len(rows),
        "epsilon_min": min(eps),
        "epsilon_min_cell": int(float(rows[min_eps_index]["i"])),
        "epsilon_min_z_m": 0.5 * (
            float(rows[min_eps_index]["z_bottom_m"])
            + float(rows[min_eps_index]["z_top_m"])
        ),
        "epsilon_max": max(eps),
        "epsilon_max_cell": int(float(rows[max_eps_index]["i"])),
        "epsilon_max_z_m": 0.5 * (
            float(rows[max_eps_index]["z_bottom_m"])
            + float(rows[max_eps_index]["z_top_m"])
        ),
        "concentration_min": min(concentration),
        "concentration_max": max(concentration),
        "epsilon_bound_violations": sum(not (0.02 < value < 0.95) for value in eps),
        "concentration_bound_violations": sum(
            not (0.0 <= value < 0.95) for value in concentration
        ),
    }


with (target / "ale_diagnostics.csv").open(newline="") as handle:
    records = list(csv.DictReader(handle))
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
