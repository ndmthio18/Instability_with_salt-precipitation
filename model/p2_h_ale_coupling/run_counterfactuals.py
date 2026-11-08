"""Run standalone ALE concentration counterfactuals inside Slurm."""
import json
import os
import pathlib
import subprocess
import sys
import time


if not os.environ.get("SLURM_JOB_ID"):
    raise RuntimeError("run_counterfactuals.py requires Slurm")
root = pathlib.Path.cwd()
ceq = 0.2649408697676898
actual = 0.2640671415653032
cases = []
for cells in (64, 128, 256):
    cases.append({
        "name": f"actual_mix_cb_eq_a{cells}",
        "cells": cells,
        "c_initial": actual,
        "c_bottom": ceq,
        "mode": "actual_mix",
    })
for delta, label in ((-1e-7, "m1e7"), (-1e-6, "m1e6"), (-5e-6, "m5e6"), (1e-6, "p1e6")):
    cases.append({
        "name": f"actual_mix_cb_{label}_a64",
        "cells": 64,
        "c_initial": actual,
        "c_bottom": ceq + delta,
        "mode": "actual_mix",
    })
cases.append({
    "name": "ceq_reconstructed_cb_eq_a64",
    "cells": 64,
    "c_initial": ceq,
    "c_bottom": ceq,
    "mode": "equilibrium_reconstruction",
})
index = int(sys.argv[1])
spec = cases[index]
out = root / "counterfactual" / spec["name"]
args = [
    str(root / "build" / "interface_counterfactual"),
    str(out), str(spec["cells"]), repr(spec["c_initial"]), repr(spec["c_bottom"]),
    "500", "0.5", spec["mode"],
]
started = time.monotonic()
completed = subprocess.run(args, cwd=root, text=True, capture_output=True)
record = {
    **spec,
    "command": args,
    "exit_code": completed.returncode,
    "stdout": completed.stdout,
    "stderr": completed.stderr,
    "elapsed_s": time.monotonic() - started,
    "slurm_job_id": os.environ["SLURM_JOB_ID"],
    "hostname": os.uname().nodename,
}
out.mkdir(parents=True, exist_ok=True)
(out / "run.json").write_text(json.dumps(record, indent=2) + "\n")
print(json.dumps(record), flush=True)
if completed.returncode:
    raise SystemExit(completed.returncode)
