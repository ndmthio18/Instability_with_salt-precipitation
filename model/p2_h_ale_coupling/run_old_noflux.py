"""Run the frozen task-11 executable for a controlled short no-flux regression."""
import json
import os
import pathlib
import subprocess

if not os.environ.get("SLURM_JOB_ID"):
    raise RuntimeError("run_old_noflux.py requires Slurm")
root = pathlib.Path.cwd()
source = root / "verification/no_flux_regression"
target = root / "verification/no_flux_regression_old"
target.mkdir(exist_ok=False)
for name in ("params.input", "case_spec.json"):
    (target / name).write_bytes((source / name).read_bytes())
args = json.loads((source / "command.json").read_text())
args[0] = "/WORK/hulm_work/zhixin/p2_h_kr_av_comparison_20261108/build/soil_solver"
# The frozen binary has no verification-only parameter; it is false anyway.
filtered = []
skip = False
for item in args:
    if skip:
        skip = False
        continue
    if item == "-P2H.VerificationZeroEvaporation":
        skip = True
        continue
    filtered.append(item)
(target / "command.json").write_text(json.dumps(filtered, indent=2) + "\n")
with (target / "solver.out").open("w") as output:
    code = subprocess.run(filtered, cwd=target, stdout=output, stderr=subprocess.STDOUT).returncode
(target / "exit_status.txt").write_text(str(code) + "\n")
raise SystemExit(code)
