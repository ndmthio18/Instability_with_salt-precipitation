"""Run the zero-source saturated-interface compatibility pilot on Slurm."""
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import time

if not os.environ.get("SLURM_JOB_ID"):
    raise RuntimeError("Slurm compute allocation required")
root = Path.cwd()
baseline = Path("/WORK/hulm_work/zhixin/p2_h_x009_no_dissolution_20261109/runs/x009_progressive_n80_a64_dt0p5_no_dissolution")
case = root / "runs" / "saturated_no_source"
case.mkdir(parents=True, exist_ok=False)
original = json.loads((baseline / "command.json").read_text())
args = list(original)
args[0] = str(root / "build" / "soil_solver")
args[args.index("-TimeLoop.TEnd") + 1] = "6750"
args += ["-P2H.InterfaceMode", "saturated_no_source_strict"]
assert args[args.index("-P2H.SuppressDissolution") + 1] == "true"
assert args[args.index("-Problem.SaltMoleFraction") + 1] == "0.09"
for name in ["params.input"]:
    (case / name).write_bytes((baseline / name).read_bytes())
(case / "command.json").write_text(json.dumps(args, indent=2))
(case / "source_command.json").write_text(json.dumps(original, indent=2))
spec = json.loads((baseline / "case_spec.json").read_text())
spec.update(name="saturated_no_source", requested_horizon_s=6750,
            interface_mode="saturated_no_source_strict", suppress_dissolution=True,
            interface_phase_source_kg_m2_s=0, expected_outcome="compatibility test")
(case / "case_spec.json").write_text(json.dumps(spec, indent=2))
started = time.monotonic()
with (case / "solver.out").open("w") as stream:
    completed = subprocess.run(args, cwd=case, stdout=stream, stderr=subprocess.STDOUT)
elapsed = time.monotonic() - started

with (case / "coupled_history.csv").open() as stream:
    history = list(csv.DictReader(stream))
with (case / "interface_compatibility.csv").open() as stream:
    attempts = list(csv.DictReader(stream))
with (case / "coupling_events.csv").open() as stream:
    events = list(csv.DictReader(stream))
active = [row for row in history if int(row["active"])]
assert active and float(active[0]["time_s"]) == 6745.0
assert float(history[-1]["time_s"]) == 6745.0
assert completed.returncode != 0
assert len(attempts) == 1 and int(attempts[0]["accepted"]) == 0
assert float(attempts[0]["interface_phase_source_kg_m2_s"]) == 0
assert float(active[-1]["internal_dissolution_kg_m2"]) == 0
assert any("saturated no-source interface bound exceeded" in row["event"] for row in events)
salt0 = float(history[0]["soil_salt_kg"])
water0 = float(history[0]["soil_water_kg"])
salt_relative = abs(float(history[-1]["salt_residual_kg"])) / salt0
water_relative = abs(float(history[-1]["water_residual_kg"])) / water0
assert max(salt_relative, water_relative) <= 1e-8
result = {
    "verdict": "INCOMPATIBLE_WITHOUT_SOURCE",
    "solver_exit_code_expected_nonzero": completed.returncode,
    "last_accepted_time_s": float(history[-1]["time_s"]),
    "attempted_time_s": float(attempts[0]["time_s"]),
    "attempted_dt_s": float(attempts[0]["dt_s"]),
    "ale_boundary_c": float(attempts[0]["ale_boundary_c"]),
    "soil_surface_c": float(attempts[0]["soil_surface_c"]),
    "concentration_residual": float(attempts[0]["residual_c"]),
    "soil_surface_sw": float(attempts[0]["soil_surface_sw"]),
    "ale_liquid_kg_m2_s": float(attempts[0]["ale_liquid_kg_m2_s"]),
    "ale_salt_kg_m2_s": float(attempts[0]["ale_salt_kg_m2_s"]),
    "interface_phase_source_kg_m2_s": 0.0,
    "internal_dissolution_kg_m2": 0.0,
    "accepted_salt_relative_residual": salt_relative,
    "accepted_water_relative_residual": water_relative,
    "elapsed_s": elapsed,
    "slurm_job_id": os.environ["SLURM_JOB_ID"],
    "hostname": os.uname().nodename,
}
(case / "result.json").write_text(json.dumps(result, indent=2))
(case / "outcome.json").write_text(json.dumps({"diagnostic_success": True, **result}, indent=2))
print(json.dumps(result, indent=2), flush=True)

