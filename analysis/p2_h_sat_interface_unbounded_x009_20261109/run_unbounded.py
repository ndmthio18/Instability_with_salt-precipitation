"""Run the no-dissolution saturated-crust control without a jump stop."""
import csv
import json
import os
from pathlib import Path
import subprocess
import time

if not os.environ.get("SLURM_JOB_ID"):
    raise RuntimeError("Slurm compute allocation required")

root = Path.cwd()
baseline = Path("/WORK/hulm_work/zhixin/p2_h_x009_no_dissolution_20261109/runs/x009_progressive_n80_a64_dt0p5_no_dissolution")
case = root / "runs" / "crust_saturated_soil_unbounded"
case.mkdir(parents=True, exist_ok=False)
original = json.loads((baseline / "command.json").read_text())
args = list(original)
args[0] = str(root / "build" / "soil_solver")
args[args.index("-TimeLoop.TEnd") + 1] = "12500"
args += [
    "-P2H.InterfaceMode", "crust_saturated_soil_undersaturated",
    "-P2H.EnforceSoilSideUndersaturationBound", "false",
]
assert args[args.index("-P2H.SuppressDissolution") + 1] == "true"
assert args[args.index("-Problem.SaltMoleFraction") + 1] == "0.09"
assert args[args.index("-Grid.Cells1") + 1] == "80"
(case / "params.input").write_bytes((baseline / "params.input").read_bytes())
(case / "command.json").write_text(json.dumps(args, indent=2))
(case / "source_command.json").write_text(json.dumps(original, indent=2))
spec = json.loads((baseline / "case_spec.json").read_text())
spec.update(
    name="crust_saturated_soil_unbounded",
    requested_horizon_s=12500,
    interface_mode="crust_saturated_soil_undersaturated",
    enforce_soil_side_undersaturation_bound=False,
    suppress_dissolution=True,
    interface_phase_source_kg_m2_s=0,
    physical_status="unbounded concentration-jump diagnostic",
)
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
accepted_attempts = [row for row in attempts if int(row["accepted"])]
assert active and float(active[0]["time_s"]) == 6745.0
assert attempts and all(int(row["enforce_undersaturation_bound"]) == 0 for row in attempts)
assert all(float(row["interface_phase_source_kg_m2_s"]) == 0 for row in attempts)
assert all(float(row["internal_dissolution_kg_m2"]) == 0 for row in active)
salt0 = float(history[0]["soil_salt_kg"])
water0 = float(history[0]["soil_water_kg"])
salt_relative = max(abs(float(row["salt_residual_kg"])) / salt0 for row in history)
water_relative = max(abs(float(row["water_residual_kg"])) / water0 for row in history)
assert max(salt_relative, water_relative) <= 1e-8
last = history[-1]
reached_horizon = completed.returncode == 0 and abs(float(last["time_s"]) - 12500.0) <= 1e-9
result = {
    "verdict": "COMPLETED_REQUESTED_HORIZON" if reached_horizon else "STOPPED_BY_MODEL_OR_SOLVER",
    "solver_exit_code": completed.returncode,
    "requested_horizon_s": 12500.0,
    "last_accepted_time_s": float(last["time_s"]),
    "activation_time_s": float(active[0]["time_s"]),
    "accepted_interface_steps": len(accepted_attempts),
    "minimum_soil_minus_ceq": min(float(row["residual_c"]) for row in accepted_attempts),
    "maximum_abs_soil_minus_ceq": max(abs(float(row["residual_c"])) for row in accepted_attempts),
    "final_soil_surface_c": float(accepted_attempts[-1]["soil_surface_c"]),
    "ale_boundary_c": float(accepted_attempts[-1]["ale_boundary_c"]),
    "final_soil_surface_sw": float(accepted_attempts[-1]["soil_surface_sw"]),
    "final_height_um": float(last["h_m"]) * 1e6,
    "final_internal_positive_g_m2": float(last["internal_positive_kg_m2"]) * 1e3,
    "final_top_new_g_m2": float(last["top_new_solid_kg_m2"]) * 1e3,
    "internal_dissolution_kg_m2": 0.0,
    "interface_phase_source_kg_m2_s": 0.0,
    "salt_max_relative_residual": salt_relative,
    "water_max_relative_residual": water_relative,
    "last_event": events[-1]["event"] if events else None,
    "elapsed_s": elapsed,
    "slurm_job_id": os.environ["SLURM_JOB_ID"],
    "hostname": os.uname().nodename,
}
(case / "result.json").write_text(json.dumps(result, indent=2))
(case / "outcome.json").write_text(json.dumps({"diagnostic_success": True, **result}, indent=2))
print(json.dumps(result, indent=2), flush=True)
