"""Run one registered M1/M2/M3 diagnostic comparison on a Slurm compute node."""
import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import time

if not os.environ.get("SLURM_JOB_ID") or os.uname().nodename.startswith(("ibcln", "login", "test")):
    raise RuntimeError("Slurm compute allocation required")

root = Path.cwd()
src = root / "repo/model/prat_analytical_boundary"
law = sys.argv[1].upper()
dt = float(sys.argv[2]) if len(sys.argv) > 2 else 0.5
if law not in {"M1", "M2", "M3"} or dt not in {0.5, 0.25}:
    raise ValueError((law, dt))

name = f"{law.lower()}_dt{str(dt).replace('.', 'p')}"
case = root / "runs" / name
case.mkdir(parents=True, exist_ok=False)
exe = root / "build/soil_solver"
binary_sha = hashlib.sha256(exe.read_bytes()).hexdigest()
(case / "params.input").write_bytes((src / "column.input").read_bytes())
args = [
    str(exe), "params.input",
    "-Grid.Cells0", "1",
    "-Grid.Cells1", "80",
    "-TimeLoop.TEnd", "9000",
    "-TimeLoop.DtInitial", str(dt),
    "-TimeLoop.MaxTimeStepSize", str(dt),
    "-Problem.SaltMoleFraction", "0.09",
    "-Problem.BottomBoundaryMode", "reservoir",
    "-PratBoundary.ActivationHeight", "3.5e-6",
    "-PratBoundary.ReactionCoefficient", "0.001",
    "-PratBoundary.GrowthLaw", law,
    "-PratBoundary.SubstrateGrainSize", "25e-6",
    "-PratBoundary.MinimumSupplySaturation", "0.99",
]
(case / "command.json").write_text(json.dumps(args, indent=2))
(case / "case_spec.json").write_text(json.dumps({
    "growth_law": law,
    "formula_replacement_diagnostic_only": True,
    "initial_and_reservoir_nacl_mole_fraction": 0.09,
    "soil_vertical_cells": 80,
    "soil_horizontal_cells": 1,
    "dt_max_s": dt,
    "requested_horizon_s": 9000,
    "activation_height_m": 3.5e-6,
    "reaction_coefficient_m_s": 0.001,
    "beta": 2.0,
    "crust_porosity": 0.43,
    "crystal_density_kg_m3": 2165.0,
    "molecular_diffusivity_m2_s": 1.3e-9,
    "substrate_grain_size_m": 25e-6,
    "specific_area_law": "36*(1-epsilon)/d_sub",
    "effective_diffusivity_law": "sqrt(epsilon)*D_s",
    "interface_relative_tolerance": 0.02,
    "binary_sha256": binary_sha,
}, indent=2))

started = time.monotonic()
with (case / "solver.out").open("w") as stream:
    try:
        code = subprocess.run(args, cwd=case, stdout=stream, stderr=subprocess.STDOUT, timeout=3500).returncode
    except subprocess.TimeoutExpired:
        code = 124
elapsed = time.monotonic() - started

rows = list(csv.DictReader((case / "history.csv").open()))
events = list(csv.DictReader((case / "events.csv").open()))
rejected = list(csv.DictReader((case / "rejected_trials.csv").open()))
last = rows[-1]
active = [row for row in rows if row["active"] == "1"]
first_precip = next((float(row["time_s"]) for row in events if row["event"] == "accepted_first_precipitation"), None)
activation = float(active[0]["activation_s"]) if active else None
completed = code == 0 and abs(float(last["time_s"]) - 9000.0) <= 1e-9
result = {
    "name": name,
    "growth_law": law,
    "solver_exit_code": code,
    "completed_requested_horizon": completed,
    "requested_horizon_s": 9000.0,
    "last_accepted_time_s": float(last["time_s"]),
    "first_precipitation_time_s": first_precip,
    "activation_time_s": activation,
    "active_accepted_steps": len(active),
    "last": last,
    "first_rejected_trial": rejected[0] if rejected else None,
    "last_rejected_trial": rejected[-1] if rejected else None,
    "max_salt_residual_relative": max(abs(float(row["salt_residual_rel"])) for row in rows),
    "max_water_residual_relative": max(abs(float(row["water_residual_rel"])) for row in rows),
    "all_m3_steps_valid_da_ge_10": all(row["m3_valid"] == "1" for row in active),
    "binary_sha256": binary_sha,
    "history_sha256": hashlib.sha256((case / "history.csv").read_bytes()).hexdigest(),
    "elapsed_s": elapsed,
    "slurm_job_id": os.environ["SLURM_JOB_ID"],
    "hostname": os.uname().nodename,
}
(case / "outcome.json").write_text(json.dumps(result, indent=2))
print(json.dumps(result, indent=2), flush=True)
