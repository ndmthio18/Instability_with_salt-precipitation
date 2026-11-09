"""Run the registered x=0.08/0.09 flux-mechanism cases on a Slurm compute node."""
import json
import os
import pathlib
import subprocess
import sys
import time


root = pathlib.Path.cwd()
if not os.environ.get("SLURM_JOB_ID"):
    raise RuntimeError("run_flux_x009.py requires a Slurm compute allocation")


def case(name, mole_fraction, dt):
    return {
        "name": name,
        "soil_cells": 80,
        "horizontal_cells": 1,
        "ale_cells": 64,
        "dt_s": dt,
        "tend_s": 25000.0,
        "salt_mole_fraction": mole_fraction,
        "bottom_mode": "reservoir",
        "water_level_m": 0.01,
        "activation_height_m": 3.5e-6,
        "reaction_coefficient_m_s": 1e-3,
        "av_multiplier": 1.0,
        "donor_mode": "progressive_top_boundary_filling",
        "precipitation_mode": "implicit_surface",
        "minimum_supply_saturation": 0.99,
    }


cases = [
    case("x008_progressive_reservoir_n80_a64_dt0p5_diagnostic", 0.08, 0.5),
    case("x009_progressive_reservoir_n80_a64_dt0p5", 0.09, 0.5),
    case("x009_progressive_reservoir_n80_a64_dt0p25_timecheck", 0.09, 0.25),
]
index = int(sys.argv[1])
spec = cases[index]

base = (root / "base.input").read_text()
base = base.replace("Name = p2f_analytical", "Name = p2h_flux_x009")
base = base.replace("EvaporationBoundaryMode = p2f_analytical", "EvaporationBoundaryMode = p2h_ale")
base = base.split("[P2F]")[0] + "[Transport]" + base.split("[Transport]", 1)[1]

case_dir = root / "runs" / spec["name"]
case_dir.mkdir(parents=True, exist_ok=False)
(case_dir / "case_spec.json").write_text(json.dumps(spec, indent=2) + "\n")
(case_dir / "params.input").write_text(base)
options = {
    "Grid.Cells0": spec["horizontal_cells"],
    "Grid.Cells1": spec["soil_cells"],
    "TimeLoop.DtInitial": spec["dt_s"],
    "TimeLoop.MaxTimeStepSize": spec["dt_s"],
    "TimeLoop.TEnd": spec["tend_s"],
    "Problem.BottomBoundaryMode": spec["bottom_mode"],
    "Problem.WaterLevel": spec["water_level_m"],
    "Problem.SaltMoleFraction": spec["salt_mole_fraction"],
    "P2H.ActivationHeight": spec["activation_height_m"],
    "P2H.DonorMode": spec["donor_mode"],
    "P2H.MaxDonorDepth": 0.005,
    "P2H.MinimumSupplySaturation": spec["minimum_supply_saturation"],
    "P2H.CrustCells": spec["ale_cells"],
    "P2H.ReactionCoefficient": spec["reaction_coefficient_m_s"],
    "P2H.AvMultiplier": spec["av_multiplier"],
    "P2H.PrecipitationMode": spec["precipitation_mode"],
    "P2H.FieldOutputInterval": 10.0,
    "P2H.GridDiagnostics": True,
}
args = [str(root / "build" / "soil_solver"), "params.input"]
for key, value in options.items():
    args += ["-" + key, str(value).lower() if isinstance(value, bool) else str(value)]
(case_dir / "command.json").write_text(json.dumps(args, indent=2) + "\n")

started = time.monotonic()
with (case_dir / "solver.out").open("w") as output:
    try:
        code = subprocess.run(args, cwd=case_dir, stdout=output, stderr=subprocess.STDOUT,
                              timeout=3.75 * 3600).returncode
    except subprocess.TimeoutExpired:
        code = 124
outcome = {
    "case": spec["name"],
    "exit_code": code,
    "elapsed_s": time.monotonic() - started,
    "slurm_job_id": os.environ["SLURM_JOB_ID"],
    "hostname": os.uname().nodename,
}
(case_dir / "exit_status.txt").write_text(str(code) + "\n")
(case_dir / "outcome.json").write_text(json.dumps(outcome, indent=2) + "\n")
print(json.dumps(outcome), flush=True)
