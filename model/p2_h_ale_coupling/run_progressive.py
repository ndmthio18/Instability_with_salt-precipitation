"""Registered P2-H progressive-activation cases; Slurm compute nodes only."""
import json
import os
import pathlib
import subprocess
import sys
import time


root = pathlib.Path.cwd()
if not os.environ.get("SLURM_JOB_ID"):
    raise RuntimeError("run_progressive.py requires a Slurm compute allocation")

phase = sys.argv[1]
index = int(sys.argv[2])
base = (root / "base.input").read_text()
base = base.replace("Name = p2f_analytical", "Name = p2h_progressive")
base = base.replace("EvaporationBoundaryMode = p2f_analytical", "EvaporationBoundaryMode = p2h_ale")
base = base.split("[P2F]")[0] + "[Transport]" + base.split("[Transport]", 1)[1]


def case(name, soil, dt, mode="progressive_top_boundary_filling", ale=64):
    return {
        "name": name,
        "soil_cells": soil,
        "dt_s": dt,
        "donor_mode": mode,
        "ale_cells": ale,
        "tend_s": 25000.0,
        "salt_mole_fraction": 0.08,
        "bottom_mode": "reservoir",
        "water_level_m": 0.01,
        "activation_height_m": 3.5e-6,
        "reaction_coefficient_m_s": 1e-3,
        "av_multiplier": 1.0,
    }


pilots = [case("pilot_progressive_n80_dt0p5_a64", 80, 0.5)]
matrix = [
    case(f"progressive_n{soil}_dt{str(dt).replace('.', 'p')}_a64", soil, dt)
    for soil in (80, 160, 320) for dt in (0.5, 0.25, 0.125)
] + [
    case("adaptive_n80_dt0p5_a64", 80, 0.5, "adaptive_multilayer"),
    case("toponly_n80_dt0p5_a64", 80, 0.5, "top_only"),
]
ale_matrix = [
    case(f"progressive_n80_dt0p5_a{ale}", 80, 0.5, ale=ale)
    for ale in (128, 256)
]
registry = {"pilot": pilots, "matrix": matrix, "ale": ale_matrix}
if phase not in registry or not 0 <= index < len(registry[phase]):
    raise IndexError((phase, index))
spec = registry[phase][index]

case_dir = root / "runs" / spec["name"]
case_dir.mkdir(parents=True, exist_ok=False)
(case_dir / "case_spec.json").write_text(json.dumps(spec, indent=2) + "\n")
(case_dir / "params.input").write_text(base)
options = {
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
    "P2H.MinimumSupplySaturation": 0.99,
    "P2H.CrustCells": spec["ale_cells"],
    "P2H.ReactionCoefficient": spec["reaction_coefficient_m_s"],
    "P2H.AvMultiplier": spec["av_multiplier"],
    "P2H.PrecipitationMode": "implicit_surface",
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
                              timeout=17.5 * 3600).returncode
    except subprocess.TimeoutExpired:
        code = 124
outcome = {
    "case": spec["name"], "phase": phase, "exit_code": code,
    "elapsed_s": time.monotonic() - started,
    "slurm_job_id": os.environ["SLURM_JOB_ID"], "hostname": os.uname().nodename,
}
(case_dir / "exit_status.txt").write_text(str(code) + "\n")
(case_dir / "outcome.json").write_text(json.dumps(outcome, indent=2) + "\n")
print(json.dumps(outcome), flush=True)
