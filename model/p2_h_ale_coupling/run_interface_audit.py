"""Short, registered P2-H interface-audit cases; execute only in Slurm."""
import json
import os
import pathlib
import subprocess
import sys
import time


root = pathlib.Path.cwd()
if not os.environ.get("SLURM_JOB_ID"):
    raise RuntimeError("run_interface_audit.py requires a Slurm compute allocation")

index = int(sys.argv[1])
base = (root / "base.input").read_text()
base = base.replace("Name = p2f_analytical", "Name = p2h_interface_audit")
base = base.replace("EvaporationBoundaryMode = p2f_analytical", "EvaporationBoundaryMode = p2h_ale")
base = base.split("[P2F]")[0] + "[Transport]" + base.split("[Transport]", 1)[1]


def case(name, x0, bottom, tend, soil=80, ale=64, activation=3.5e-6):
    return {
        "name": name,
        "salt_mole_fraction": x0,
        "bottom_mode": bottom,
        "tend_s": tend,
        "soil_cells": soil,
        "ale_cells": ale,
        "dt_s": 0.5,
        "activation_height_m": activation,
        "reaction_coefficient_m_s": 1e-3,
        "av_multiplier": 1.0,
    }


cases = [
    case("x09999_noflux_n80_a64", 0.09999, "no_flux", 2500.0),
    case("x09999_reservoir_n80_a64", 0.09999, "reservoir", 2500.0),
    case("x008_noflux_n80_a64", 0.08, "no_flux", 25000.0),
    case("x008_reservoir_n80_a64", 0.08, "reservoir", 25000.0),
    case("x008_reservoir_n80_a128", 0.08, "reservoir", 25000.0, ale=128),
    case("x008_reservoir_n80_a256", 0.08, "reservoir", 25000.0, ale=256),
    case("x008_reservoir_n160_a64", 0.08, "reservoir", 25000.0, soil=160),
    case("x008_reservoir_h0p1_n80_a64", 0.08, "reservoir", 21000.0, activation=1e-7),
]
if not 0 <= index < len(cases):
    raise IndexError(index)
spec = cases[index]

case_dir = root / "coupled" / spec["name"]
case_dir.mkdir(parents=True, exist_ok=False)
(case_dir / "case_spec.json").write_text(json.dumps(spec, indent=2) + "\n")
(case_dir / "params.input").write_text(base)
options = {
    "Grid.Cells1": spec["soil_cells"],
    "TimeLoop.DtInitial": spec["dt_s"],
    "TimeLoop.MaxTimeStepSize": spec["dt_s"],
    "TimeLoop.TEnd": spec["tend_s"],
    "Problem.BottomBoundaryMode": spec["bottom_mode"],
    "Problem.WaterLevel": 0.01,
    "Problem.SaltMoleFraction": spec["salt_mole_fraction"],
    "P2H.ActivationHeight": spec["activation_height_m"],
    "P2H.DonorMode": "adaptive_multilayer",
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
        code = subprocess.run(
            args,
            cwd=case_dir,
            stdout=output,
            stderr=subprocess.STDOUT,
            timeout=17.5 * 3600,
        ).returncode
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
