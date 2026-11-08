"""Registered P2-H k_r/a_v comparison cases; Slurm compute nodes only."""
import json
import os
import pathlib
import subprocess
import sys
import time

root = pathlib.Path.cwd()
stage = sys.argv[1]
index = int(sys.argv[2]) if len(sys.argv) > 2 else None
if not os.environ.get("SLURM_JOB_ID"):
    raise RuntimeError("run_kr_av.py requires a Slurm allocation")

base = (root / "repo/model/p2_f_analytical_growth/integration/p2f_1d.input").read_text()
base = base.replace("p2f_analytical", "p2h_ale")
base = base.split("[P2F]")[0] + "[Transport]" + base.split("[Transport]", 1)[1]

def case(name, kr, av, ale=64, tend=180000.0):
    return {
        "name": name,
        "reaction_coefficient_m_s": kr,
        "av_multiplier": av,
        "soil_cells": 80,
        "ale_cells": ale,
        "dt_s": 0.5,
        "tend_s": tend,
    }

if stage == "pilot":
    cases = [
        case("pilot_ref_ale64", 1e-3, 1.0, ale=64, tend=20000.0),
        case("pilot_ref_ale128", 1e-3, 1.0, ale=128, tend=20000.0),
    ]
elif stage in ("production", "production_corrected"):
    # The reference belongs to both one-factor groups and is computed once.
    cases = [
        case("kr0p1_av1", 1e-4, 1.0),
        case("kr1_av1_ref", 1e-3, 1.0),
        case("kr10_av1", 1e-2, 1.0),
        case("kr1_av0p5", 1e-3, 0.5),
        case("kr1_av2", 1e-3, 2.0),
    ]
else:
    raise ValueError(f"unknown stage {stage}")

if index is None:
    selected = cases
else:
    if not 0 <= index < len(cases):
        raise IndexError(index)
    selected = [cases[index]]

stage_dir = root / stage
stage_dir.mkdir(exist_ok=True)
manifest_tmp = stage_dir / (
    "cases." + os.environ["SLURM_JOB_ID"] + "." + str(os.getpid()) + ".tmp"
)
manifest_tmp.write_text(json.dumps(cases, indent=2) + "\n")
os.replace(manifest_tmp, stage_dir / "cases.json")

for spec in selected:
    path = stage_dir / spec["name"]
    path.mkdir(exist_ok=False)
    (path / "case_spec.json").write_text(json.dumps(spec, indent=2) + "\n")
    (path / "params.input").write_text(base)
    options = {
        "Grid.Cells1": spec["soil_cells"],
        "TimeLoop.DtInitial": spec["dt_s"],
        "TimeLoop.MaxTimeStepSize": spec["dt_s"],
        "TimeLoop.TEnd": spec["tend_s"],
        "Problem.SaltMoleFraction": 0.05,
        "P2H.ActivationHeight": 1e-7,
        "P2H.DonorMode": "adaptive_multilayer",
        "P2H.MaxDonorDepth": 0.005,
        "P2H.MinimumSupplySaturation": 0.99,
        "P2H.CrustCells": spec["ale_cells"],
        "P2H.ReactionCoefficient": spec["reaction_coefficient_m_s"],
        "P2H.AvMultiplier": spec["av_multiplier"],
        "P2H.PrecipitationMode": "implicit_surface",
        "P2H.FieldOutputInterval": 100.0,
        "P2H.GridDiagnostics": False,
    }
    args = [str(root / "build/soil_solver"), "params.input"]
    for key, value in options.items():
        args += ["-" + key, str(value).lower() if isinstance(value, bool) else str(value)]
    (path / "command.json").write_text(json.dumps(args, indent=2) + "\n")
    started = time.monotonic()
    with (path / "solver.out").open("w") as output:
        try:
            code = subprocess.run(
                args,
                cwd=path,
                stdout=output,
                stderr=subprocess.STDOUT,
                timeout=11 * 3600,
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
    (path / "exit_status.txt").write_text(str(code) + "\n")
    (path / "outcome.json").write_text(json.dumps(outcome, indent=2) + "\n")
    print(json.dumps(outcome), flush=True)
