"""Registered reservoir-fed P2-H cases; execute only inside Slurm."""
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
    raise RuntimeError("run_reservoir.py requires a Slurm compute allocation")
salt_mole_fraction = float(os.environ.get("P2H_SALT_MOLE_FRACTION", "0.05"))
if not 0.0 <= salt_mole_fraction < 1.0:
    raise ValueError("P2H_SALT_MOLE_FRACTION must be in [0, 1)")

base = (root / "repo/model/p2_f_analytical_growth/integration/p2f_1d.input").read_text()
base = base.replace("Name = p2f_analytical", "Name = p2h_ale_reservoir")
base = base.replace("EvaporationBoundaryMode = p2f_analytical", "EvaporationBoundaryMode = p2h_ale")
base = base.split("[P2F]")[0] + "[Transport]" + base.split("[Transport]", 1)[1]


def case(name, kr=1e-3, av=1.0, ale=64, tend=180000.0, dt=0.5, mode="reservoir", zero_evap=False):
    return {
        "name": name,
        "reaction_coefficient_m_s": kr,
        "av_multiplier": av,
        "soil_cells": 80,
        "ale_cells": ale,
        "dt_s": dt,
        "tend_s": tend,
        "bottom_mode": mode,
        "salt_mole_fraction": salt_mole_fraction,
        "zero_evaporation_verification": zero_evap,
    }


if stage == "verification":
    cases = [
        case("hydrostatic_reservoir", tend=100.0, zero_evap=True),
        case("evaporating_reservoir", tend=300.0),
        case("no_flux_regression", tend=300.0, mode="no_flux"),
    ]
elif stage == "pilot":
    cases = [case("pilot_ref_ale64")]
elif stage == "ale_pilot":
    with (root / "pilot/pilot_ref_ale64/coupled_history.csv").open() as handle:
        header = handle.readline().strip().split(",")
        active_index = header.index("active")
        activation_index = header.index("activation_s")
        activation = None
        for line in handle:
            values = line.strip().split(",")
            if int(values[active_index]) == 1:
                activation = float(values[activation_index])
                break
    if activation is None:
        raise RuntimeError("ALE pilot forbidden: reference pilot did not activate")
    gate_end = min(180000.0, activation + 3600.0)
    cases = [
        case("pilot_grid_ale128", tend=gate_end, ale=128),
        case("pilot_dt_half_ale64", tend=gate_end, ale=64, dt=0.25),
    ]
elif stage == "matrix":
    cases = [
        case("kr0p1_av1", kr=1e-4),
        case("kr1_av1_ref"),
        case("kr10_av1", kr=1e-2),
        case("kr1_av0p5", av=0.5),
        case("kr1_av2", av=2.0),
    ]
else:
    raise ValueError(f"unknown stage {stage}")

selected = cases if index is None else [cases[index]]
stage_dir = root / stage
stage_dir.mkdir(exist_ok=True)
manifest = stage_dir / "cases.json"
manifest.write_text(json.dumps(cases, indent=2) + "\n")

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
        "Problem.BottomBoundaryMode": spec["bottom_mode"],
        "Problem.WaterLevel": 0.01,
        "Problem.SaltMoleFraction": spec["salt_mole_fraction"],
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
        "P2H.VerificationZeroEvaporation": spec["zero_evaporation_verification"],
    }
    args = [str(root / "build/soil_solver"), "params.input"]
    for key, value in options.items():
        args += ["-" + key, str(value).lower() if isinstance(value, bool) else str(value)]
    (path / "command.json").write_text(json.dumps(args, indent=2) + "\n")
    started = time.monotonic()
    with (path / "solver.out").open("w") as output:
        try:
            code = subprocess.run(args, cwd=path, stdout=output, stderr=subprocess.STDOUT, timeout=18 * 3600).returncode
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
