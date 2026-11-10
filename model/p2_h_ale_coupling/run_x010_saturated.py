"""Run the recorded 80-layer Marc-ALE case from saturated NaCl pore liquid.

The verified P2-H executable and the x=0.09 reference command are reused.
Only the initial/reservoir NaCl mole fraction and requested horizon are changed.
The Marc finite-rate top closure, reversible ALE internal reaction, normal soil
diffusion, 3.5 um activation threshold, grid, and time-step controls remain the
recorded reference values.
"""

import argparse
import csv
import hashlib
import json
import os
import pathlib
import subprocess
import time


BASE = pathlib.Path(
    "/WORK/hulm_work/zhixin/p2_h_flux_mechanism_x009_20261109/"
    "runs/x009_progressive_reservoir_n80_a64_dt0p5"
)
EXPECTED_BINARY_SHA256 = (
    "32f2c0dd1f6265859489f9c882b899bb6119a0c6571c6127e92f475c75498a93"
)
SATURATED_MOLE_FRACTION = 0.10


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def set_option(args, key, value):
    flag = "-" + key
    if flag in args:
        args[args.index(flag) + 1] = str(value)
    else:
        args.extend([flag, str(value)])


parser = argparse.ArgumentParser()
parser.add_argument("--name", required=True)
parser.add_argument("--horizon", required=True, type=float)
options = parser.parse_args()

if not os.environ.get("SLURM_JOB_ID"):
    raise RuntimeError("Slurm compute allocation required")
if os.uname().nodename.startswith(("ibcln", "login", "test")):
    raise RuntimeError("Refusing to calculate on a login or test node")

root = pathlib.Path.cwd()
original = json.loads((BASE / "command.json").read_text())
binary = pathlib.Path(original[0])
if sha256(binary) != EXPECTED_BINARY_SHA256:
    raise RuntimeError("Verified soil_solver identity changed")

args = list(original)
set_option(args, "TimeLoop.TEnd", options.horizon)
set_option(args, "Problem.SaltMoleFraction", SATURATED_MOLE_FRACTION)

# Encode the intended single-factor comparison explicitly.
assert args[args.index("-Grid.Cells1") + 1] == "80"
assert args[args.index("-TimeLoop.MaxTimeStepSize") + 1] == "0.5"
assert args[args.index("-Problem.BottomBoundaryMode") + 1] == "reservoir"
assert args[args.index("-Problem.SaltMoleFraction") + 1] == "0.1"
assert args[args.index("-P2H.ActivationHeight") + 1] == "3.5e-06"
assert args[args.index("-P2H.CrustCells") + 1] == "64"
assert args[args.index("-P2H.ReactionCoefficient") + 1] == "0.001"
assert "-Problem.SaltDiffusivity" not in args
assert "-P2H.TopBoundaryMode" not in args
assert "-P2H.SuppressDissolution" not in args

case = root / "runs" / options.name
case.mkdir(parents=True, exist_ok=False)
(case / "params.input").write_bytes((BASE / "params.input").read_bytes())

spec = json.loads((BASE / "case_spec.json").read_text())
spec.update(
    name=options.name,
    tend_s=options.horizon,
    salt_mole_fraction=SATURATED_MOLE_FRACTION,
    reservoir_salt_mole_fraction=SATURATED_MOLE_FRACTION,
    soil_salt_diffusivity="unchanged_reference_value",
    ale_diffusion="unchanged_verified_binary",
    top_boundary_closure="Marc finite kinetics (verified binary default)",
    interface_concentration="free evolution (no imposed C_eq)",
    internal_reaction="reversible precipitation/dissolution",
    bottom_solid_dissolution_source=False,
    comparison="x=0.09 reference with only initial/reservoir salinity changed",
)
(case / "case_spec.json").write_text(json.dumps(spec, indent=2) + "\n")
(case / "command.json").write_text(json.dumps(args, indent=2) + "\n")
(case / "source_command.json").write_text(json.dumps(original, indent=2) + "\n")
(case / "binary_sha256.txt").write_text(
    EXPECTED_BINARY_SHA256 + "  " + str(binary) + "\n"
)

started = time.monotonic()
with (case / "solver.out").open("w") as output:
    try:
        code = subprocess.run(
            args,
            cwd=case,
            stdout=output,
            stderr=subprocess.STDOUT,
            timeout=3.5 * 3600,
        ).returncode
    except subprocess.TimeoutExpired:
        code = 124

outcome = {
    "case": options.name,
    "requested_horizon_s": options.horizon,
    "exit_code": code,
    "elapsed_s": time.monotonic() - started,
    "slurm_job_id": os.environ["SLURM_JOB_ID"],
    "hostname": os.uname().nodename,
}

history_path = case / "coupled_history.csv"
if history_path.exists():
    with history_path.open() as stream:
        history = list(csv.DictReader(stream))
    if history:
        outcome["last_accepted_time_s"] = float(history[-1]["time_s"])
        active = [row for row in history if int(row["active"])]
        outcome["activation_time_s"] = (
            float(active[0]["time_s"]) if active else None
        )
        outcome["internal_positive_kg_m2"] = (
            float(active[-1]["internal_positive_kg_m2"]) if active else 0.0
        )
        outcome["internal_dissolution_kg_m2"] = (
            float(active[-1]["internal_dissolution_kg_m2"]) if active else 0.0
        )

(case / "exit_status.txt").write_text(str(code) + "\n")
(case / "outcome.json").write_text(json.dumps(outcome, indent=2) + "\n")
print(json.dumps(outcome), flush=True)
raise SystemExit(code)
