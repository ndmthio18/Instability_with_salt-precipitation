"""Prepare the frozen WP12 A=1e-4 three-grid campaign without launching it."""
from __future__ import annotations

import configparser
import hashlib
import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path


HERE = Path(__file__).resolve().parent
PAPER2 = HERE.parent
OUT = PAPER2 / "analysis" / "wp12_amplitude_convergence_20261106"
PROFILE = PAPER2 / "analysis" / "wp10_r24_long_20261104" / "run" / "eigenprofile.dat"
EXECUTABLE = PAPER2 / "analysis" / "wp10_r24_long_20261104" / "run" / "bundle" / "bin" / "paper2_transport"
DT = 113.68588647695695
STEPS = 7200
T_END = 818538.38263409004
T_SCALE = DT * 1200
AMPLITUDE = 1.0e-4
X0 = 0.0036
CAMPAIGN_TAG = "WP12"
CAMPAIGN_NAME = "WP12 enhanced initial perturbation convergence test"
EXPECTED_PROFILE = "d721e605d126e555f92b43200aa73b999ed4e2b61397edd1e98c0d9fddb36575"
EXPECTED_EXECUTABLE = "112076759fdff0942000307f729c2bd3c7e96963229b43c1b67c1f12bbf9df81"

SPECS = {
    "R24": {
        "source": PAPER2 / "analysis" / "wp10_r24_long_20261104" / "run" / "params.input",
        "name": "WP12-R24-A1e4-6T",
        "timeout_s": 9000,
        "wall": "03:00:00",
        "early_seconds_per_step": 0.128739,
        "reference_wall_s": 3742.261341990903,
        "reference_steps": 25502,
        "reference_field_bytes": 139006081,
        "reference_frames": 512,
    },
    "R48": {
        "source": PAPER2 / "analysis" / "wp11_long_validation_20261105" / "R48_run" / "params.input",
        "name": "WP12-R48-A1e4-6T",
        "timeout_s": 14400,
        "wall": "05:00:00",
        "early_seconds_per_step": 0.286219,
        "reference_wall_s": 15817.97,
        "reference_steps": 30000,
        "reference_field_bytes": 317804376,
        "reference_frames": 601,
    },
    "R96": {
        "source": PAPER2 / "analysis" / "wp11_long_validation_20261105" / "inputs" / "R96" / "params.input",
        "name": "WP12-R96-A1e4-6T",
        "timeout_s": 25200,
        "wall": "08:00:00",
        "early_seconds_per_step": 0.686535,
        "reference_wall_s": 32233.45,
        "reference_steps": 30000,
        "reference_field_bytes": 615550172,
        "reference_frames": 601,
    },
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_ini(path: Path) -> dict[str, dict[str, str]]:
    parser = configparser.ConfigParser(interpolation=None)
    parser.optionxform = str
    parser.read(path)
    return {section: dict(parser[section]) for section in parser.sections()}


def replace_unique(text: str, key: str, value: str) -> str:
    pattern = rf"(?m)^({re.escape(key)}\s*=\s*).*$"
    changed, count = re.subn(pattern, rf"\g<1>{value}", text)
    if count != 1:
        raise RuntimeError(f"expected one {key}, got {count}")
    return changed


def amplitude_text() -> str:
    return f"{AMPLITUDE:.0e}".replace("e-0", "e-").replace("e+0", "e+")


def amplitude_name() -> str:
    return amplitude_text().replace("-", "")


def job_script(grid: str, wall: str) -> str:
    return f"""#!/bin/bash
#SBATCH --job-name={CAMPAIGN_TAG}-{grid}-A{amplitude_name()}
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time={wall}
#SBATCH --output=production.stdout
#SBATCH --error=production.stderr
set -euo pipefail
export PATH=/rmprog/slurm/v22.05.7/bin:/usr/local/bin:/usr/bin:/bin
case "$(hostname)" in ibcln*|test0*|gln0*) echo "refusing login/test node"; exit 2;; esac
root="$1"; export WP11_THREADS=4
srun -n 1 -c 4 /WORK/hulm_work/ys/miniconda3/bin/python3 "$root/wp11_hpcc_worker.py" "$root" {grid}
"""


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    profile_hash = sha256(PROFILE)
    executable_hash = sha256(EXECUTABLE)
    if profile_hash != EXPECTED_PROFILE or executable_hash != EXPECTED_EXECUTABLE:
        raise RuntimeError("profile or executable hash drift")

    profile_values = []
    for line in PROFILE.read_text().splitlines():
        fields = line.split()
        if len(fields) >= 2:
            try:
                profile_values.append(float(fields[1]))
            except ValueError:
                continue
    profile_min, profile_max = min(profile_values), max(profile_values)

    comparisons = []
    estimates = []
    for grid, spec in SPECS.items():
        target_dir = OUT / "inputs" / grid
        target_dir.mkdir(parents=True, exist_ok=True)
        source_text = spec["source"].read_text()
        text = replace_unique(source_text, "TEnd", format(T_END, ".17g"))
        text = replace_unique(text, "Name", spec["name"])
        text = replace_unique(text, "Amplitude", amplitude_text())
        target = target_dir / "params.input"
        target.write_text(text)

        before, after = load_ini(spec["source"]), load_ini(target)
        differing = []
        for section in sorted(set(before) | set(after)):
            for key in sorted(set(before.get(section, {})) | set(after.get(section, {}))):
                if before.get(section, {}).get(key) != after.get(section, {}).get(key):
                    differing.append(f"{section}.{key}")
        allowed = ["Problem.Amplitude", "Problem.Name", "TimeLoop.TEnd"]
        if differing != allowed:
            raise RuntimeError(f"{grid}: unexpected differences {differing}")

        nx = int(after["Grid"]["Cells0"])
        width = float(after["Grid"]["Positions0"].split()[1])
        sampled = [
            X0 + AMPLITUDE * aval * math.cos(2.0 * math.pi * (width * i / nx) / width)
            for aval in profile_values for i in range(nx + 1)
        ]
        initial_range = [min(sampled), max(sampled)]
        if not (math.isclose(initial_range[0], X0 - AMPLITUDE, abs_tol=2e-15) and math.isclose(initial_range[1], X0 + AMPLITUDE, abs_tol=2e-15)):
            raise RuntimeError(f"{grid}: unexpected initial range {initial_range}")

        plan = {
            "schema_version": 1,
            "campaign": CAMPAIGN_NAME,
            "grid": grid,
            "run_name": spec["name"],
            "threads": 4,
            "input_sha256": sha256(target),
            "profile_sha256": profile_hash,
            "executable_sha256": executable_hash,
            "solver_timeout_s": spec["timeout_s"],
            "minimum_scratch_free_bytes": 32 * 1024**3,
            "minimum_available_memory_bytes": 16 * 1024**3,
            "authorized_steps": STEPS,
            "authorized_horizon_T": 6.0,
            "authorized_end_s": T_END,
            "authorized_amplitude": AMPLITUDE,
            "initial_salt_mole_fraction_range": initial_range,
            "allowed_input_differences_from_reference": allowed,
            "reference_input": str(spec["source"].relative_to(PAPER2)),
            "reference_input_sha256": sha256(spec["source"]),
        }
        (target_dir / "plan.json").write_text(json.dumps(plan, indent=2))
        (target_dir / "job.sbatch").write_text(job_script(grid, spec["wall"]), newline="\n")
        comparisons.append({"grid": grid, "differences": differing, "initial_range": initial_range, "input_sha256": plan["input_sha256"]})
        expected_frames = 1 + STEPS // 50
        estimates.append({
            "grid": grid,
            "early_linear_wall_s": spec["early_seconds_per_step"] * STEPS,
            "long_run_linear_wall_s": spec["reference_wall_s"] * STEPS / spec["reference_steps"],
            "allocated_wall": spec["wall"],
            "expected_frames": expected_frames,
            "estimated_fields_archive_bytes": round(spec["reference_field_bytes"] * expected_frames / spec["reference_frames"]),
        })

    protocol = {
        "schema_version": 1,
        "frozen_utc": datetime.now(timezone.utc).isoformat(),
        "scope": f"three grids at amplitude {amplitude_text()} through at most 6T; no extension is authorized",
        "fixed": {
            "case": 36877,
            "amplitude": AMPLITUDE,
            "reference_amplitude": 1.0e-6,
            "profile": "A28",
            "profile_sha256": profile_hash,
            "physical_wavelength_m": 0.021428571428571432,
            "domain_m": [0.021428571428571432, 0.2],
            "vertical_cells": 200,
            "horizontal_grids": [24, 48, 96],
            "dt_s": DT,
            "steps": STEPS,
            "end_s": T_END,
            "end_T": 6.0,
            "output_interval_steps": 50,
            "stop_at_surface_solubility": True,
        },
        "source_interpretation": "Problem.Amplitude is an additive salt-mole-fraction coefficient multiplying A28(y)*cos(2*pi*x/lambda); it does not alter wavelength or profile.",
        "comparison_axes": {
            "amplitude_sensitivity": f"1e-6 versus {amplitude_text()} on the same grid only",
            "grid_sensitivity": f"R24/R48/R96 at {amplitude_text()} only",
        },
        "registered_metrics": {
            "nonlinear_harmonic_signal": "H23=sqrt(harmonic2_top^2+harmonic3_top^2)",
            "nonlinear_milestone": "first accepted step of 50 consecutive accepted steps with H23 >= 1e-5",
            "sustained_positive_advective_flux": "first accepted step of 50 consecutive accepted steps with conservative SCVF net advective salt flux > 1e-15 m/s",
            "local_downward_velocity": "first two consecutive saved fields with max downward water velocity > 1e-10 m/s and downward-area fraction > 0.001",
            "mirror_error": "L2(x(x,z)-x(W-x,z))/L2(x(x,z)-horizontal_mean(z))",
            "odd_mode_fraction": "sum odd nonzero horizontal Fourier energy / sum all nonzero horizontal Fourier energy",
            "milestone_alignment": "plot against t-t_nl only when the registered milestone is observed; never optimize shifts",
        },
        "claim_boundary": "accelerated nonlinear convergence benchmark; not a replacement for A=1e-6 and not validation of solubility-crossing time",
    }
    (OUT / "protocol.json").write_text(json.dumps(protocol, indent=2))
    preflight = {
        "schema_version": 1,
        "profile": {"minimum": profile_min, "maximum": profile_max, "sha256": profile_hash},
        "executable_sha256": executable_hash,
        "initial_range_checks": comparisons,
        "cost_and_storage_estimates": estimates,
        "estimated_total_fields_archive_bytes": sum(item["estimated_fields_archive_bytes"] for item in estimates),
        "authorization": f"user explicitly authorized R24/R48/R96, A={amplitude_text()}, maximum 6T",
    }
    (OUT / "preflight_local.json").write_text(json.dumps(preflight, indent=2))
    print(json.dumps(preflight, indent=2))


if __name__ == "__main__":
    main()
