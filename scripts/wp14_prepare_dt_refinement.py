"""Prepare the single authorized WP14 R48 A=1e-4 dt/2 production input."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path


PAPER2 = Path(__file__).resolve().parents[1]
REFERENCE = PAPER2 / "analysis" / "wp12_amplitude_convergence_20261106" / "inputs" / "R48" / "params.input"
OUT = PAPER2 / "analysis" / "wp14_nonlinear_timestep_20261106" / "phase_b" / "inputs"
REFERENCE_SHA = "ea7993e1b65fc5329fe1da9d3725641bfcdf6d0a715c8b227e469fcb9f986ca2"
ELF_SHA = "112076759fdff0942000307f729c2bd3c7e96963229b43c1b67c1f12bbf9df81"
PROFILE_SHA = "d721e605d126e555f92b43200aa73b999ed4e2b61397edd1e98c0d9fddb36575"
DT_TEXT = "56.842943238478475"
DT = float(DT_TEXT)
T_END = 818538.38263409
STEPS = 14400


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def replace_unique(text: str, old: str, new: str) -> str:
    if text.count(old) != 1:
        raise AssertionError(f"expected one occurrence: {old!r}")
    return text.replace(old, new)


def parse_ini(text: str) -> dict[str, str]:
    section = ""
    result = {}
    for raw in text.splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if line.startswith("[") and line.endswith("]"):
            section = line[1:-1]
        elif "=" in line:
            key, value = line.split("=", 1)
            result[f"{section}.{key.strip()}"] = value.strip()
    return result


def main() -> None:
    if sha256(REFERENCE) != REFERENCE_SHA:
        raise AssertionError("frozen WP12 R48 input hash mismatch")
    original = REFERENCE.read_bytes().decode("utf-8")
    modified = original
    modified = replace_unique(modified, "Name = WP12-R48-A1e4-6T", "Name = WP14-R48-A1e4-dt2-6T")
    modified = replace_unique(modified, "DtInitial = 113.68588647695695", f"DtInitial = {DT_TEXT}")
    modified = replace_unique(modified, "MaxTimeStepSize = 113.68588647695695", f"MaxTimeStepSize = {DT_TEXT}")
    modified = replace_unique(modified, "OutputInterval = 50", "OutputInterval = 100")
    before, after = parse_ini(original), parse_ini(modified)
    changed = {key: {"reference": before[key], "new": after[key]} for key in before if before[key] != after[key]}
    allowed = {
        "Problem.Name", "TimeLoop.DtInitial", "TimeLoop.MaxTimeStepSize", "TimeLoop.OutputInterval"
    }
    if set(changed) != allowed or set(before) != set(after):
        raise AssertionError({"changed": changed, "missing_or_added": set(before) ^ set(after)})
    if float(after["TimeLoop.TEnd"]) != float(before["TimeLoop.TEnd"]):
        raise AssertionError("TEnd changed")
    OUT.mkdir(parents=True, exist_ok=True)
    params = OUT / "params.input"
    params.write_bytes(modified.encode("utf-8"))
    input_sha = sha256(params)

    plan = {
        "schema_version": 1,
        "campaign": "WP14 controlled R48 timestep refinement",
        "grid": "R48",
        "run_name": "WP14-R48-A1e4-dt2-6T",
        "threads": 4,
        "reference_commit": "213771645aa2cb45dad7e4f46e24c18637679ed3",
        "reference_input_sha256": REFERENCE_SHA,
        "input_sha256": input_sha,
        "profile_sha256": PROFILE_SHA,
        "executable_sha256": ELF_SHA,
        "nominal_dt_s": DT,
        "nominal_dt_text": DT_TEXT,
        "output_interval_steps": 100,
        "expected_frames_if_horizon": 145,
        "authorized_steps": STEPS,
        "authorized_horizon_T": 6.0,
        "authorized_end_s": T_END,
        "authorized_amplitude": 1e-4,
        "solver_timeout_s": 14400,
        "minimum_scratch_free_bytes": 32 * 1024**3,
        "minimum_available_memory_bytes": 16 * 1024**3,
        "estimated_wall_s_from_wp12_dt_scaling": 4686.640,
        "requested_wall_time": "05:00:00",
        "estimated_durable_output_bytes": 85 * 1024**2,
        "allowed_input_differences": changed,
        "initial_salt_mole_fraction_range": [0.0035, 0.0037],
        "production_simulation_limit": 1,
        "automatic_retry_authorized": False,
    }
    (OUT / "plan.json").write_text(json.dumps(plan, indent=2))

    job = f"""#!/bin/bash
#SBATCH --job-name=WP14-R48-dt2
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=05:00:00
#SBATCH --output=production.stdout
#SBATCH --error=production.stderr
set -euo pipefail
export PATH=/rmprog/slurm/v22.05.7/bin:/usr/local/bin:/usr/bin:/bin
case "$(hostname)" in ibcln*|test0*|gln0*) echo "refusing login/test node"; exit 2;; esac
root="$1"; export WP11_THREADS=4
srun -n 1 -c 4 /WORK/hulm_work/ys/miniconda3/bin/python3 "$root/wp11_hpcc_worker.py" "$root" R48
"""
    (OUT / "job.sbatch").write_bytes(job.encode("utf-8"))

    preflight_job = """#!/bin/bash
#SBATCH --job-name=WP14-dt2-preflight
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:05:00
#SBATCH --output=preflight.stdout
#SBATCH --error=preflight.stderr
set -euo pipefail
export PATH=/rmprog/slurm/v22.05.7/bin:/usr/local/bin:/usr/bin:/bin
case "$(hostname)" in ibcln*|test0*|gln0*) echo "refusing login/test node"; exit 2;; esac
root="$1"
srun -n 1 -c 1 /WORK/hulm_work/ys/miniconda3/bin/python3 "$root/wp14_compute_preflight.py" "$root"
"""
    (OUT / "preflight.sbatch").write_bytes(preflight_job.encode("utf-8"))

    audit = {
        "reference_input": str(REFERENCE.relative_to(PAPER2)),
        "reference_sha256": REFERENCE_SHA,
        "new_input_sha256": input_sha,
        "changed_keys": changed,
        "unchanged_key_count": len(before) - len(changed),
        "all_other_parameters_byte_preserved": True,
        "output_interval_rationale": "100 dt/2 steps equals the original 50 dt steps at matching physical output times",
    }
    (OUT / "input_diff_audit.json").write_text(json.dumps(audit, indent=2))
    print(json.dumps({"input_sha256": input_sha, "changed_keys": sorted(changed), "output": str(OUT)}, indent=2))


if __name__ == "__main__":
    main()
