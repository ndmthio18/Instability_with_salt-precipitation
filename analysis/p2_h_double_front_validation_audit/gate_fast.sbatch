#!/bin/bash
#SBATCH --job-name=p2h28-fast-gate
#SBATCH --partition=cnmix
#SBATCH --ntasks=1
#SBATCH --time=00:10:00
#SBATCH --output=logs/fastgate.%j.out
#SBATCH --error=logs/fastgate.%j.err
set -euo pipefail
cd /WORK/hulm_work/zhixin/p2_h_double_front_audit_20261112
test -n "${SLURM_JOB_ID:-}"
CASE_PREFIX=fast_ python3 gate_p01.py
