#!/bin/bash
#SBATCH --job-name=p2h28-grid-check
#SBATCH --partition=cnmix
#SBATCH --ntasks=1
#SBATCH --time=00:10:00
#SBATCH --output=logs/gridcheck.%j.out
#SBATCH --error=logs/gridcheck.%j.err
set -euo pipefail
cd /WORK/hulm_work/zhixin/p2_h_double_front_audit_20261112
test -n "${SLURM_JOB_ID:-}"
python3 check_grid.py
