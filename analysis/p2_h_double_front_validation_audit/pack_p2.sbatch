#!/bin/bash
#SBATCH --job-name=p2h28-p2-pack
#SBATCH --partition=cnmix
#SBATCH --ntasks=1
#SBATCH --time=00:10:00
#SBATCH --output=logs/p2pack.%j.out
#SBATCH --error=logs/p2pack.%j.err
set -euo pipefail
cd /WORK/hulm_work/zhixin/p2_h_double_front_audit_20261112
test -n "${SLURM_JOB_ID:-}"
python3 pack_p2.py
