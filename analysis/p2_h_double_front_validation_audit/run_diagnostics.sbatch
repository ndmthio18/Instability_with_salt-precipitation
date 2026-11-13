#!/bin/bash
#SBATCH --job-name=p2h28-diagnostic
#SBATCH --partition=cnmix
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=0-8%9
#SBATCH --time=00:15:00
#SBATCH --output=logs/diagnostic.%A_%a.out
#SBATCH --error=logs/diagnostic.%A_%a.err
set -euo pipefail
test -n "${SLURM_JOB_ID:-}"
case "$(hostname)" in ibcln*|test*|ibtest*) exit 90;; esac
cd /WORK/hulm_work/zhixin/p2_h_double_front_audit_20261112
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
line=$((SLURM_ARRAY_TASK_ID+2))
row=$(sed -n "${line}p" diagnostic_cases.csv)
IFS=, read -r id nc nl kappa gap rtol maxstep u end weight jac scaling <<< "$row"
export WEIGHT_MODE="$weight" JAC_MODE="$jac"
if [ "$scaling" = none ]; then export SCALE_NONE=1; fi
./run_case "$id" "$nc" "$nl" "$kappa" "$gap" "$rtol" "$maxstep" "$u" "$end" forcing.csv
