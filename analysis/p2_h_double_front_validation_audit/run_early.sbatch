#!/bin/bash
#SBATCH --job-name=p2h28-early
#SBATCH --partition=cnmix
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=0-29%16
#SBATCH --time=01:00:00
#SBATCH --output=logs/early.%A_%a.out
#SBATCH --error=logs/early.%A_%a.err
set -euo pipefail
test -n "${SLURM_JOB_ID:-}"
case "$(hostname)" in ibcln*|test*|ibtest*) exit 90;; esac
cd /WORK/hulm_work/zhixin/p2_h_double_front_audit_20261112
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
line=$((SLURM_ARRAY_TASK_ID+2))
row=$(sed -n "${line}p" early_cases.csv)
IFS=, read -r id nc nl kappa gap rtol maxstep u end init div <<< "$row"
export INIT_MODE="$init" EARLY_DIV="$div"
./run_case_early_v3 "$id" "$nc" "$nl" "$kappa" "$gap" "$rtol" "$maxstep" "$u" "$end" forcing.csv
