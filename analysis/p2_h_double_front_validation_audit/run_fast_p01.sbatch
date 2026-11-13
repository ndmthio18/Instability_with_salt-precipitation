#!/bin/bash
#SBATCH --job-name=p2h28-fast-p01
#SBATCH --partition=cnmix
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=0-35%16
#SBATCH --time=02:00:00
#SBATCH --output=logs/fastp01.%A_%a.out
#SBATCH --error=logs/fastp01.%A_%a.err
set -euo pipefail
test -n "${SLURM_JOB_ID:-}"
case "$(hostname)" in ibcln*|test*|ibtest*) exit 90;; esac
cd /WORK/hulm_work/zhixin/p2_h_double_front_audit_20261112
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
line=$((SLURM_ARRAY_TASK_ID+2))
row=$(sed -n "${line}p" fast_p01_cases.csv)
IFS=, read -r id nc nl kappa gap rtol maxstep u end init div phase <<< "$row"
export INIT_MODE="$init" FRESH_JAC=setup
if [ "$phase" = early ]; then export EARLY_DIV="$div"; fi
./run_case_sqrt_v7 "$id" "$nc" "$nl" "$kappa" "$gap" "$rtol" "$maxstep" "$u" "$end" forcing.csv
