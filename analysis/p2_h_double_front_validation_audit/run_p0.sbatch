#!/bin/bash
#SBATCH --job-name=p2h28-p0-full
#SBATCH --partition=cnmix
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=0-5%6
#SBATCH --time=02:00:00
#SBATCH --output=logs/p0full.%A_%a.out
#SBATCH --error=logs/p0full.%A_%a.err
set -euo pipefail
test -n "${SLURM_JOB_ID:-}"
case "$(hostname)" in ibcln*|test*|ibtest*) exit 90;; esac
cd /WORK/hulm_work/zhixin/p2_h_double_front_audit_20261112
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
line=$((SLURM_ARRAY_TASK_ID+2))
row=$(sed -n "${line}p" p0_cases.csv)
IFS=, read -r id nc nl kappa gap rtol maxstep u end weight jac scaling <<< "$row"
export WEIGHT_MODE="$weight" JAC_MODE="$jac"
if [ "$scaling" = none ]; then export SCALE_NONE=1; fi
./run_case_ld_v2 "$id" "$nc" "$nl" "$kappa" "$gap" "$rtol" "$maxstep" "$u" "$end" forcing.csv
