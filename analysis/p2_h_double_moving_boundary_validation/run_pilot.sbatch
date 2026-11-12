#!/bin/bash
#SBATCH --job-name=p2h27-pilot
#SBATCH --partition=cnmix
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=0-5%6
#SBATCH --time=00:20:00
#SBATCH --output=logs/pilot.%A_%a.out
#SBATCH --error=logs/pilot.%A_%a.err
set -euo pipefail
test -n "${SLURM_JOB_ID:-}"
case "$(hostname)" in ibcln*|test*|ibtest*) exit 90;; esac
cd /WORK/hulm_work/zhixin/p2_h_double_moving_boundary_20261112
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
grep -q '^PASS' verified_gate.txt
line=$((SLURM_ARRAY_TASK_ID+2))
row=$(sed -n "${line}p" pilot_cases.csv)
IFS=, read -r id nc nl kappa gap rtol maxstep u end <<< "$row"
./run_case "$id" "$nc" "$nl" "$kappa" "$gap" "$rtol" "$maxstep" "$u" "$end" forcing.csv
