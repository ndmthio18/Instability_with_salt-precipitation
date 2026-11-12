#!/bin/bash
#SBATCH --job-name=p2h26-refine
#SBATCH --partition=cnmix
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=0-4%5
#SBATCH --time=00:20:00
#SBATCH --output=logs/refine.%A_%a.out
#SBATCH --error=logs/refine.%A_%a.err
set -euo pipefail
test -n "${SLURM_JOB_ID:-}"
case "$(hostname)" in ibcln*|test*|ibtest*) exit 90;; esac
cd /WORK/hulm_work/zhixin/p2_h_corrected_law_longterm_20261112
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
line=$((SLURM_ARRAY_TASK_ID+2))
row=$(sed -n "${line}p" cases_refine.csv)
IFS=, read -r id n power rtol step u beta av jfactor <<< "$row"
./longterm "$id" "$n" "$power" "$rtol" "$step" "$u" "$beta" "$av" "$jfactor" forcing.csv
