#!/bin/bash
#SBATCH --job-name=p2h26-events
#SBATCH --partition=cnmix
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=0-7%8
#SBATCH --time=00:10:00
#SBATCH --output=logs/events.%A_%a.out
#SBATCH --error=logs/events.%A_%a.err
set -euo pipefail
test -n "${SLURM_JOB_ID:-}"
case "$(hostname)" in ibcln*|test*|ibtest*) exit 90;; esac
cd /WORK/hulm_work/zhixin/p2_h_corrected_law_longterm_20261112
line=$((SLURM_ARRAY_TASK_ID+2))
row=$(sed -n "${line}p" cases_events.csv)
IFS=, read -r id n power rtol step u beta av jfactor <<< "$row"
./longterm "$id" "$n" "$power" "$rtol" "$step" "$u" "$beta" "$av" "$jfactor" forcing.csv
