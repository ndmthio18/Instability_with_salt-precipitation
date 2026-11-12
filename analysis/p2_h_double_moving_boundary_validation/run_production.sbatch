#!/bin/bash
#SBATCH --job-name=p2h27-main
#SBATCH --partition=cnmix
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=0-22%12
#SBATCH --time=01:30:00
#SBATCH --output=logs/main.%A_%a.out
#SBATCH --error=logs/main.%A_%a.err
set -euo pipefail
test -n "${SLURM_JOB_ID:-}"
case "$(hostname)" in ibcln*|test*|ibtest*) exit 90;; esac
cd /WORK/hulm_work/zhixin/p2_h_double_moving_boundary_20261112
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
grep -q '^PASS' verified_gate.txt
test -f pilot_gate.json
python3 - <<'PY'
import json
assert json.load(open('pilot_gate.json'))['pass'], 'pilot gate failed'
PY
line=$((SLURM_ARRAY_TASK_ID+2))
row=$(sed -n "${line}p" cases.csv)
IFS=, read -r id nc nl kappa gap rtol maxstep u end <<< "$row"
./run_case "$id" "$nc" "$nl" "$kappa" "$gap" "$rtol" "$maxstep" "$u" "$end" forcing.csv
