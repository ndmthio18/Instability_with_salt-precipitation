#!/bin/bash
#SBATCH --job-name=p2h28-final
#SBATCH --partition=cnmix
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=0-14%15
#SBATCH --time=02:00:00
#SBATCH --output=logs/final.%A_%a.out
#SBATCH --error=logs/final.%A_%a.err
set -euo pipefail
test -n "${SLURM_JOB_ID:-}"
case "$(hostname)" in ibcln*|test*|ibtest*) exit 90;; esac
cd /WORK/hulm_work/zhixin/p2_h_double_front_audit_20261112
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python3 -c 'import json; a=json.load(open("fast_p01_gate.json")); b=json.load(open("p2_kinetic_gate.json")); assert a["P0_pass"] and a["P1_pass"] and b["passed"]'
line=$((SLURM_ARRAY_TASK_ID+2))
row=$(sed -n "${line}p" final_cases.csv)
IFS=, read -r id nc nl kappa gap rtol maxstep u end <<< "$row"
export FRESH_JAC=setup
./run_case_sqrt_v7 "$id" "$nc" "$nl" "$kappa" "$gap" "$rtol" "$maxstep" "$u" "$end" forcing.csv
