#!/bin/bash
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --job-name=p2f_analyze
#SBATCH --output=logs/analyze-%j.out
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || { echo "Slurm allocation required" >&2; exit 90; }
case "$(hostname)" in ibcln*|login*|test*) echo "Refusing login-node computation" >&2; exit 91;; esac
ROOT=/WORK/hulm_work/zhixin/p2_f_analytical_growth_20261107
python3 "$ROOT/source/model/p2_f_analytical_growth/verification/analyze_coupled.py" \
  "$ROOT/results/${P2F_RESULT_SET:-coupled}"
echo P2F_ANALYSIS_COMPLETE
