#!/bin/bash
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=01:00:00
#SBATCH --array=0-7%8
#SBATCH --job-name=p2h-if-cf
#SBATCH --output=logs/counterfactual-%A_%a.out
#SBATCH --error=logs/counterfactual-%A_%a.err
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || exit 90
case "$(hostname)" in ibcln*|login*|test*) exit 91;; esac
module purge
module load compilers/gcc/v12.2.0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
ROOT=/WORK/hulm_work/zhixin/p2_h_interface_concentration_audit_20261108
PREFIX=/WORK/hulm_work/zhixin/p2_h_kr_av_comparison_20261108/dependencies/install
cd "$ROOT"
test -x build/interface_counterfactual
python3 run_counterfactuals.py "$SLURM_ARRAY_TASK_ID"
