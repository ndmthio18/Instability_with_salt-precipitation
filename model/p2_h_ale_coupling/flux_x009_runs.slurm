#!/bin/bash
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=04:00:00
#SBATCH --array=0-2%3
#SBATCH --job-name=p2h-x009-runs
#SBATCH --output=logs/run-%A_%a.out
#SBATCH --error=logs/run-%A_%a.err
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || exit 90
case "$(hostname)" in ibcln*|login*|test*) exit 91;; esac
module purge
module load compilers/gcc/v12.2.0 mpi/openmpi-gcc12/v4.1.4 tools/cmake/v3.25.2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 DUMUX_NUM_THREADS=1
ROOT=/WORK/hulm_work/zhixin/p2_h_flux_mechanism_x009_20261109
cd "$ROOT"
test -x build/soil_solver
python3 run_flux_x009.py "$SLURM_ARRAY_TASK_ID"
