#!/bin/bash
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=18:00:00
#SBATCH --job-name=p2h-prog-pilot
#SBATCH --output=logs/pilot-%j.out
#SBATCH --error=logs/pilot-%j.err
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || exit 90
case "$(hostname)" in ibcln*|login*|test*) exit 91;; esac
module purge
module load compilers/gcc/v12.2.0 mpi/openmpi-gcc12/v4.1.4 tools/cmake/v3.25.2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 DUMUX_NUM_THREADS=1
ROOT=/WORK/hulm_work/zhixin/p2_h_progressive_activation_20261108
cd "$ROOT"
test -x build/soil_solver
python3 run_progressive.py pilot 0
