#!/bin/bash
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=04:00:00
#SBATCH --job-name=p2h-x010-sat
#SBATCH --output=logs/run-%j.out
#SBATCH --error=logs/run-%j.err

set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || exit 90
case "$(hostname)" in ibcln*|login*|test*) exit 91;; esac

module purge
module load compilers/gcc/v12.2.0 mpi/openmpi-gcc12/v4.1.4 tools/cmake/v3.25.2
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export DUMUX_NUM_THREADS=1

ROOT=/WORK/hulm_work/zhixin/p2_h_x010_saturated_20261110
cd "$ROOT"
test -x /WORK/hulm_work/zhixin/p2_h_flux_mechanism_x009_20261109/build/soil_solver

# Minimal runtime acceptance before the requested production calculation.
python3 run_x010_saturated.py --name smoke_x010_2s --horizon 2
python3 run_x010_saturated.py --name x010_saturated_marc_n80_12500s --horizon 12500
