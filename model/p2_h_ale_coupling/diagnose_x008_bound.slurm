#!/bin/bash
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=01:00:00
#SBATCH --job-name=p2h-x008-bound
#SBATCH --output=logs/bound-%j.out
#SBATCH --error=logs/bound-%j.err
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || exit 90
case "$(hostname)" in ibcln*|login*|test*) exit 91;; esac
module purge
module load compilers/gcc/v12.2.0 mpi/openmpi-gcc12/v4.1.4 tools/cmake/v3.25.2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 DUMUX_NUM_THREADS=1
export P2H_SALT_MOLE_FRACTION=0.08
ROOT=/WORK/hulm_work/zhixin/p2_h_kr_av_reservoir_x008_20261108
cd "$ROOT"
python3 repo/model/p2_h_ale_coupling/diagnose_x008_bound.py
