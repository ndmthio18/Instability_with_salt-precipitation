#!/bin/bash
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:20:00
#SBATCH --job-name=prat-analysis
#SBATCH --output=logs/analysis-%j.out
#SBATCH --error=logs/analysis-%j.err
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || exit 90
case "$(hostname)" in ibcln*|login*|test*) exit 91;; esac
module purge
module load compilers/gcc/v12.2.0 mpi/openmpi-gcc12/v4.1.4
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 DUMUX_NUM_THREADS=1
cd /WORK/hulm_work/zhixin/prat_analytical_boundary_20261109
export MPLCONFIGDIR="$PWD/.mplconfig"
mkdir -p "$MPLCONFIGDIR"
python3 repo/model/prat_analytical_boundary/analyze.py
tar -czf deliverable.tar.gz deliverable
sha256sum deliverable.tar.gz > deliverable.tar.gz.sha256
