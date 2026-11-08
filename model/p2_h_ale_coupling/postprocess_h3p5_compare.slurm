#!/bin/bash
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --job-name=p2h-h3p5-post
#SBATCH --output=logs/post-%j.out
#SBATCH --error=logs/post-%j.err
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || exit 90
case "$(hostname)" in ibcln*|login*|test*) exit 91;; esac
module purge
module load compilers/gcc/v12.2.0 mpi/openmpi-gcc12/v4.1.4 tools/cmake/v3.25.2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 DUMUX_NUM_THREADS=1
ROOT=/WORK/hulm_work/zhixin/p2_h_reservoir_x008_h3p5_20261108
cd "$ROOT"
python3 repo/model/p2_h_ale_coupling/postprocess_h3p5_compare.py > analysis_h3p5_summary.json
