#!/bin/bash
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=02:30:00
#SBATCH --job-name=p2h-crystal-grid
#SBATCH --output=matrix-%j.out
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || exit 90
case "$(hostname)" in ibcln*|login*|test*) exit 91;; esac
module purge
module load compilers/gcc/v12.2.0 mpi/openmpi-gcc12/v4.1.4
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /WORK/hulm_work/zhixin/p2_h_crystal_boundary_20261108
sha256sum -c smoke/binary.sha256
hostname > matrix_node.txt
echo "$SLURM_JOB_ID" > matrix_job.txt
python3 repo/model/p2_h_ale_coupling/verify_crystal.py smoke
python3 repo/model/p2_h_ale_coupling/run_crystal.py matrix
