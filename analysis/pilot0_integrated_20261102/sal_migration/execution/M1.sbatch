#!/bin/bash
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=12:00:00
#SBATCH --job-name=p2-M1
#SBATCH --output=/tmp/zhixin_paper2_pilot0_20261102_r01/M1/slurm-%j.log
set -eu
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
srun --mpi=none -n 1 -c 1 python3 /tmp/zhixin_paper2_pilot0_20261102_r01/pilot0_sal_worker.py /tmp/zhixin_paper2_pilot0_20261102_r01 M1
