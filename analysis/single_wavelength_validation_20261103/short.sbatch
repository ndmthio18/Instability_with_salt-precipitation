#!/bin/bash
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --hint=nomultithread
#SBATCH --time=12:00:00
#SBATCH --array=0-3%2
#SBATCH --job-name=P2-WP5-short
#SBATCH --output=/tmp/zhixin_paper2_wp5_20261103_r01/short-%A_%a.log
set -eu
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
runs=(R24-M0-short R24-M1-short R48-M0-short R48-M1-short)
srun --mpi=none -n 1 -c 1 --cpu-bind=cores python3 /tmp/zhixin_paper2_wp5_20261103_r01/single_wavelength_worker.py /tmp/zhixin_paper2_wp5_20261103_r01 "${runs[$SLURM_ARRAY_TASK_ID]}"
