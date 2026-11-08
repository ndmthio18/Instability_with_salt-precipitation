#!/bin/bash
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:15:00
#SBATCH --job-name=p2h-grid-visual-qa
#SBATCH --output=visual-qa-%j.out
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || exit 90
case "$(hostname)" in ibcln*|login*|test*) exit 91;; esac
module purge
module load soft/anaconda3/config
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /WORK/hulm_work/zhixin/p2_h_grid_precipitation_fix_20261108
export MPLCONFIGDIR="$PWD/build/matplotlib-cache"
hostname
echo "$SLURM_JOB_ID"
python repo/model/p2_h_ale_coupling/plot_grid.py
P2H_DELIVERY_TAG=v2 python repo/model/p2_h_ale_coupling/collect_grid.py
