#!/bin/bash
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:05:00
#SBATCH --job-name=p2h-figure-qa
#SBATCH --output=figure-qa-%j.out
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || exit 90
case "$(hostname)" in ibcln*|login*|test*) exit 91;; esac
module purge
module load soft/anaconda3/config
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /WORK/hulm_work/zhixin/p2_h_adaptive_donor_20261108
export MPLCONFIGDIR="$PWD/build/matplotlib-cache"
hostname
python repo/model/p2_h_ale_coupling/plot_adaptive.py
tar -czf figure_revision.tar.gz audit/figures audit/column_metadata.json repo/model/p2_h_ale_coupling/plot_adaptive.py repo/model/p2_h_ale_coupling/plot_revision.slurm
sha256sum figure_revision.tar.gz > figure_revision.sha256
