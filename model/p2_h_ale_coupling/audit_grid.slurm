#!/bin/bash
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:25:00
#SBATCH --job-name=p2h-grid-verify
#SBATCH --output=verification-%j.out
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || exit 90
case "$(hostname)" in ibcln*|login*|test*) exit 91;; esac
module purge
module load soft/anaconda3/config
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /WORK/hulm_work/zhixin/p2_h_grid_precipitation_fix_20261108
export MPLCONFIGDIR="$PWD/build/matplotlib-cache"
mkdir -p "$MPLCONFIGDIR"
hostname
echo "$SLURM_JOB_ID"
python --version
module list
python repo/model/p2_h_ale_coupling/analyze_grid.py diagnostic ale_diagnostic
(cd correction; python repo/model/p2_h_ale_coupling/analyze_grid.py fixed_smoke fixed_tests fixed_pre fixed_coupled)
python repo/model/p2_h_ale_coupling/analyze_ale_failure.py
python repo/model/p2_h_ale_coupling/verify_grid.py
if [[ "${1:-partial}" == final ]]; then
    python repo/model/p2_h_ale_coupling/compare_pre.py
    python repo/model/p2_h_ale_coupling/plot_grid.py
fi
