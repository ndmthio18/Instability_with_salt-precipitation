#!/bin/bash
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --job-name=p2h-x009-post
#SBATCH --output=logs/post-%j.out
#SBATCH --error=logs/post-%j.err
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || exit 90
case "$(hostname)" in ibcln*|login*|test*) exit 91;; esac
module purge
module load soft/anaconda3
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
ROOT=/WORK/hulm_work/zhixin/p2_h_flux_mechanism_x009_20261109
cd "$ROOT"
python3 repo/model/p2_h_ale_coupling/analyze_flux_mechanism.py \
  --x008-progressive runs/x008_progressive_reservoir_n80_a64_dt0p5_diagnostic \
  --x008-adaptive /WORK/hulm_work/zhixin/p2_h_reservoir_x008_h3p5_20261108/diagnostic_bound \
  --x009 runs/x009_progressive_reservoir_n80_a64_dt0p5 \
  --x009-half runs/x009_progressive_reservoir_n80_a64_dt0p25_timecheck \
  --out analysis_output > analysis_summary.json
sha256sum analysis_output/*.csv analysis_output/verification.json analysis_output/figures/*.png \
  > analysis_output/output_sha256.txt
