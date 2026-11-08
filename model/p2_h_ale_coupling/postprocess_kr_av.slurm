#!/bin/bash
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=01:00:00
#SBATCH --job-name=p2h-krav-post
#SBATCH --output=logs/postprocess-%j.out
#SBATCH --error=logs/postprocess-%j.err
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || exit 90
case "$(hostname)" in ibcln*|login*|test*) exit 91;; esac
module purge
module load compilers/gcc/v12.2.0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export MPLCONFIGDIR="${SLURM_TMPDIR:-/tmp}/matplotlib-${SLURM_JOB_ID}"
mkdir -p "$MPLCONFIGDIR"
cd /WORK/hulm_work/zhixin/p2_h_kr_av_comparison_20261108
python3 repo/model/p2_h_ale_coupling/analyze_kr_av.py
python3 repo/model/p2_h_ale_coupling/plot_kr_av.py
python3 repo/model/p2_h_ale_coupling/verify_kr_av.py
find analysis/results -type f | sort | xargs sha256sum > analysis/delivery_sha256.txt
sha256sum analysis/provenance.json analysis/verification.json analysis/figures/* >> analysis/delivery_sha256.txt
hostname > analysis/postprocess_node.txt
printf '%s\n' "$SLURM_JOB_ID" > analysis/postprocess_job.txt
