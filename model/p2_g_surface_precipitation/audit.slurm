#!/bin/bash
#SBATCH --partition=cnmix
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:05:00
#SBATCH --job-name=p2g-audit
#SBATCH --output=audit-%j.out
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || exit 90
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
cd /WORK/hulm_work/zhixin/p2_g_surface_precipitation_20261107
python3 audit.py
sha256sum audit.py > results/audit_source.sha256
echo "$SLURM_JOB_ID" > results/audit_job.txt
tar -czf accepted_histories.tar.gz results/replay
tar --exclude=results/replay -czf compact_results.tar.gz results
sha256sum accepted_histories.tar.gz compact_results.tar.gz > transfer.sha256
