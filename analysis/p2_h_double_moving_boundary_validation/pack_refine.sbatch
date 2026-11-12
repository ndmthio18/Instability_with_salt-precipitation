#!/bin/bash
#SBATCH --job-name=p2h27-refpack
#SBATCH --partition=cnmix
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:05:00
#SBATCH --output=logs/refpack.%j.out
#SBATCH --error=logs/refpack.%j.err
set -euo pipefail
test -n "${SLURM_JOB_ID:-}"
cd /WORK/hulm_work/zhixin/p2_h_double_moving_boundary_20261112
sha256sum raw/*n512* > provenance/refine_raw_sha256.txt
sacct -S 2026-11-12 -u hulm --name=p2h27-verify,p2h27-pilot,p2h27-gate,p2h27-main,p2h27-pack,p2h27-refine,p2h27-refpack -X -P -o JobID,JobName,State,ExitCode,Start,End,Elapsed,AllocCPUS,NodeList > provenance/slurm_jobs_final.csv
tar -czf task27_refine.tar.gz raw/*n512* logs/refine.* logs/refpack.* provenance/refine_raw_sha256.txt provenance/slurm_jobs_final.csv cases_refine.csv run_refine.sbatch pack_refine.sbatch
sha256sum task27_refine.tar.gz > task27_refine.sha256
