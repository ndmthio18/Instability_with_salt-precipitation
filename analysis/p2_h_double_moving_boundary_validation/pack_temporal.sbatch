#!/bin/bash
#SBATCH --job-name=p2h27-timepack
#SBATCH --partition=cnmix
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:02:00
#SBATCH --output=logs/timepack.%j.out
#SBATCH --error=logs/timepack.%j.err
set -euo pipefail
test -n "${SLURM_JOB_ID:-}"
cd /WORK/hulm_work/zhixin/p2_h_double_moving_boundary_20261112
sha256sum raw/u0.6_strict_tighter* > provenance/temporal_raw_sha256.txt
sacct -j 16572125,16572120 -X -P -o JobID,JobName,State,ExitCode,Start,End,Elapsed,AllocCPUS,NodeList > provenance/slurm_temporal.csv
tar -czf task27_temporal.tar.gz raw/u0.6_strict_tighter* logs/refine.16572125* provenance/temporal_raw_sha256.txt provenance/slurm_temporal.csv cases_refine.csv
sha256sum task27_temporal.tar.gz > task27_temporal.sha256
