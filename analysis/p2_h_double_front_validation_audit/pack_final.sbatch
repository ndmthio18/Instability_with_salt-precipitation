#!/bin/bash
#SBATCH --job-name=p2h28-final-pack
#SBATCH --partition=cnmix
#SBATCH --ntasks=1
#SBATCH --time=00:10:00
#SBATCH --output=logs/finalpack.%j.out
#SBATCH --error=logs/finalpack.%j.err
set -euo pipefail
cd /WORK/hulm_work/zhixin/p2_h_double_front_audit_20261112
test -n "${SLURM_JOB_ID:-}"
mkdir -p provenance
hostname > provenance/pack_node.txt
sacct -S 2026-11-12 -u hulm --name=p2h28-build,p2h28-p0-full,p2h28-early,p2h28-early-build,p2h28-diagnostic,p2h28-setup-matrix,p2h28-setup,p2h28-fresh-jac,p2h28-nonlinear,p2h28-p2-grid,p2h28-production,p2h28-final,p2h28-fast-p01,p2h28-stress-repro,p2h28-sqrt -X -P -o JobID,JobName,State,ExitCode,Start,End,Elapsed,AllocCPUS,NodeList > provenance/slurm_jobs.csv
find raw -type f -print0 | sort -z | xargs -0 sha256sum > provenance/raw_sha256.txt
sha256sum source*/* run_case* verify forcing.csv > provenance/sources_sha256.txt
python3 pack_final.py
