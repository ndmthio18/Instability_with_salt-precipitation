#!/bin/bash
#SBATCH --job-name=p2h27-pack
#SBATCH --partition=cnmix
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --output=logs/pack.%j.out
#SBATCH --error=logs/pack.%j.err
set -euo pipefail
test -n "${SLURM_JOB_ID:-}"
cd /WORK/hulm_work/zhixin/p2_h_double_moving_boundary_20261112
mkdir -p provenance
date -u > provenance/pack_time.txt
hostname > provenance/pack_node.txt
sha256sum source/* verify run_case forcing.csv > provenance/final_source_sha256.txt
find raw -type f -print0 | sort -z | xargs -0 sha256sum > provenance/raw_sha256.txt
sacct -S 2026-11-12 -u hulm --name=p2h27-verify,p2h27-pilot,p2h27-gate,p2h27-main,p2h27-pack -X -P -o JobID,JobName,State,ExitCode,Start,End,Elapsed,AllocCPUS,NodeList > provenance/slurm_jobs.csv
tar -czf task27_evidence.tar.gz raw logs source source_* provenance *_sha256.txt *.sbatch *_cases.csv cases.csv forcing.csv pilot_gate.json pilot_q_gate_failed.json verified_gate.txt
sha256sum task27_evidence.tar.gz > task27_evidence.sha256
