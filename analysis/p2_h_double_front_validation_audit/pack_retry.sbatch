#!/bin/bash
#SBATCH --job-name=p2h28-gate-pack
#SBATCH --partition=cnmix
#SBATCH --ntasks=1
#SBATCH --time=00:10:00
#SBATCH --output=logs/retrypack.%j.out
#SBATCH --error=logs/retrypack.%j.err
set -euo pipefail
cd /WORK/hulm_work/zhixin/p2_h_double_front_audit_20261112
test -n "${SLURM_JOB_ID:-}"
python3 gate_p01.py
mkdir -p provenance
sha256sum raw/setup* raw/*fresh* raw/*nonlinear* > provenance/retry_raw_sha256.txt
sacct -S 2026-11-12 -u hulm --name=p2h28-setup-matrix,p2h28-setup,p2h28-fresh-jac,p2h28-nonlinear -X -P -o JobID,JobName,State,ExitCode,Start,End,Elapsed,AllocCPUS,NodeList > provenance/retry_slurm.csv
tar --exclude='*_profiles.csv' -czf retry_light.tar.gz raw/setup* raw/*fresh* raw/*nonlinear* logs source_setup_v6 source_nonlinear_v5 provenance *sha256.txt p01_gate.json
sha256sum retry_light.tar.gz > retry_light.sha256
