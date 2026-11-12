#!/bin/bash
#SBATCH --job-name=p2h27-gate
#SBATCH --partition=cnmix
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:05:00
#SBATCH --output=logs/pilotgate.%j.out
#SBATCH --error=logs/pilotgate.%j.err
set -euo pipefail
test -n "${SLURM_JOB_ID:-}"
case "$(hostname)" in ibcln*|test*|ibtest*) exit 90;; esac
cd /WORK/hulm_work/zhixin/p2_h_double_moving_boundary_20261112
python3 pilot_gate.py
tar -czf pilot_evidence.tar.gz raw/pilot* logs verification*_sha256.txt production_sha256.txt pilot_gate.json verified_gate.txt
sha256sum pilot_evidence.tar.gz > pilot_evidence.sha256
