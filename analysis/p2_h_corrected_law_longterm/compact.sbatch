#!/bin/bash
#SBATCH --job-name=p2h26-compact
#SBATCH --partition=cnmix
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:05:00
#SBATCH --output=logs/compact.%j.out
#SBATCH --error=logs/compact.%j.err
set -euo pipefail
test -n "${SLURM_JOB_ID:-}"
case "$(hostname)" in ibcln*|test*|ibtest*) exit 90;; esac
cd /WORK/hulm_work/zhixin/p2_h_corrected_law_longterm_20261112
python3 compact.py
sha256sum compact.tar.gz > compact.sha256
