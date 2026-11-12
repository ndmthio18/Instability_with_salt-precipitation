#!/bin/bash
#SBATCH --job-name=p2h26-snapshot
#SBATCH --partition=cnmix
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:05:00
#SBATCH --output=logs/snapshot.%j.out
#SBATCH --error=logs/snapshot.%j.err
set -euo pipefail
test -n "${SLURM_JOB_ID:-}"
case "$(hostname)" in ibcln*|test*|ibtest*) exit 90;; esac
cd /WORK/hulm_work/zhixin/p2_h_corrected_law_longterm_20261112
tar --exclude='*n512*' -czf initial_complete.tar.gz raw
sha256sum initial_complete.tar.gz > initial_complete.sha256
