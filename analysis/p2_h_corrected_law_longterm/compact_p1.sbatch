#!/bin/bash
#SBATCH --job-name=p2h26-p1pack
#SBATCH --partition=cnmix
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:05:00
#SBATCH --output=logs/p1pack.%j.out
#SBATCH --error=logs/p1pack.%j.err
set -euo pipefail
test -n "${SLURM_JOB_ID:-}"
case "$(hostname)" in ibcln*|test*|ibtest*) exit 90;; esac
cd /WORK/hulm_work/zhixin/p2_h_corrected_law_longterm_20261112
sacct -j 16571663,16571666,16571705,16571730,16571732 --format=JobID,State,ExitCode,Elapsed,NodeList -P > accounting.txt
python3 compact.py p1
sha256sum compact_p1.tar.gz > compact_p1.sha256
