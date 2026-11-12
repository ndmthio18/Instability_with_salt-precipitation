#!/bin/bash
#SBATCH --job-name=p2h26-pack
#SBATCH --partition=cnmix
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --output=logs/pack.%j.out
#SBATCH --error=logs/pack.%j.err
set -euo pipefail
test -n "${SLURM_JOB_ID:-}"
case "$(hostname)" in ibcln*|test*|ibtest*) exit 90;; esac
cd /WORK/hulm_work/zhixin/p2_h_corrected_law_longterm_20261112
sacct -j 16571666,16571705 --format=JobID,State,ExitCode,Elapsed,NodeList -P > accounting.txt
sha256sum source/* raw/* > final_sha256.txt
tar -czf evidence.tar.gz raw logs source cases*.csv *.sbatch *_sha256.txt accounting.txt
sha256sum evidence.tar.gz > evidence.sha256
