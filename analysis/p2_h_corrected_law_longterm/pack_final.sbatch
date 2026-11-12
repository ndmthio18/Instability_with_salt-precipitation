#!/bin/bash
#SBATCH --job-name=p2h26-finalpack
#SBATCH --partition=cnmix
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --output=logs/finalpack.%j.out
#SBATCH --error=logs/finalpack.%j.err
set -euo pipefail
test -n "${SLURM_JOB_ID:-}"
case "$(hostname)" in ibcln*|test*|ibtest*) exit 90;; esac
cd /WORK/hulm_work/zhixin/p2_h_corrected_law_longterm_20261112
sacct -j 16571663,16571666,16571705,16571730,16571732,16571762,16571763,16571786 --format=JobID,State,ExitCode,Elapsed,NodeList -P > accounting.txt
sha256sum source/* source_p0_frozen/* source_p1_frozen/* raw/* cases*.csv forcing.csv longterm > final_sha256.txt
tar -czf evidence_final.tar.gz raw logs source source_p0_frozen source_p1_frozen cases*.csv forcing.csv *.sbatch *_sha256.txt accounting.txt
sha256sum evidence_final.tar.gz > evidence_final.sha256
python3 compact.py 'p1_*n512' compact_refine
sha256sum compact_refine.tar.gz > compact_refine.sha256
