#!/bin/bash
#SBATCH --job-name=p2h28-full-archive
#SBATCH --partition=cnmix
#SBATCH --ntasks=1
#SBATCH --time=00:15:00
#SBATCH --output=logs/fullarchive.%j.out
#SBATCH --error=logs/fullarchive.%j.err
set -euo pipefail
cd /WORK/hulm_work/zhixin/p2_h_double_front_audit_20261112
test -n "${SLURM_JOB_ID:-}"
find raw -type f -print0 | sort -z | xargs -0 sha256sum > provenance/final_raw_sha256.txt
sha256sum stage_evidence.tar.gz > stage_evidence.sha256
tar -czf task28_complete_raw.tar.gz raw source_legacy source_smoke_v0 source_double_v1 source_longdouble_v2 source_early_v3 source_nonlinear_v5 source_setup_v6 source_sqrt_v7 source_parallel_v8 logs provenance *.sbatch *_cases.csv *.py *.cc *.hh *gate* *sha256.txt forcing.csv
sha256sum task28_complete_raw.tar.gz > task28_complete_raw.sha256
