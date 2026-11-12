#!/bin/bash
#SBATCH --job-name=p2h25-pack
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --output=logs/pack_stdout.%j
#SBATCH --error=logs/pack_stderr.%j

set -euo pipefail
if [[ -z "${SLURM_JOB_ID:-}" ]]; then echo "Refusing to run outside Slurm" >&2; exit 90; fi
case "$(hostname)" in ibcln*|ibtest*|test*) echo "Refusing login/test-node work" >&2; exit 91;; esac
run_dir="/WORK/hulm_work/zhixin/p2_h_corrected_law_validation_20261112"
cd "$run_dir"
sacct -j 16570915 --format=JobID,JobName,State,ExitCode,Elapsed,NodeList -P > results/sacct_16570915.txt
sha256sum source/* configs/* build/* results/raw/* > results/final_sha256.txt
tar -czf results/p2_h_corrected_law_validation_16570915.tar.gz source configs *.sbatch staged_sha256.txt retry_16570777_fix_sha256.txt results/raw results/operator_test.log results/build_sha256.txt results/sacct_16570915.txt results/final_sha256.txt failed_attempts
sha256sum results/p2_h_corrected_law_validation_16570915.tar.gz > results/p2_h_corrected_law_validation_16570915.tar.gz.sha256
