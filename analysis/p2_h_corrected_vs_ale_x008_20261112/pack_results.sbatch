#!/bin/bash
#SBATCH --job-name=p2h_pack
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:05:00
#SBATCH --output=logs/pack_stdout.%j
#SBATCH --error=logs/pack_stderr.%j

set -euo pipefail
if [[ -z "${SLURM_JOB_ID:-}" ]]; then echo "Refusing to run outside Slurm" >&2; exit 90; fi
case "$(hostname)" in ibcln*|ibtest*|test*) echo "Refusing login/test-node work" >&2; exit 91;; esac

run_dir="/WORK/hulm_work/zhixin/p2_h_corrected_vs_ale_x008_20261112"
cd "$run_dir"
tar -czf results/final_results_16569032.tar.gz \
  results/ale_history.csv results/ale_status.txt results/ale_run.log \
  results/operator_test.log results/sacct_16569032.txt \
  results/run_sha256_16569032.txt logs/stdout.16569032 logs/stderr.16569032
sha256sum results/final_results_16569032.tar.gz > results/final_results_16569032.tar.gz.sha256
