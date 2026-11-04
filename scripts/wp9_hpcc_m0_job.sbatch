#!/bin/bash
#SBATCH --job-name=WP9-R96-M0
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=08:00:00
#SBATCH --output=stdout.%j
#SBATCH --error=stderr.%j
set -euo pipefail

if [ -z "${SLURM_JOB_ID:-}" ]; then echo "WP9: no SLURM_JOB_ID; refusing"; exit 2; fi
case "$(hostname)" in
  ibcln*|test0*|gln0*) echo "WP9: login or test node $(hostname); refusing"; exit 2;;
esac

run_root="$1"
run_name="$2"
cd "${run_root}"
export WP9_MANAGED=1
export WP9_THREADS="${SLURM_CPUS_PER_TASK:-1}"
python3 "${run_root}/wp9_hpcc_worker.py" "${run_root}" "${run_name}"
