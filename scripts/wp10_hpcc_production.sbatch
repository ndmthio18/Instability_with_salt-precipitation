#!/bin/bash
#SBATCH --job-name=WP10-R24-long
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=05:00:00
#SBATCH --output=production.stdout.%j
#SBATCH --error=production.stderr.%j
set -euo pipefail
export PATH=/usr/local/bin:/usr/bin:/bin

case "$(hostname)" in ibcln*|test0*|gln0*) echo "refusing login/test node"; exit 2;; esac
root="$1"
export WP10_THREADS="${SLURM_CPUS_PER_TASK}"
python3 "$root/wp10_hpcc_worker.py" "$root"

