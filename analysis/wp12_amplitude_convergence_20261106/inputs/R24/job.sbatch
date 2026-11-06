#!/bin/bash
#SBATCH --job-name=WP12-R24-A1e4
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=03:00:00
#SBATCH --output=production.stdout
#SBATCH --error=production.stderr
set -euo pipefail
export PATH=/rmprog/slurm/v22.05.7/bin:/usr/local/bin:/usr/bin:/bin
case "$(hostname)" in ibcln*|test0*|gln0*) echo "refusing login/test node"; exit 2;; esac
root="$1"; export WP11_THREADS=4
srun -n 1 -c 4 /WORK/hulm_work/ys/miniconda3/bin/python3 "$root/wp11_hpcc_worker.py" "$root" R24
