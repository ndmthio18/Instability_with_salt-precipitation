#!/bin/bash
#SBATCH --job-name=WP14-dt2-preflight
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:05:00
#SBATCH --output=preflight.stdout
#SBATCH --error=preflight.stderr
set -euo pipefail
export PATH=/rmprog/slurm/v22.05.7/bin:/usr/local/bin:/usr/bin:/bin
case "$(hostname)" in ibcln*|test0*|gln0*) echo "refusing login/test node"; exit 2;; esac
root="$1"
srun -n 1 -c 1 /WORK/hulm_work/ys/miniconda3/bin/python3 "$root/wp14_compute_preflight.py" "$root"
