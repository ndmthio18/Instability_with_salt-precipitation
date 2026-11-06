#!/bin/bash
#SBATCH --job-name=WP12-reduce
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=02:00:00
#SBATCH --output=reduction.stdout
#SBATCH --error=reduction.stderr
set -euo pipefail
export PATH=/rmprog/slurm/v22.05.7/bin:/usr/local/bin:/usr/bin:/bin
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
case "$(hostname)" in ibcln*|test0*|gln0*) echo "refusing login/test node"; exit 2;; esac
root="$1"; fields="$2"; products="$3"; output="$4"; manifest="$5"; label="$6"; amplitude="$7"
srun -n 1 -c 2 /WORK/hulm_work/ys/miniconda3/bin/python3 "$root/wp12_reduce_archive.py" \
  "$root" "$fields" "$products" "$output" "$manifest" "$label" "$amplitude"
