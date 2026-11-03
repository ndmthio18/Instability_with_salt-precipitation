#!/bin/bash
#SBATCH --partition=debug
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --threads-per-core=1
set -euo pipefail
run_root="$1"
run_name="$2"
cd "${run_root}/runs/${run_name}"
srun --ntasks=1 --cpus-per-task=1 --cpu-bind=cores \
  env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python3 "${run_root}/wp7_sal_worker.py" "${run_root}" "${run_name}"
