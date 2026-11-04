#!/bin/bash
#SBATCH --job-name=WP9-portability
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:20:00
#SBATCH --output=stdout.%j
#SBATCH --error=stderr.%j
set -euo pipefail

# Refuse to compute outside a Slurm allocation or on a login/test node.
if [ -z "${SLURM_JOB_ID:-}" ]; then echo "WP9: no SLURM_JOB_ID; refusing"; exit 2; fi
case "$(hostname)" in
  ibcln*|test0*|gln0*) echo "WP9: login or test node $(hostname); refusing"; exit 2;;
esac

run_root="${SLURM_SUBMIT_DIR}"
cd "${run_root}"

{
  echo "hostname=$(hostname)"
  echo "job=${SLURM_JOB_ID}"
  echo "partition=${SLURM_JOB_PARTITION:-}"
  echo "cpus=${SLURM_CPUS_PER_TASK:-}"
  echo "kernel=$(uname -r)"
  echo "host_glibc=$(ldd --version 2>/dev/null | head -1)"
  echo "python=$(command -v python3) $(python3 --version 2>&1)"
  echo "elf_sha256=$(sha256sum bundle/bin/paper2_transport | cut -d' ' -f1)"
  echo "input_sha256=$(sha256sum params.input | cut -d' ' -f1)"
  echo "profile_sha256=$(sha256sum eigenprofile.dat | cut -d' ' -f1)"
  echo "loader_sha256=$(sha256sum bundle/lib/ld-linux-x86-64.so.2 | cut -d' ' -f1)"
  echo "libc_sha256=$(sha256sum bundle/lib/libc.so.6 | cut -d' ' -f1)"
} > environment.txt
module list > modules.txt 2>&1 || true

start=$(date +%s.%N)
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python3 bundle/launch_solver.py params.input > solver.log 2>&1
rc=$?
end=$(date +%s.%N)
echo "returncode=${rc}" >> environment.txt
echo "solver_wall_s=$(echo "${end} - ${start}" | bc)" >> environment.txt
echo "vtu_files=$(ls -1 *.vtu 2>/dev/null | wc -l)" >> environment.txt
exit ${rc}
