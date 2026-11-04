#!/bin/bash
#SBATCH --job-name=WP10-R24-bench
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:15:00
#SBATCH --output=benchmark.stdout.%j
#SBATCH --error=benchmark.stderr.%j
set -euo pipefail

case "$(hostname)" in ibcln*|test0*|gln0*) echo "refusing login/test node"; exit 2;; esac
root="$1"
scratch="${TMPDIR:-/tmp}/wp10_benchmark_${SLURM_JOB_ID}"
mkdir -p "$scratch"
trap 'rm -rf "$scratch"' EXIT
cp -aL "$root/bundle" "$scratch/bundle"
cp "$root/params_benchmark.input" "$scratch/params.input"
cp "$root/eigenprofile.dat" "$scratch/eigenprofile.dat"

for threads in 1 4; do
  run="$scratch/t${threads}"
  mkdir "$run"
  cp "$scratch/params.input" "$scratch/eigenprofile.dat" "$run/"
  start=$(date +%s.%N)
  (cd "$run" && env DUMUX_NUM_THREADS="$threads" OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 \
    python3 "$scratch/bundle/launch_solver.py" params.input > solver.log 2>&1)
  end=$(date +%s.%N)
  printf '%s\n' "$start" "$end" > "$run/times.txt"
done

python3 "$root/wp10_benchmark_reduce.py" "$scratch" "$root" "${SLURM_JOB_ID}"

