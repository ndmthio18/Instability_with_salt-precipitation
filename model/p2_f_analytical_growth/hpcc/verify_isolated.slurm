#!/bin/bash
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=00:15:00
#SBATCH --job-name=p2f_iso
#SBATCH --output=logs/verify-isolated-%j.out
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || { echo "Slurm allocation required" >&2; exit 90; }
case "$(hostname)" in ibcln*|login*|test*) echo "Refusing login-node computation" >&2; exit 91;; esac
module purge
module load compilers/gcc/v12.2.0
module load tools/cmake/v3.25.2
ROOT=/WORK/hulm_work/zhixin/p2_f_analytical_growth_20261107
SCRATCH="${SLURM_TMPDIR:-/tmp/p2f_iso_${SLURM_JOB_ID}}"
mkdir -p "$SCRATCH" "$ROOT/results/isolated" "$ROOT/logs"
cp -a "$ROOT/source/model/p2_f_analytical_growth" "$SCRATCH/p2f"
cd "$SCRATCH/p2f"
python3 -m unittest -v test_analytical_growth.py |& tee "$ROOT/results/isolated/python_tests.txt"
python3 run_isolated_verification.py
cmake -S cpp -B cpp/build -DCMAKE_BUILD_TYPE=Release > "$ROOT/results/isolated/cmake_configure.txt"
cmake --build cpp/build --parallel 2 > "$ROOT/results/isolated/cmake_build.txt"
cpp/build/verify_analytical | tee "$ROOT/results/isolated/cpp_verification.txt"
python3 verification/compare_cpp_python.py \
  "$ROOT/results/isolated/cpp_verification.txt" \
  "$ROOT/results/isolated/cpp_python_comparison.json"
cp -a results/. "$ROOT/results/isolated/"
{
  echo "job_id=$SLURM_JOB_ID"
  echo "host=$(hostname)"
  echo "compiler=$(g++ --version | head -1)"
  echo "cmake=$(cmake --version | head -1)"
  sha256sum analytical_growth.py cpp/analyticalgrowth.hh
} > "$ROOT/results/isolated/identity.txt"
echo P2F_ISOLATED_COMPLETE
