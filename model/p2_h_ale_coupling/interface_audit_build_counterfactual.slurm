#!/bin/bash
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=00:20:00
#SBATCH --job-name=p2h-if-build
#SBATCH --output=logs/build-%j.out
#SBATCH --error=logs/build-%j.err
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || exit 90
case "$(hostname)" in ibcln*|login*|test*) exit 91;; esac
module purge
module load compilers/gcc/v12.2.0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
ROOT=/WORK/hulm_work/zhixin/p2_h_interface_concentration_audit_20261108
PREFIX=/WORK/hulm_work/zhixin/p2_h_kr_av_comparison_20261108/dependencies/install
cd "$ROOT"
mkdir -p build
g++ -std=c++17 -O2 -Wall -Wextra \
  -I "$PREFIX/include" -I source \
  source/interface_counterfactual.cc \
  "$PREFIX/lib64/libsundials_cvode.a" \
  "$PREFIX/lib64/libsundials_nvecserial.a" \
  "$PREFIX/lib64/libsundials_sunlinsoldense.a" \
  "$PREFIX/lib64/libsundials_sunmatrixdense.a" \
  -lm -o build/interface_counterfactual
sha256sum build/interface_counterfactual source/interface_counterfactual.cc source/ale.hh source/bdf.hh \
  source/pratcrust.hh > build/counterfactual_build.sha256
hostname > build/counterfactual_build.node
printf '%s\n' "$SLURM_JOB_ID" > build/counterfactual_build.job
