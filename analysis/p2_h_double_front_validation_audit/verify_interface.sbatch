#!/bin/bash
#SBATCH --job-name=p2h28-interface-derivative
#SBATCH --partition=cnmix
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:15:00
#SBATCH --output=logs/interface.%j.out
#SBATCH --error=logs/interface.%j.err
set -euo pipefail
test -n "${SLURM_JOB_ID:-}"
case "$(hostname)" in ibcln*|test*|ibtest*) exit 90;; esac
cd /WORK/hulm_work/zhixin/p2_h_double_front_audit_20261112
module purge
module load compilers/gcc/v12.2.0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
deps=/WORK/hulm_work/zhixin/p2_h_ale_coupling_20261107/dependencies/install
g++ -O2 -std=c++17 -Wall -Wextra verify_interface.cc -Isource_sqrt_v7 -o verify_interface
./verify_interface
sha256sum verify_interface.cc verify_interface source_sqrt_v7/* > interface_test_sha256.txt
