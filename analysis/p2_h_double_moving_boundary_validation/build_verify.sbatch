#!/bin/bash
#SBATCH --job-name=p2h27-verify
#SBATCH --partition=cnmix
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:15:00
#SBATCH --output=logs/verify.%j.out
#SBATCH --error=logs/verify.%j.err
set -euo pipefail
test -n "${SLURM_JOB_ID:-}"
case "$(hostname)" in ibcln*|test*|ibtest*) exit 90;; esac
cd /WORK/hulm_work/zhixin/p2_h_double_moving_boundary_20261112
printf 'NOT_VERIFIED: build and checks in progress\n' > verified_gate.txt
module purge
module load compilers/gcc/v12.2.0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
deps=/WORK/hulm_work/zhixin/p2_h_ale_coupling_20261107/dependencies/install
g++ -O2 -std=c++17 -Wall -Wextra source/verify.cc -Isource -I"$deps/include" -L"$deps/lib64" -lsundials_cvode -lsundials_nvecserial -lsundials_sunmatrixdense -lsundials_sunlinsoldense -lsundials_generic -lm -o verify
sha256sum source/* verify > verification_sha256.txt
./verify
g++ -O2 -std=c++17 -Wall -Wextra source/run_case.cc -Isource -I"$deps/include" -L"$deps/lib64" -lsundials_cvode -lsundials_nvecserial -lsundials_sunmatrixdense -lsundials_sunlinsoldense -lsundials_generic -lm -o run_case
sha256sum source/* verify run_case forcing.csv > production_sha256.txt
printf 'PASS transport, recovery, dissolution, remap, independence\n' > verified_gate.txt
