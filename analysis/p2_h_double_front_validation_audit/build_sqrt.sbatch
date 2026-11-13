#!/bin/bash
#SBATCH --job-name=p2h28-sqrt
#SBATCH --partition=cnmix
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:15:00
#SBATCH --output=logs/sqrt.%j.out
#SBATCH --error=logs/sqrt.%j.err
set -euo pipefail
test -n "${SLURM_JOB_ID:-}"
case "$(hostname)" in ibcln*|test*|ibtest*) exit 90;; esac
cd /WORK/hulm_work/zhixin/p2_h_double_front_audit_20261112
module purge
module load compilers/gcc/v12.2.0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
deps=/WORK/hulm_work/zhixin/p2_h_ale_coupling_20261107/dependencies/install
g++ -O2 -std=c++17 -Wall -Wextra source_sqrt_v7/run_case.cc -Isource_sqrt_v7 -I"$deps/include" -L"$deps/lib64" -lsundials_cvode -lsundials_nvecserial -lsundials_sunmatrixdense -lsundials_sunlinsoldense -lsundials_generic -lm -o run_case_sqrt_v7
sha256sum source_sqrt_v7/* run_case_sqrt_v7 forcing.csv > sqrt_build_sha256.txt
g++ -O2 -std=c++17 -Wall -Wextra source_sqrt_v7/verify.cc -Isource_sqrt_v7 -I"$deps/include" -L"$deps/lib64" -lsundials_cvode -lsundials_nvecserial -lsundials_sunmatrixdense -lsundials_sunlinsoldense -lsundials_generic -lm -o verify_sqrt_v7
set +e
FRESH_JAC=setup DIAG_ID=verify_sqrt ./verify_sqrt_v7 > logs/verify_sqrt.out 2> logs/verify_sqrt.err
printf 'verify_exit=%s\n' "$?" > sqrt_status.txt
FRESH_JAC=setup EARLY_DIV=1 ./run_case_sqrt_v7 sqrt_early_u06_b01 128 32 1 1e-10 1e-9 1 .6 100 forcing.csv > logs/sqrt_early.out 2> logs/sqrt_early.err
printf 'early_exit=%s\n' "$?" >> sqrt_status.txt
sha256sum source_sqrt_v7/* run_case_sqrt_v7 verify_sqrt_v7 > sqrt_sha256.txt
