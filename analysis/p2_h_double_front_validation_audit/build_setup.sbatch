#!/bin/bash
#SBATCH --job-name=p2h28-setup
#SBATCH --partition=cnmix
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:15:00
#SBATCH --output=logs/setup.%j.out
#SBATCH --error=logs/setup.%j.err
set -euo pipefail
test -n "${SLURM_JOB_ID:-}"
case "$(hostname)" in ibcln*|test*|ibtest*) exit 90;; esac
cd /WORK/hulm_work/zhixin/p2_h_double_front_audit_20261112
module purge
module load compilers/gcc/v12.2.0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
deps=/WORK/hulm_work/zhixin/p2_h_ale_coupling_20261107/dependencies/install
g++ -O2 -std=c++17 -Wall -Wextra source/run_case.cc -Isource -I"$deps/include" -L"$deps/lib64" -lsundials_cvode -lsundials_nvecserial -lsundials_sunmatrixdense -lsundials_sunlinsoldense -lsundials_generic -lm -o run_case
sha256sum source/* run_case forcing.csv > build_sha256.txt
g++ -O2 -std=c++17 -Wall -Wextra source/verify.cc -Isource -I"$deps/include" -L"$deps/lib64" -lsundials_cvode -lsundials_nvecserial -lsundials_sunmatrixdense -lsundials_sunlinsoldense -lsundials_generic -lm -o verify
set +e
FRESH_JAC=setup DIAG_ID=verify_setup ./verify > logs/verify_setup.out 2> logs/verify_setup.err
printf 'verify_exit=%s\n' "$?" > setup_status.txt
FRESH_JAC=setup EARLY_DIV=1 ./run_case setup_early_u06_b01 128 32 1 1e-10 1e-9 1 .6 100 forcing.csv > logs/setup_early.out 2> logs/setup_early.err
printf 'early_exit=%s\n' "$?" >> setup_status.txt
sha256sum source/* run_case verify > setup_sha256.txt
