#!/bin/bash
#SBATCH --job-name=p2h28-early-build
#SBATCH --partition=cnmix
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:15:00
#SBATCH --output=logs/earlybuild.%j.out
#SBATCH --error=logs/earlybuild.%j.err
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
DIAG_ID=verify_remaining SKIP_STRESS=1 ./verify
sha256sum source/* run_case verify forcing.csv > verified_sha256.txt
printf 'PASS non-stress tests only; stress unresolved\n' > remaining_tests_gate.txt
