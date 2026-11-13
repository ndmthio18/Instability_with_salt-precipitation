#!/bin/bash
#SBATCH --job-name=p2h28-parallel
#SBATCH --partition=cnmix
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:15:00
#SBATCH --output=logs/parallel.%j.out
#SBATCH --error=logs/parallel.%j.err
set -euo pipefail
test -n "${SLURM_JOB_ID:-}"
case "$(hostname)" in ibcln*|test*|ibtest*) exit 90;; esac
cd /WORK/hulm_work/zhixin/p2_h_double_front_audit_20261112
module purge
module load compilers/gcc/v12.2.0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
deps=/WORK/hulm_work/zhixin/p2_h_ale_coupling_20261107/dependencies/install
g++ -O2 -fopenmp -std=c++17 -Wall -Wextra source_parallel_v8/run_case.cc -Isource_parallel_v8 -I"$deps/include" -L"$deps/lib64" -lsundials_cvode -lsundials_nvecserial -lsundials_sunmatrixdense -lsundials_sunlinsoldense -lsundials_generic -lm -o run_case_parallel_v8
sha256sum source_parallel_v8/* run_case_parallel_v8 forcing.csv > parallel_build_sha256.txt
g++ -O2 -fopenmp -std=c++17 -Wall -Wextra source_parallel_v8/verify.cc -Isource_parallel_v8 -I"$deps/include" -L"$deps/lib64" -lsundials_cvode -lsundials_nvecserial -lsundials_sunmatrixdense -lsundials_sunlinsoldense -lsundials_generic -lm -o verify_parallel_v8
FRESH_JAC=setup OMP_NUM_THREADS=8 DIAG_ID=verify_parallel ./verify_parallel_v8 > logs/verify_parallel.out 2> logs/verify_parallel.err
FRESH_JAC=setup OMP_NUM_THREADS=1 ./run_case_parallel_v8 parallel_one 128 32 1 1e-9 1e-9 30 -1 1000 forcing.csv
FRESH_JAC=setup OMP_NUM_THREADS=8 ./run_case_parallel_v8 parallel_eight 128 32 1 1e-9 1e-9 30 -1 1000 forcing.csv
FRESH_JAC=setup OMP_NUM_THREADS=1 ./run_case_sqrt_v7 parallel_reference 128 32 1 1e-9 1e-9 30 -1 1000 forcing.csv
cmp raw/parallel_one_history.csv raw/parallel_eight_history.csv
cmp raw/parallel_one_history.csv raw/parallel_reference_history.csv
printf 'PASS bitwise 1-thread / 8-thread / v7 history equality\n' > parallel_gate.txt
sha256sum source_parallel_v8/* run_case_parallel_v8 verify_parallel_v8 > parallel_sha256.txt
