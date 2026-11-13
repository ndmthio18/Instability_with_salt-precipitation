#!/bin/bash
#SBATCH --job-name=p2h28-reproduce
#SBATCH --partition=cnmix
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --output=logs/reproduce.%j.out
#SBATCH --error=logs/reproduce.%j.err
set -euo pipefail
test -n "${SLURM_JOB_ID:-}"
case "$(hostname)" in ibcln*|test*|ibtest*) exit 90;; esac
cd /WORK/hulm_work/zhixin/p2_h_double_front_audit_20261112
module purge
module load compilers/gcc/v12.2.0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
deps=/WORK/hulm_work/zhixin/p2_h_ale_coupling_20261107/dependencies/install
g++ -O2 -std=c++17 -Wall -Wextra source_legacy/run_case.cc -Isource_legacy -I"$deps/include" -L"$deps/lib64" -lsundials_cvode -lsundials_nvecserial -lsundials_sunmatrixdense -lsundials_sunlinsoldense -lsundials_generic -lm -o legacy
sha256sum source_legacy/* legacy forcing.csv > legacy_sha256.txt
set +e
./legacy legacy_r9 128 32 1 1e-9 1e-9 15 .6 100 forcing.csv
./legacy legacy_r10 128 32 1 1e-9 1e-10 15 .6 100 forcing.csv
./legacy legacy_r8 128 32 1 1e-9 1e-8 30 .6 100 forcing.csv
set -e
tar -czf legacy_evidence.tar.gz raw/legacy* logs legacy_sha256.txt
sha256sum legacy_evidence.tar.gz > legacy_evidence.sha256
