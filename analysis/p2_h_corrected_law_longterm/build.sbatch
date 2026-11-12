#!/bin/bash
#SBATCH --job-name=p2h26-build
#SBATCH --partition=cnmix
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --output=logs/build.%j.out
#SBATCH --error=logs/build.%j.err
set -euo pipefail
test -n "${SLURM_JOB_ID:-}"
case "$(hostname)" in ibcln*|test*|ibtest*) exit 90;; esac
cd /WORK/hulm_work/zhixin/p2_h_corrected_law_longterm_20261112
module purge
module load compilers/gcc/v12.2.0
deps=/WORK/hulm_work/zhixin/p2_h_ale_coupling_20261107/dependencies/install
g++ -O2 -std=c++17 -Wall -Wextra source/prescribed_flux_test.cc -Isource -o operator_test
./operator_test
g++ -O2 -std=c++17 -Wall -Wextra source/longterm.cc -Isource -I"$deps/include" -L"$deps/lib64" -lsundials_cvode -lsundials_nvecserial -lsundials_sunmatrixdense -lsundials_sunlinsoldense -lsundials_generic -lm -o longterm
sha256sum source/* longterm cases.csv forcing.csv > build_sha256.txt
