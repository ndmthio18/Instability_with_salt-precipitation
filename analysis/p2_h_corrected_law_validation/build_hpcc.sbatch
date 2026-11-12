#!/bin/bash
#SBATCH --job-name=p2h25-build
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=00:10:00
#SBATCH --output=logs/build_stdout.%j
#SBATCH --error=logs/build_stderr.%j

set -euo pipefail
if [[ -z "${SLURM_JOB_ID:-}" ]]; then echo "Refusing to compute outside Slurm" >&2; exit 90; fi
case "$(hostname)" in ibcln*|ibtest*|test*) echo "Refusing login/test-node work" >&2; exit 91;; esac
run_dir="/WORK/hulm_work/zhixin/p2_h_corrected_law_validation_20261112"
deps="/WORK/hulm_work/zhixin/p2_h_ale_coupling_20261107/dependencies/install"
cd "$run_dir"
mkdir -p build logs results/raw
module purge
module load compilers/gcc/v12.2.0
g++ -std=c++17 -O2 -Wall -Wextra source/prescribed_flux_test.cc -Isource \
  -o build/prescribed_flux_test
./build/prescribed_flux_test | tee results/operator_test.log
g++ -std=c++17 -O2 -Wall -Wextra source/ale_validation_driver.cc -Isource \
  -I"$deps/include" -L"$deps/lib64" -lsundials_cvode -lsundials_nvecserial \
  -lsundials_sunmatrixdense -lsundials_sunlinsoldense -lsundials_generic -lm \
  -o build/ale_validation_driver
sha256sum source/* configs/* build/* > results/build_sha256.txt
