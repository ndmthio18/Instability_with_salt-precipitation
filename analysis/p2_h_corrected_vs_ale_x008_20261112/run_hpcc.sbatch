#!/bin/bash
#SBATCH --job-name=p2h_x008_cmp
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --output=logs/stdout.%j
#SBATCH --error=logs/stderr.%j

set -euo pipefail
if [[ -z "${SLURM_JOB_ID:-}" ]]; then echo "Refusing to compute outside Slurm" >&2; exit 90; fi
case "$(hostname)" in ibcln*|ibtest*|test*) echo "Refusing login/test-node computation" >&2; exit 91;; esac

run_dir="/WORK/hulm_work/zhixin/p2_h_corrected_vs_ale_x008_20261112"
deps="/WORK/hulm_work/zhixin/p2_h_ale_coupling_20261107/dependencies/install"
cd "$run_dir"
mkdir -p build results logs
module purge
module load compilers/gcc/v12.2.0
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1

g++ -std=c++17 -O2 -Wall -Wextra source/prescribed_flux_test.cc -Isource -o build/prescribed_flux_test
./build/prescribed_flux_test | tee results/operator_test.log

g++ -std=c++17 -O2 -Wall -Wextra source/x008_prescribed_flux_ale.cc -Isource -I"$deps/include" \
  -L"$deps/lib64" -lsundials_cvode -lsundials_nvecserial -lsundials_sunmatrixdense \
  -lsundials_sunlinsoldense -lsundials_generic -lm -o build/x008_prescribed_flux_ale

./build/x008_prescribed_flux_ale results/x008_forcing.csv results/ale_history.csv 64 10 \
  | tee results/ale_run.log
cp ale_status.txt results/ale_status.txt
sha256sum source/* results/x008_forcing.csv build/* > results/run_sha256.txt
