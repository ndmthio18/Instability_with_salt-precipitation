#!/bin/bash
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:30:00
#SBATCH --job-name=p2h-av-smoke
#SBATCH --output=smoke-%j.out
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || exit 90
case "$(hostname)" in ibcln*|login*|test*) exit 91;; esac
module purge
module load compilers/gcc/v12.2.0 mpi/openmpi-gcc12/v4.1.4 tools/cmake/v3.25.2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 DUMUX_NUM_THREADS=1
cd /WORK/hulm_work/zhixin/p2_h_av_sensitivity_20261108
PREFIX="$PWD/dependencies/install"
LIBS="-lsundials_cvode -lsundials_nvecserial -lsundials_sunlinsoldense -lsundials_sunmatrixdense"
g++ -std=c++17 -O2 -Wall -Wextra repo/model/p2_h_ale_coupling/av_parameter_test.cc -o build/av_parameter_test_r02
g++ -std=c++17 -O3 -Wall -Wextra -I"$PREFIX/include" repo/model/p2_h_ale_coupling/av_mechanism.cc -L"$PREFIX/lib64" $LIBS -o build/av_mechanism_r02
build/av_parameter_test_r02 > unit/av_parameter_test_r02.log
for AV in 0.1 1 10; do
  TARGET="mechanism_r02/av${AV}"
  mkdir -p "$TARGET"
  set +e
  (cd "$TARGET"; ../../build/av_mechanism_r02 "$AV" 64 5 > solver.out 2>&1)
  CODE=$?
  set -e
  echo "$CODE" > "$TARGET/exit_status.txt"
done
sha256sum build/av_parameter_test_r02 build/av_mechanism_r02 > unit/r02_binaries.sha256
python3 repo/model/p2_h_ale_coupling/run_av.py smoke_r02
