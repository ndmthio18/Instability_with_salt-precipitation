#!/bin/bash
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:30:00
#SBATCH --job-name=prat-rebuild
#SBATCH --output=logs/rebuild-%j.out
#SBATCH --error=logs/rebuild-%j.err
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" && -n "${PRAT_PREVIOUS_JOB:-}" ]] || exit 90
case "$(hostname)" in ibcln*|login*|test*) exit 91;; esac
module purge
module load compilers/gcc/v12.2.0 mpi/openmpi-gcc12/v4.1.4 tools/cmake/v3.25.2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 DUMUX_NUM_THREADS=1
cd /WORK/hulm_work/zhixin/prat_analytical_boundary_20261109
STACK="/tmp/prat_${PRAT_STACK_JOB:-$PRAT_PREVIOUS_JOB}/stack"
[[ -d "$STACK/kiemle2025a/build-cmake" ]] || exit 92
mkdir -p "attempts/${PRAT_PREVIOUS_JOB}"
cp -a unit "attempts/${PRAT_PREVIOUS_JOB}/"
if [[ -d runs ]]; then mv runs "attempts/${PRAT_PREVIOUS_JOB}/"; fi
cp build/soil_solver "attempts/${PRAT_PREVIOUS_JOB}/soil_solver"
SRC=repo/model/prat_analytical_boundary
g++ -std=c++17 -O2 -Wall -Wextra -pedantic "$SRC/test.cc" -I "$SRC" -I repo/model/p2_f_analytical_growth/cpp -I repo/model/p2_h_ale_coupling -o build/unit_test
build/unit_test > unit/reference.csv 2> unit/tests.log
python3 "$SRC/check_reference.py" >> unit/tests.log
python3 "$SRC/prepare.py" "$STACK/kiemle2025a/appl"
cmake --build "$STACK/kiemle2025a/build-cmake" --target paper2_transport -- -j8 > "build/rebuild_${SLURM_JOB_ID}.log" 2>&1
cp "$STACK/kiemle2025a/build-cmake/appl/simpleRichardsnc/paper2_transport" build/soil_solver
cp "$STACK/kiemle2025a/appl/simpleRichardsnc/main.cc" build/compiled_main.cc
sha256sum build/soil_solver build/unit_test build/compiled_main.cc > unit/binary.sha256
hostname > unit/node.txt
printf '%s\n' "$SLURM_JOB_ID" > unit/job.txt
python3 "$SRC/run.py" verification
