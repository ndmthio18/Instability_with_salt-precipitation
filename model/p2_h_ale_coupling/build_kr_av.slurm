#!/bin/bash
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:45:00
#SBATCH --job-name=p2h-krav-build
#SBATCH --output=logs/build-%j.out
#SBATCH --error=logs/build-%j.err
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || exit 90
case "$(hostname)" in ibcln*|login*|test*) exit 91;; esac
module purge
module load compilers/gcc/v12.2.0 mpi/openmpi-gcc12/v4.1.4 tools/cmake/v3.25.2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 DUMUX_NUM_THREADS=1
ROOT=/WORK/hulm_work/zhixin/p2_h_kr_av_comparison_20261108
cd "$ROOT"
mkdir -p logs build unit
test -f dependencies/install/lib64/libsundials_cvode.a
STACK="${SLURM_TMPDIR:-/tmp/p2h_krav_${SLURM_JOB_ID}}/stack"
mkdir -p "$STACK"
ARCHIVE=/WORK/hulm_work/zhixin/Greal81_coupling_regime_20261106/staging/exact_stack_sources_r04_full.tar.gz
echo '2a29d942b87c157eb0bd58dd35577c1088932a048411e069837bf36b80cad684' "$ARCHIVE" | sha256sum -c -
tar -xzf "$ARCHIVE" -C "$STACK"
sed -i 's/\r$//' "$STACK/dune-uggrid/dune.module"
python3 repo/model/p2_h_ale_coupling/prepare_dumux.py "$STACK/kiemle2025a/appl" "$ROOT/dependencies/install"
cp repo/model/p2_g_surface_precipitation/hpcc.opts "$STACK/hpcc.opts"
(cd "$STACK"; bash dune-common/bin/dunecontrol --opts=hpcc.opts all) > build/soil_build.log 2>&1
cmake --build "$STACK/kiemle2025a/build-cmake" --target paper2_transport -- -j8 >> build/soil_build.log 2>&1
cp "$STACK/kiemle2025a/build-cmake/appl/simpleRichardsnc/paper2_transport" build/soil_solver
cp "$STACK/kiemle2025a/appl/simpleRichardsnc/main.cc" build/compiled_main.cc
g++ -std=c++17 -O2 -Wall -Wextra repo/model/p2_h_ale_coupling/av_parameter_test.cc -o build/av_parameter_test
build/av_parameter_test > unit/av_parameter_test.log
sha256sum build/soil_solver build/av_parameter_test > unit/binaries.sha256
sha256sum repo/model/p2_h_ale_coupling/ale.hh repo/model/p2_h_ale_coupling/soil_methods.inc repo/model/p2_h_ale_coupling/run_kr_av.py > unit/source.sha256
hostname > unit/node.txt
printf '%s\n' "$SLURM_JOB_ID" > unit/job.txt
