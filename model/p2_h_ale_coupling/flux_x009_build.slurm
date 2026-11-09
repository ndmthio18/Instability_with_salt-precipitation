#!/bin/bash
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:45:00
#SBATCH --job-name=p2h-x009-build
#SBATCH --output=logs/build-%j.out
#SBATCH --error=logs/build-%j.err
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || exit 90
case "$(hostname)" in ibcln*|login*|test*) exit 91;; esac
module purge
module load compilers/gcc/v12.2.0 mpi/openmpi-gcc12/v4.1.4 tools/cmake/v3.25.2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 DUMUX_NUM_THREADS=1
ROOT=/WORK/hulm_work/zhixin/p2_h_flux_mechanism_x009_20261109
PREFIX=/WORK/hulm_work/zhixin/p2_h_kr_av_comparison_20261108/dependencies/install
cd "$ROOT"
mkdir -p build unit logs runs
g++ -std=c++17 -O2 -Wall -Wextra -pedantic source/progressive_test.cc -I source -o build/progressive_test
build/progressive_test > unit/progressive_test.log
g++ -std=c++17 -O2 -Wall -Wextra -pedantic source/flux_quadrature_test.cc -I source \
  -I "$PREFIX/include" "$PREFIX/lib64/libsundials_cvode.a" \
  "$PREFIX/lib64/libsundials_nvecserial.a" "$PREFIX/lib64/libsundials_sunlinsoldense.a" \
  "$PREFIX/lib64/libsundials_sunmatrixdense.a" -lm -o build/flux_quadrature_test
build/flux_quadrature_test > unit/flux_quadrature_test.log
STACK="${SLURM_TMPDIR:-/tmp/p2h_x009_${SLURM_JOB_ID}}/stack"
mkdir -p "$STACK"
ARCHIVE=/WORK/hulm_work/zhixin/Greal81_coupling_regime_20261106/staging/exact_stack_sources_r04_full.tar.gz
echo '2a29d942b87c157eb0bd58dd35577c1088932a048411e069837bf36b80cad684' "$ARCHIVE" | sha256sum -c -
tar -xzf "$ARCHIVE" -C "$STACK"
sed -i 's/\r$//' "$STACK/dune-uggrid/dune.module"
python3 repo/model/p2_h_ale_coupling/prepare_dumux.py "$STACK/kiemle2025a/appl" "$PREFIX"
cp repo/model/p2_g_surface_precipitation/hpcc.opts "$STACK/hpcc.opts"
(cd "$STACK"; bash dune-common/bin/dunecontrol --opts=hpcc.opts all) > build/soil_build.log 2>&1
cmake --build "$STACK/kiemle2025a/build-cmake" --target paper2_transport -- -j8 >> build/soil_build.log 2>&1
cp "$STACK/kiemle2025a/build-cmake/appl/simpleRichardsnc/paper2_transport" build/soil_solver
cp "$STACK/kiemle2025a/appl/simpleRichardsnc/main.cc" build/compiled_main.cc
sha256sum build/soil_solver build/progressive_test build/flux_quadrature_test source/*.hh source/*.inc \
  > unit/source_binary.sha256
hostname > unit/node.txt
printf '%s\n' "$SLURM_JOB_ID" > unit/job.txt
