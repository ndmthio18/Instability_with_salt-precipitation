#!/bin/bash
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:45:00
#SBATCH --job-name=p2h-av-exact
#SBATCH --output=build-exact-%j.out
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || exit 90
case "$(hostname)" in ibcln*|login*|test*) exit 91;; esac
module purge
module load compilers/gcc/v12.2.0 mpi/openmpi-gcc12/v4.1.4 tools/cmake/v3.25.2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 DUMUX_NUM_THREADS=1
cd /WORK/hulm_work/zhixin/p2_h_av_sensitivity_20261108
test -f dependencies/install/lib64/libsundials_cvode.a
STACK="${SLURM_TMPDIR:-/tmp/p2h_av_exact_${SLURM_JOB_ID}}/stack"
mkdir -p "$STACK" build_exact unit_exact stage_a_exact/n640 mechanism_exact
ARCHIVE=/WORK/hulm_work/zhixin/Greal81_coupling_regime_20261106/staging/exact_stack_sources_r04_full.tar.gz
echo '2a29d942b87c157eb0bd58dd35577c1088932a048411e069837bf36b80cad684' "$ARCHIVE" | sha256sum -c -
tar -xzf "$ARCHIVE" -C "$STACK"
sed -i 's/\r$//' "$STACK/dune-uggrid/dune.module"
python3 repo/model/p2_h_ale_coupling/prepare_dumux.py "$STACK/kiemle2025a/appl" "$PWD/dependencies/install"
cp repo/model/p2_g_surface_precipitation/hpcc.opts "$STACK/hpcc.opts"
(cd "$STACK"; bash dune-common/bin/dunecontrol --opts=hpcc.opts all) > build_exact/soil_build.log 2>&1
cmake --build "$STACK/kiemle2025a/build-cmake" --target paper2_transport -- -j8 >> build_exact/soil_build.log 2>&1
cp "$STACK/kiemle2025a/build-cmake/appl/simpleRichardsnc/paper2_transport" build_exact/soil_solver
cp "$STACK/kiemle2025a/appl/simpleRichardsnc/main.cc" build_exact/compiled_main.cc
PREFIX="$PWD/dependencies/install"
LIBS="-lsundials_cvode -lsundials_nvecserial -lsundials_sunlinsoldense -lsundials_sunmatrixdense"
g++ -std=c++17 -O2 -Wall -Wextra repo/model/p2_h_ale_coupling/av_parameter_test.cc -o build_exact/av_parameter_test
g++ -std=c++17 -O3 -Wall -Wextra -I"$PREFIX/include" repo/model/p2_h_ale_coupling/regression.cc -L"$PREFIX/lib64" $LIBS -o build_exact/regression
g++ -std=c++17 -O3 -Wall -Wextra -I"$PREFIX/include" repo/model/p2_h_ale_coupling/av_mechanism.cc -L"$PREFIX/lib64" $LIBS -o build_exact/av_mechanism
build_exact/av_parameter_test > unit_exact/av_parameter_test.log
(cd stage_a_exact/n640; ../../build_exact/regression 640 21600 > solver.out 2>&1)
cmp stage_a_exact/n640/history.csv stage_a/n640/history.csv
cmp stage_a_exact/n640/fields.csv stage_a/n640/fields.csv
echo 'PASS exact historical M3 histories byte-identical' > stage_a_exact/verification.log
for AV in 0.1 1 10; do
  TARGET="mechanism_exact/av${AV}"
  mkdir -p "$TARGET"
  set +e
  (cd "$TARGET"; ../../build_exact/av_mechanism "$AV" 64 5 > solver.out 2>&1)
  CODE=$?
  set -e
  echo "$CODE" > "$TARGET/exit_status.txt"
done
sha256sum build_exact/soil_solver build_exact/av_parameter_test build_exact/regression build_exact/av_mechanism > unit_exact/binaries.sha256
sha256sum repo/model/p2_h_ale_coupling/ale.hh repo/model/p2_h_ale_coupling/soil_methods.inc repo/model/p2_h_ale_coupling/grid_diagnostics.inc > unit_exact/source.sha256
hostname > unit_exact/node.txt
echo "$SLURM_JOB_ID" > unit_exact/job.txt
python3 repo/model/p2_h_ale_coupling/run_av.py smoke_exact
