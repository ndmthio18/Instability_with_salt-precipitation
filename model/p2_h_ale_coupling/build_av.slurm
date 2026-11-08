#!/bin/bash
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:45:00
#SBATCH --job-name=p2h-av-build
#SBATCH --output=build-%j.out
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || exit 90
case "$(hostname)" in ibcln*|login*|test*) exit 91;; esac
module purge
module load compilers/gcc/v12.2.0 mpi/openmpi-gcc12/v4.1.4 tools/cmake/v3.25.2 soft/anaconda3/config
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 DUMUX_NUM_THREADS=1
cd /WORK/hulm_work/zhixin/p2_h_av_sensitivity_20261108
test -f dependencies/install/lib64/libsundials_cvode.a
STACK="${SLURM_TMPDIR:-/tmp/p2h_av_${SLURM_JOB_ID}}/stack"
mkdir -p "$STACK" build unit mechanism stage_a/n640 smoke
ARCHIVE=/WORK/hulm_work/zhixin/Greal81_coupling_regime_20261106/staging/exact_stack_sources_r04_full.tar.gz
echo '2a29d942b87c157eb0bd58dd35577c1088932a048411e069837bf36b80cad684' "$ARCHIVE" | sha256sum -c -
tar -xzf "$ARCHIVE" -C "$STACK"
sed -i 's/\r$//' "$STACK/dune-uggrid/dune.module"
python3 repo/model/p2_h_ale_coupling/prepare_dumux.py "$STACK/kiemle2025a/appl" "$PWD/dependencies/install"
cp repo/model/p2_g_surface_precipitation/hpcc.opts "$STACK/hpcc.opts"
(cd "$STACK"; bash dune-common/bin/dunecontrol --opts=hpcc.opts all) > build/soil_build.log 2>&1
cmake --build "$STACK/kiemle2025a/build-cmake" --target paper2_transport -- -j8 >> build/soil_build.log 2>&1
cp "$STACK/kiemle2025a/build-cmake/appl/simpleRichardsnc/paper2_transport" build/soil_solver
cp "$STACK/kiemle2025a/appl/simpleRichardsnc/main.cc" build/compiled_main.cc
PREFIX="$PWD/dependencies/install"
LIBS="-lsundials_cvode -lsundials_nvecserial -lsundials_sunlinsoldense -lsundials_sunmatrixdense"
g++ -std=c++17 -O2 -Wall -Wextra repo/model/p2_h_ale_coupling/av_parameter_test.cc -o build/av_parameter_test
g++ -std=c++17 -O2 -Wall -Wextra repo/model/p2_h_ale_coupling/donor_test.cc -o build/donor_test
g++ -std=c++17 -O3 -Wall -Wextra -I"$PREFIX/include" repo/model/p2_h_ale_coupling/regression.cc -L"$PREFIX/lib64" $LIBS -o build/regression
g++ -std=c++17 -O2 -Wall -Wextra repo/model/p2_h_ale_coupling/operator_check.cc -o build/operator_check
g++ -std=c++17 -O3 -Wall -Wextra -I"$PREFIX/include" repo/model/p2_h_ale_coupling/av_mechanism.cc -L"$PREFIX/lib64" $LIBS -o build/av_mechanism
build/av_parameter_test > unit/av_parameter_test.log
build/donor_test > unit/donor_test.log
(cd stage_a/n640; ../../build/regression 640 21600 > solver.out 2>&1)
python3 repo/model/p2_h_ale_coupling/verify_reference.py > stage_a/verification.log
for AV in 0.1 1 10; do
  mkdir -p "mechanism/av${AV}"
  (cd "mechanism/av${AV}"; ../../build/av_mechanism "$AV" 64 5 > solver.out 2>&1)
done
sha256sum build/soil_solver build/av_parameter_test build/regression build/av_mechanism > unit/binaries.sha256
sha256sum repo/model/p2_h_ale_coupling/ale.hh repo/model/p2_h_ale_coupling/soil_methods.inc repo/model/p2_h_ale_coupling/grid_diagnostics.inc > unit/source.sha256
hostname > unit/node.txt
echo "$SLURM_JOB_ID" > unit/job.txt
python3 repo/model/p2_h_ale_coupling/run_av.py smoke
