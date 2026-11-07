#!/bin/bash
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --job-name=p2f_build
#SBATCH --output=logs/build-%j.out
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || { echo "Slurm allocation required" >&2; exit 90; }
case "$(hostname)" in ibcln*|login*|test*) echo "Refusing login-node build" >&2; exit 91;; esac
module purge
module load compilers/gcc/v12.2.0
module load mpi/openmpi-gcc12/v4.1.4
module load tools/cmake/v3.25.2
ROOT=/WORK/hulm_work/zhixin/p2_f_analytical_growth_20261107
ARCHIVE=/WORK/hulm_work/zhixin/Greal81_coupling_regime_20261106/staging/exact_stack_sources_r04_full.tar.gz
EXPECTED=2a29d942b87c157eb0bd58dd35577c1088932a048411e069837bf36b80cad684
OUTPUT="$ROOT/build"
SCRATCH="${SLURM_TMPDIR:-/tmp/p2f_build_${SLURM_JOB_ID}}"
STACK="$SCRATCH/stack"
mkdir -p "$STACK" "$OUTPUT" "$ROOT/logs"
actual=$(sha256sum "$ARCHIVE" | awk '{print $1}')
[[ "$actual" == "$EXPECTED" ]] || { echo "stack archive SHA mismatch" >&2; exit 92; }
tar -xzf "$ARCHIVE" -C "$STACK"
sed -i 's/\r$//' "$STACK/dune-uggrid/dune.module"
cp -a "$ROOT/source/model/p1_r2_conservative_richards/source/appl/." "$STACK/kiemle2025a/appl/"
patch -d "$STACK/kiemle2025a" -p1 < "$ROOT/source/model/p2_f_analytical_growth/integration/p2f_main.patch"
cp "$ROOT/source/model/p2_f_analytical_growth/cpp/analyticalgrowth.hh" \
   "$STACK/kiemle2025a/appl/simpleRichardsnc/analyticalgrowth.hh"
cat > "$STACK/hpcc.opts" <<EOF
CMAKE_FLAGS="-DCMAKE_C_COMPILER=/apps/compilers/gcc/v12.2.0/bin/gcc -DCMAKE_CXX_COMPILER=/apps/compilers/gcc/v12.2.0/bin/g++ -DCMAKE_BUILD_TYPE=Release -DCMAKE_CXX_FLAGS_RELEASE='-O3 -DNDEBUG -march=x86-64-v3 -DNONISOTHERMAL=1 -DModelP -DREDUCEDOUTPUT=1' -DDUNE_ENABLE_PYTHONBINDINGS=OFF"
MAKE_FLAGS="-j8"
BUILDDIR=build-cmake
EOF
cd "$STACK"
bash "$STACK/dune-common/bin/dunecontrol" --opts="$STACK/hpcc.opts" all > "$OUTPUT/dunecontrol.log" 2>&1
cmake --build "$STACK/kiemle2025a/build-cmake" --target paper2_transport -- -j8 >> "$OUTPUT/dunecontrol.log" 2>&1
candidate="$STACK/kiemle2025a/build-cmake/appl/simpleRichardsnc/paper2_transport"
[[ -x "$candidate" ]] || { echo "solver not produced" >&2; exit 93; }
cp "$candidate" "$OUTPUT/paper2_transport_p2f"
chmod 755 "$OUTPUT/paper2_transport_p2f"
sha256sum "$OUTPUT/paper2_transport_p2f" > "$OUTPUT/paper2_transport_p2f.sha256"
ldd "$OUTPUT/paper2_transport_p2f" > "$OUTPUT/paper2_transport_p2f.ldd.txt"
module list 2> "$OUTPUT/modules.txt"
{
  echo "job_id=$SLURM_JOB_ID"
  echo "host=$(hostname)"
  echo "stack_archive_sha256=$actual"
  echo "compiler=$(g++ --version | head -1)"
  echo "cmake=$(cmake --version | head -1)"
} > "$OUTPUT/build_identity.txt"
echo P2F_DUMUX_BUILD_COMPLETE
