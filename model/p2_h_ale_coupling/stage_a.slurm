#!/bin/bash
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:40:00
#SBATCH --job-name=p2h-A
#SBATCH --output=stage-a-%j.out
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || exit 90
case "$(hostname)" in ibcln*|login*|test*) exit 91;; esac
module purge
module load compilers/gcc/v12.2.0 tools/cmake/v3.25.2 soft/anaconda3/config
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /WORK/hulm_work/zhixin/p2_h_ale_coupling_20261107
mkdir -p dependencies build reference stage_a/n80 stage_a/n640
hostname > stage_a/node.txt
echo "$SLURM_JOB_ID" > stage_a/job.txt
echo 'd5ef5a9227d1ecb002c784def3975a6f099c6bd88fd0d6f802f124f43420e3ae  sundials-source.tar.gz' | sha256sum -c -
tar -xzf sundials-source.tar.gz -C dependencies
tar -xf m3_reference.tar -C reference
PREFIX="$PWD/dependencies/install"
cmake -S dependencies/sundials-6.7.0 -B dependencies/build -DCMAKE_INSTALL_PREFIX="$PREFIX" \
  -DCMAKE_BUILD_TYPE=Release -DBUILD_ARKODE=OFF -DBUILD_CVODES=OFF -DBUILD_IDA=OFF \
  -DBUILD_IDAS=OFF -DBUILD_KINSOL=OFF -DEXAMPLES_ENABLE_C=OFF -DEXAMPLES_ENABLE_CXX=OFF \
  -DBUILD_SHARED_LIBS=OFF -DBUILD_STATIC_LIBS=ON > build/sundials.log 2>&1
cmake --build dependencies/build -j4 >> build/sundials.log 2>&1
cmake --install dependencies/build >> build/sundials.log 2>&1
g++ -std=c++17 -O3 -Wall -Wextra -I"$PREFIX/include" source/regression.cc -L"$PREFIX/lib64" \
  -lsundials_cvode -lsundials_nvecserial -lsundials_sunlinsoldense -lsundials_sunmatrixdense \
  -o build/regression > build/regression.log 2>&1
g++ -std=c++17 -O2 -Wall -Wextra source/operator_check.cc -o build/operator_check > build/operator.log 2>&1
sha256sum build/regression build/operator_check > stage_a/binaries.sha256
(cd stage_a/n80; ../../build/regression 80 300 > solver.out 2>&1)
(cd stage_a/n640; ../../build/regression 640 21600 > solver.out 2>&1)
python source/verify_reference.py > stage_a/verification.log 2>&1
cat stage_a/verification.log
