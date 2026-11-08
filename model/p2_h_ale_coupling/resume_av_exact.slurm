#!/bin/bash
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=00:30:00
#SBATCH --job-name=p2h-av-exact2
#SBATCH --output=exact-resume-%j.out
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || exit 90
case "$(hostname)" in ibcln*|login*|test*) exit 91;; esac
module purge
module load compilers/gcc/v12.2.0 mpi/openmpi-gcc12/v4.1.4 tools/cmake/v3.25.2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 DUMUX_NUM_THREADS=1
cd /WORK/hulm_work/zhixin/p2_h_av_sensitivity_20261108
test -x build_exact/soil_solver
g++ -std=c++17 -O2 -Wall -Wextra repo/model/p2_h_ale_coupling/av_parameter_test.cc -o build_exact/av_parameter_test_r02
build_exact/av_parameter_test_r02 > unit_exact/av_parameter_test_r02.log
(cd stage_a_exact/n640; ../../build_exact/regression 640 21600 > solver_r02.out 2>&1)
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
sha256sum build_exact/soil_solver build_exact/av_parameter_test_r02 build_exact/regression build_exact/av_mechanism > unit_exact/binaries_r02.sha256
hostname > unit_exact/resume_node.txt
echo "$SLURM_JOB_ID" > unit_exact/resume_job.txt
python3 repo/model/p2_h_ale_coupling/run_av.py smoke_exact
