#!/bin/bash
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:45:00
#SBATCH --job-name=p2h-C
#SBATCH --output=stage-c-%j.out
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || exit 90
case "$(hostname)" in ibcln*|login*|test*) exit 91;; esac
module purge
module load compilers/gcc/v12.2.0 mpi/openmpi-gcc12/v4.1.4 tools/cmake/v3.25.2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /WORK/hulm_work/zhixin/p2_h_ale_coupling_20261107
python3 repo/model/p2_h_ale_coupling/verify_soil.py b
mkdir -p stage_c
hostname > stage_c/node.txt
echo "$SLURM_JOB_ID" > stage_c/job.txt
sha256sum build/soil_solver > stage_c/binary.sha256
python3 repo/model/p2_h_ale_coupling/run_soil.py c
python3 repo/model/p2_h_ale_coupling/verify_soil.py c
