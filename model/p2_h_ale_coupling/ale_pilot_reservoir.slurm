#!/bin/bash
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=04:00:00
#SBATCH --job-name=p2h-res-ale
#SBATCH --output=logs/ale-pilot-%j.out
#SBATCH --error=logs/ale-pilot-%j.err
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || exit 90
case "$(hostname)" in ibcln*|login*|test*) exit 91;; esac
module purge
module load compilers/gcc/v12.2.0 mpi/openmpi-gcc12/v4.1.4 tools/cmake/v3.25.2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 DUMUX_NUM_THREADS=1
cd /WORK/hulm_work/zhixin/p2_h_kr_av_reservoir_20261108
python3 repo/model/p2_h_ale_coupling/run_reservoir.py ale_pilot 0 & p0=$!
python3 repo/model/p2_h_ale_coupling/run_reservoir.py ale_pilot 1 & p1=$!
wait "$p0"
wait "$p1"
python3 repo/model/p2_h_ale_coupling/verify_ale_pilot.py > ale_pilot_gate.json
