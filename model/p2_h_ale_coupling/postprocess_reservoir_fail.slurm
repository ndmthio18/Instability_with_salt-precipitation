#!/bin/bash
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --job-name=p2h-res-failpost
#SBATCH --output=logs/failpost-%j.out
#SBATCH --error=logs/failpost-%j.err
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || exit 90
case "$(hostname)" in ibcln*|login*|test*) exit 91;; esac
module purge
module load compilers/gcc/v12.2.0 mpi/openmpi-gcc12/v4.1.4 tools/cmake/v3.25.2
cd /WORK/hulm_work/zhixin/p2_h_kr_av_reservoir_20261108
python3 repo/model/p2_h_ale_coupling/postprocess_reservoir_fail.py > analysis_fail_summary.json
