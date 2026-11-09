#!/bin/bash
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --job-name=p2h-sat-nodiss-run
#SBATCH --output=logs/run-%j.out
#SBATCH --error=logs/run-%j.err
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || exit 90
case "$(hostname)" in ibcln*|login*|test*) exit 91;; esac
module purge
module load compilers/gcc/v12.2.0 mpi/openmpi-gcc12/v4.1.4 tools/cmake/v3.25.2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 DUMUX_NUM_THREADS=1
cd /WORK/hulm_work/zhixin/p2_h_sat_interface_nodiss_x009_20261109
test -x build/soil_solver
python3 run_sat_interface_nodiss.py
sha256sum runs/saturated_no_source/* > run_outputs.sha256
tar -czf compact.tar.gz runs/saturated_no_source unit logs input.sha256 build_submission.txt run_submission.txt run_outputs.sha256
sha256sum compact.tar.gz > compact.sha256

