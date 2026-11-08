#!/bin/bash
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=01:30:00
#SBATCH --job-name=p2h-donor-complete
#SBATCH --output=resume-%j.out
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || exit 90
case "$(hostname)" in ibcln*|login*|test*) exit 91;; esac
module purge
module load compilers/gcc/v12.2.0 mpi/openmpi-gcc12/v4.1.4
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /WORK/hulm_work/zhixin/p2_h_adaptive_donor_20261108
sha256sum -c matrix_binary.sha256
hostname > resume_node.txt
echo "$SLURM_JOB_ID" > resume_job.txt
python3 repo/model/p2_h_ale_coupling/resume_adaptive.py
