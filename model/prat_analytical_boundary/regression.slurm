#!/bin/bash
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:15:00
#SBATCH --job-name=prat-m1-regression
#SBATCH --output=logs/regression-%j.out
#SBATCH --error=logs/regression-%j.err
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || exit 90
case "$(hostname)" in ibcln*|login*|test*) exit 91;; esac
module purge
module load compilers/gcc/v12.2.0 mpi/openmpi-gcc12/v4.1.4
cd /WORK/hulm_work/zhixin/prat_growth_law_comparison_20261109
python3 repo/model/prat_analytical_boundary/verify_m1_regression.py
