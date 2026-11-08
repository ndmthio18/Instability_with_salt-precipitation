#!/bin/bash
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:20:00
#SBATCH --job-name=p2h-donor-audit
#SBATCH --output=audit-%j.out
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || exit 90
case "$(hostname)" in ibcln*|login*|test*) exit 91;; esac
module purge
module load soft/anaconda3/config
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /WORK/hulm_work/zhixin/p2_h_adaptive_donor_20261108
export MPLCONFIGDIR="$PWD/build/matplotlib-cache"
mkdir -p "$MPLCONFIGDIR"
hostname > audit_node.txt
echo "$SLURM_JOB_ID" > audit_job.txt
extra=""
if [[ -f resume_job.txt ]]; then extra=",$(cat resume_job.txt)"; fi
if [[ -f boundary_job.txt ]]; then extra="$extra,$(cat boundary_job.txt)"; fi
sacct -j "16371528,16371535,16371543,16371564,$SLURM_JOB_ID$extra" --format=JobID,JobName,State,ExitCode,Elapsed,AllocCPUS,NodeList -P > accounting.txt
python --version > analysis_environment.txt
module list >> analysis_environment.txt 2>&1
python repo/model/p2_h_ale_coupling/verify_adaptive.py
python repo/model/p2_h_ale_coupling/plot_adaptive.py
python repo/model/p2_h_ale_coupling/collect_adaptive.py
