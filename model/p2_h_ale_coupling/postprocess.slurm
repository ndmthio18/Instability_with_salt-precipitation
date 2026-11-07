#!/bin/bash
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --job-name=p2h-audit
#SBATCH --output=postprocess-%j.out
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || exit 90
case "$(hostname)" in ibcln*|login*|test*) exit 91;; esac
module purge
module load soft/anaconda3/config
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /WORK/hulm_work/zhixin/p2_h_ale_coupling_20261107
echo "$SLURM_JOB_ID" > postprocess_job.txt
hostname > postprocess_node.txt
python -c 'import sys,numpy,scipy,matplotlib; print(sys.version); print(numpy.__version__,scipy.__version__,matplotlib.__version__)'
sacct -j 16371392,16371399,16371406,16371415,16371422,16371429,16371431 --format=JobID,State,ExitCode,Elapsed,AllocCPUS,NodeList -P > accounting.txt
python repo/model/p2_h_ale_coupling/verify_soil.py b
python repo/model/p2_h_ale_coupling/verify_soil.py c
python repo/model/p2_h_ale_coupling/plot_results.py
python repo/model/p2_h_ale_coupling/collect.py
