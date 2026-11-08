#!/bin/bash
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=00:40:00
#SBATCH --job-name=p2h-av-audit
#SBATCH --output=audit-%j.out
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || exit 90
case "$(hostname)" in ibcln*|login*|test*) exit 91;; esac
module purge
module load soft/anaconda3/config
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /WORK/hulm_work/zhixin/p2_h_av_sensitivity_20261108
python3 repo/model/p2_h_ale_coupling/analyze_av.py
python3 repo/model/p2_h_ale_coupling/verify_av.py
python3 repo/model/p2_h_ale_coupling/plot_av.py
find primary_exact dt_half_exact ale_refine_exact ale256_exact smoke_exact mechanism_exact unit_exact stage_a_exact primary dt_half ale_refine smoke_r02 mechanism_r02 unit stage_a -type f -print0 | sort -z | xargs -0 sha256sum > analysis/raw_output_sha256.txt
sha256sum build_exact/soil_solver build_exact/av_parameter_test_r02 build_exact/av_mechanism > analysis/executable_sha256.txt
hostname > analysis/audit_node.txt
echo "$SLURM_JOB_ID" > analysis/audit_job.txt
