#!/bin/bash
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=00:30:00
#SBATCH --job-name=p2h-av-package
#SBATCH --output=package-%j.out
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || exit 90
case "$(hostname)" in ibcln*|login*|test*) exit 91;; esac
cd /WORK/hulm_work/zhixin/p2_h_av_sensitivity_20261108
test -f analysis/verification_final.json
python3 repo/model/p2_h_ale_coupling/collect_av.py
tar -czf full_av_evidence.tar.gz primary_exact dt_half_exact ale_refine_exact ale256_exact smoke_exact mechanism_exact \
  primary dt_half ale_refine smoke_r02 mechanism_r02 unit unit_exact stage_a stage_a_exact analysis \
  build-*.out smoke-*.out matrix-*.out build-exact-*.out exact-resume-*.out audit-*.out \
  source_av_r01.sha256 overlay_av_r02.sha256 overlay_av_r04.sha256 overlay_av_r05.sha256 \
  source_av_r06_exact.sha256 overlay_av_r07.sha256 overlay_av_r08_analysis.sha256 \
  overlay_av_r09_auditfix.sha256 overlay_av_r10_verifyfix.sha256 overlay_av_r11_ale256.sha256 \
  overlay_av_r12_commonstate.sha256 overlay_av_r13_figure9.sha256
tar -czf av_delivery.tar.gz delivery
sha256sum full_av_evidence.tar.gz av_delivery.tar.gz > delivery_archives.sha256
echo "$SLURM_JOB_ID" > delivery/package_job.txt
