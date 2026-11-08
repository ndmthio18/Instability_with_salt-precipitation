#!/bin/bash
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:10:00
#SBATCH --job-name=p2h-crystal-pack
#SBATCH --output=package-%j.out
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || exit 90
case "$(hostname)" in ibcln*|login*|test*) exit 91;; esac
cd /WORK/hulm_work/zhixin/p2_h_crystal_boundary_20261108
hostname
python3 repo/model/p2_h_ale_coupling/collect_crystal.py
