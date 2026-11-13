#!/bin/bash
#SBATCH --job-name=p2h28-delta-pack
#SBATCH --partition=cnmix
#SBATCH --ntasks=1
#SBATCH --time=00:10:00
#SBATCH --output=logs/deltapack.%j.out
#SBATCH --error=logs/deltapack.%j.err
set -euo pipefail
cd /WORK/hulm_work/zhixin/p2_h_double_front_audit_20261112
test -n "${SLURM_JOB_ID:-}"
module purge
module load compilers/gcc/v12.2.0
hostname > provenance/final_environment.txt
g++ --version >> provenance/final_environment.txt
uname -a >> provenance/final_environment.txt
sacct -S 2026-11-12 -u hulm -X -P -o JobID,JobName,State,ExitCode,Start,End,Elapsed,AllocCPUS,NodeList | awk -F'|' 'NR==1 || $2 ~ /^p2h28/' > provenance/all_task28_jobs.csv
tar -czf final_delta.tar.gz raw/fast_grid_x008_n512_history.csv raw/fast_grid_x008_n512_diagnostics.csv raw/fast_grid_x008_n512_status.txt profiles_selected/fast_grid_x008_n512_profiles.csv raw/parallel* raw/stress_repro* raw/grid_x008_n512_history.csv raw/grid_x008_n512_diagnostics.csv raw/grid_x008_n512_status.txt logs source_sqrt_v7 source_parallel_v8 provenance *gate*json *gate.txt *sha256.txt final_completion.json final_light.sha256
sha256sum final_delta.tar.gz > final_delta.sha256
