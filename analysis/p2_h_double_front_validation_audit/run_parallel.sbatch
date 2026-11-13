#!/bin/bash
#SBATCH --job-name=p2h28-parallel-x008
#SBATCH --partition=cnmix
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --output=logs/parallel_x008.%j.out
#SBATCH --error=logs/parallel_x008.%j.err
set -euo pipefail
test -n "${SLURM_JOB_ID:-}"
case "$(hostname)" in ibcln*|test*|ibtest*) exit 90;; esac
cd /WORK/hulm_work/zhixin/p2_h_double_front_audit_20261112
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python3 -c 'import json; a=json.load(open("fast_p01_gate.json")); b=json.load(open("p2_kinetic_gate.json")); assert a["P0_pass"] and a["P1_pass"] and b["passed"]'
grep -q '^PASS bitwise' parallel_gate.txt
module purge
module load compilers/gcc/v12.2.0
export OMP_NUM_THREADS=8 OMP_PROC_BIND=close OMP_PLACES=cores FRESH_JAC=setup
./run_case_parallel_v8 parallel_x008_n512 512 128 1 1e-9 1e-9 60 -1 160951.5 forcing.csv
