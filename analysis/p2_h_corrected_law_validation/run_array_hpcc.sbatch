#!/bin/bash
#SBATCH --job-name=p2h25-matrix
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:20:00
#SBATCH --output=logs/case_%A_%a.out
#SBATCH --error=logs/case_%A_%a.err

set -euo pipefail
if [[ -z "${SLURM_JOB_ID:-}" ]]; then echo "Refusing to compute outside Slurm" >&2; exit 90; fi
case "$(hostname)" in ibcln*|ibtest*|test*) echo "Refusing login/test-node work" >&2; exit 91;; esac
run_dir="/WORK/hulm_work/zhixin/p2_h_corrected_law_validation_20261112"
cd "$run_dir"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
line_number=$((SLURM_ARRAY_TASK_ID+2))
row=$(awk -F, -v line="$line_number" 'NR==line {print; exit}' configs/cases.csv)
if [[ -z "$row" ]]; then echo "Missing case row $line_number" >&2; exit 92; fi
row=${row%$'\r'}
IFS=, read -r id suite da u cells power rtol max_step_div end_tau outputs_per_tau mode <<< "$row"
build/ale_validation_driver "$id" "$suite" "$da" "$u" "$cells" "$power" \
  "$rtol" "$max_step_div" "$end_tau" "$outputs_per_tau" "$mode" results/raw
