#!/bin/bash
#SBATCH --job-name=WP12-ref-R96
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=01:00:00
#SBATCH --output=reference_extract.stdout
#SBATCH --error=reference_extract.stderr
set -euo pipefail
export PATH=/rmprog/slurm/v22.05.7/bin:/usr/local/bin:/usr/bin:/bin
case "$(hostname)" in ibcln*|test0*|gln0*) echo "refusing login/test node"; exit 2;; esac
root="$1"
source=/WORK/hulm_work/zhixin/paper2_wp11_long_20261105/R96/fields.tar.gz
srun -n 1 -c 1 /WORK/hulm_work/ys/miniconda3/bin/python3 "$root/wp12_extract_reference.py" \
  "$source" "$root/reference_R96_0_6T_fields.tar.gz" "$root/reference_R96_0_6T_manifest.json"
