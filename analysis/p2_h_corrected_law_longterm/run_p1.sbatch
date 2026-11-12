#!/bin/bash
#SBATCH --job-name=p2h26-p1
#SBATCH --partition=cnmix
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=0-23%12
#SBATCH --time=01:00:00
#SBATCH --output=logs/p1.%A_%a.out
#SBATCH --error=logs/p1.%A_%a.err
set -euo pipefail
test -n "${SLURM_JOB_ID:-}"
case "$(hostname)" in ibcln*|test*|ibtest*) exit 90;; esac
cd /WORK/hulm_work/zhixin/p2_h_corrected_law_longterm_20261112
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python3 - <<'PY'
import csv
def read(n):
    with open('raw/u0_n%d_p2_r1e-08_dt30_history.csv'%n) as f:
        return {float(r['t']):r for r in csv.DictReader(f)}
a,b=read(256),read(512)
error=max(abs(float(r['v'])/float(b[t]['v'])-1) for t,r in a.items() if t>=1 and t in b)
assert 180000 in a and 180000 in b, 'P0 incomplete'
assert error<.01, ('P0 spatial gate',error)
assert max(abs(float(r[k])) for r in b.values() for k in ['salt_res','water_res','diss_res','local_res'])<2e-9
print('P0 gate passed; common-time max speed difference:',error)
PY
line=$((SLURM_ARRAY_TASK_ID+2))
row=$(sed -n "${line}p" cases_p1.csv)
IFS=, read -r id n power rtol step u beta av jfactor <<< "$row"
./longterm "$id" "$n" "$power" "$rtol" "$step" "$u" "$beta" "$av" "$jfactor" forcing.csv
