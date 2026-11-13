#!/bin/bash
#SBATCH --job-name=p2h28-stress-repro
#SBATCH --partition=cnmix
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=0-1%2
#SBATCH --time=00:10:00
#SBATCH --output=logs/stressrepro.%A_%a.out
#SBATCH --error=logs/stressrepro.%A_%a.err
set -euo pipefail
cd /WORK/hulm_work/zhixin/p2_h_double_front_audit_20261112
test -n "${SLURM_JOB_ID:-}"
module purge
module load compilers/gcc/v12.2.0
deps=/WORK/hulm_work/zhixin/p2_h_ale_coupling_20261107/dependencies/install
src=source_double_v1
if [ "$SLURM_ARRAY_TASK_ID" = 1 ]; then src=source_longdouble_v2; fi
id=stress_repro_${SLURM_ARRAY_TASK_ID}
g++ -O2 -std=c++17 -Wall -Wextra "$src/verify.cc" -I"$src" -I"$deps/include" -L"$deps/lib64" -lsundials_cvode -lsundials_nvecserial -lsundials_sunmatrixdense -lsundials_sunlinsoldense -lsundials_generic -lm -o "verify_$id"
sha256sum "$src"/* "verify_$id" > "${id}_sha256.txt"
set +e
DIAG_ID="$id" ./"verify_$id"
code=$?
printf 'exit=%s\n' "$code" > "${id}_status.txt"
exit 0
