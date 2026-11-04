#!/bin/bash
#SBATCH --job-name=WP9-threads
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:20:00
#SBATCH --output=stdout.%j
#SBATCH --error=stderr.%j
set -euo pipefail

# Shared-memory thread scaling of the frozen executable. Ten steps only; this is an
# infrastructure measurement, not a scientific run. The executable, input and
# eigenprofile are the ones already staged and hash-checked by the portability job.
if [ -z "${SLURM_JOB_ID:-}" ]; then echo "WP9: no SLURM_JOB_ID; refusing"; exit 2; fi
case "$(hostname)" in
  ibcln*|test0*|gln0*) echo "WP9: login or test node $(hostname); refusing"; exit 2;;
esac

threads="$1"
root="${SLURM_SUBMIT_DIR}"
work="${root}/threads_${threads}"
rm -rf "${work}"; mkdir -p "${work}"
cd "${work}"
ln -s ../bundle bundle
cp ../params.input ../eigenprofile.dat .

{
  echo "hostname=$(hostname)"
  echo "job=${SLURM_JOB_ID}"
  echo "threads=${threads}"
  echo "cpus_per_task=${SLURM_CPUS_PER_TASK:-}"
  echo "affinity=$(taskset -cp $$ 2>/dev/null | sed 's/.*: //')"
  echo "elf_sha256=$(sha256sum bundle/bin/paper2_transport | cut -d' ' -f1)"
} > environment.txt

start=$(date +%s.%N)
env DUMUX_NUM_THREADS="${threads}" OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  python3 bundle/launch_solver.py params.input > solver.log 2>&1
rc=$?
end=$(date +%s.%N)
echo "returncode=${rc}" >> environment.txt
echo "wall_s=$(echo "${end} - ${start}" | bc)" >> environment.txt
grep -m1 "Simulation took" solver.log >> environment.txt || true
grep "Assemble/solve/update time" solver.log | tail -5 >> environment.txt || true
exit ${rc}
