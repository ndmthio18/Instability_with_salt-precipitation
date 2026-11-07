#!/bin/bash
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --job-name=p2f_runs
#SBATCH --output=logs/run-%A_%a.out
#SBATCH --array=0-6
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || { echo "Slurm allocation required" >&2; exit 90; }
case "$(hostname)" in ibcln*|login*|test*) echo "Refusing login-node computation" >&2; exit 91;; esac
module purge
module load compilers/gcc/v12.2.0
module load mpi/openmpi-gcc12/v4.1.4
ROOT=/WORK/hulm_work/zhixin/p2_f_analytical_growth_20261107
SOLVER="$ROOT/build/paper2_transport_p2f"
INPUT="$ROOT/source/model/p2_f_analytical_growth/integration/p2f_1d.input"
[[ -x "$SOLVER" ]] || { echo "missing solver" >&2; exit 92; }
cases=(
  grid40_seed10_dt0025
  grid80_seed10_dt0025
  grid160_seed10_dt0025
  grid80_seed5_dt0025
  grid80_seed20_dt0025
  grid80_seed10_dt00125
  inadequate_supply
)
cells=(40 80 160 80 80 80 80)
seeds=(1e-5 1e-5 1e-5 5e-6 2e-5 1e-5 1e-5)
dts=(0.025 0.025 0.025 0.025 0.025 0.0125 0.025)
water=(0.01 0.01 0.01 0.01 0.01 0.01 0.25)
fluxes=(1e-6 1e-6 1e-6 1e-6 1e-6 1e-6 1.08e-5)
i=$SLURM_ARRAY_TASK_ID
name=${cases[$i]}
run="$ROOT/results/${P2F_RESULT_SET:-coupled}/$name"
[[ ! -e "$run" ]] || { echo "refusing existing run directory: $run" >&2; exit 93; }
mkdir -p "$run"
cd "$run"
set +e
"$SOLVER" "$INPUT" \
  -Grid.Cells1 "${cells[$i]}" \
  -P2F.SeedHeight "${seeds[$i]}" \
  -TimeLoop.DtInitial "${dts[$i]}" \
  -TimeLoop.MaxTimeStepSize "${dts[$i]}" \
  -Problem.WaterLevel "${water[$i]}" \
  -P2F.PrescribedEvaporationMassFlux "${fluxes[$i]}" \
  > solver.out 2>&1
status=$?
set -e
echo "$status" > exit_status.txt
if [[ "$name" == inadequate_supply ]]; then
  [[ "$status" -ne 0 ]] || { echo "inadequate supply was not rejected" >&2; exit 94; }
  grep -q "P2-F physical validity gate failed" solver.out || {
    echo "inadequate case failed for the wrong reason" >&2; exit 95; }
  echo P2F_EXPECTED_DOMAIN_EXIT
else
  [[ "$status" -eq 0 ]] || { tail -100 solver.out; exit "$status"; }
  echo P2F_COUPLED_COMPLETE
fi
