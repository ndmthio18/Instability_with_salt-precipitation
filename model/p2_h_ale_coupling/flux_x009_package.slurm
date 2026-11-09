#!/bin/bash
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:15:00
#SBATCH --job-name=p2h-x009-pack
#SBATCH --output=logs/package-%j.out
#SBATCH --error=logs/package-%j.err
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || exit 90
case "$(hostname)" in ibcln*|login*|test*) exit 91;; esac
ROOT=/WORK/hulm_work/zhixin/p2_h_flux_mechanism_x009_20261109
cd "$ROOT"
mkdir -p delivery
tar -czf delivery/p2_h_flux_mechanism_x009_compact.tar.gz \
  analysis_output base.input input_sha256_build4.txt run_input_sha256.txt \
  post2_input_sha256.txt unit/source_binary.sha256 \
  build_submission.txt build2_submission.txt build3_submission.txt build4_submission.txt \
  run_submission.txt post_submission.txt post2_submission.txt \
  logs/build-16445764.out logs/build-16445764.err \
  logs/build-16445772.out logs/build-16445772.err \
  logs/build-16445803.out logs/build-16445803.err \
  logs/build-16445808.out logs/build-16445808.err \
  logs/run-16445813_0.out logs/run-16445813_0.err \
  logs/run-16445813_1.out logs/run-16445813_1.err \
  logs/run-16445813_2.out logs/run-16445813_2.err \
  logs/post-16446287.out logs/post-16446287.err \
  runs/x008_progressive_reservoir_n80_a64_dt0p5_diagnostic/case_spec.json \
  runs/x008_progressive_reservoir_n80_a64_dt0p5_diagnostic/command.json \
  runs/x008_progressive_reservoir_n80_a64_dt0p5_diagnostic/params.input \
  runs/x008_progressive_reservoir_n80_a64_dt0p5_diagnostic/reservoir_config.csv \
  runs/x008_progressive_reservoir_n80_a64_dt0p5_diagnostic/outcome.json \
  runs/x008_progressive_reservoir_n80_a64_dt0p5_diagnostic/stop_record.json \
  runs/x009_progressive_reservoir_n80_a64_dt0p5/case_spec.json \
  runs/x009_progressive_reservoir_n80_a64_dt0p5/command.json \
  runs/x009_progressive_reservoir_n80_a64_dt0p5/params.input \
  runs/x009_progressive_reservoir_n80_a64_dt0p5/reservoir_config.csv \
  runs/x009_progressive_reservoir_n80_a64_dt0p5/outcome.json \
  runs/x009_progressive_reservoir_n80_a64_dt0p5/stop_record.json \
  runs/x009_progressive_reservoir_n80_a64_dt0p25_timecheck/case_spec.json \
  runs/x009_progressive_reservoir_n80_a64_dt0p25_timecheck/command.json \
  runs/x009_progressive_reservoir_n80_a64_dt0p25_timecheck/params.input \
  runs/x009_progressive_reservoir_n80_a64_dt0p25_timecheck/reservoir_config.csv \
  runs/x009_progressive_reservoir_n80_a64_dt0p25_timecheck/outcome.json \
  runs/x009_progressive_reservoir_n80_a64_dt0p25_timecheck/stop_record.json
sha256sum delivery/p2_h_flux_mechanism_x009_compact.tar.gz \
  > delivery/p2_h_flux_mechanism_x009_compact.tar.gz.sha256
