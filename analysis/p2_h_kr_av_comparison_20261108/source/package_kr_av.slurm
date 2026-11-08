#!/bin/bash
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:20:00
#SBATCH --job-name=p2h-krav-pack
#SBATCH --output=logs/package-%j.out
#SBATCH --error=logs/package-%j.err
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || exit 90
case "$(hostname)" in ibcln*|login*|test*) exit 91;; esac
export PATH=/rmprog/slurm/v22.05.7/bin:$PATH
cd /WORK/hulm_work/zhixin/p2_h_kr_av_comparison_20261108
mkdir -p analysis/logs analysis/evidence/pilot analysis/evidence/production analysis/source
find analysis/evidence -type f \( -name donor_attempts.csv -o -name coupling_events.csv \) -delete
cp logs/*.out logs/*.err analysis/logs/ 2>/dev/null || true
cp input_sha256.txt *submission*.txt analysis/logs/ 2>/dev/null || true
cp unit/*.log unit/*.sha256 unit/node.txt unit/job.txt analysis/logs/ 2>/dev/null || true
for directory in pilot/pilot_ref_ale64 pilot/pilot_ref_ale128; do
    [ -d "$directory" ] || continue
    target="analysis/evidence/pilot/$(basename "$directory")"
    mkdir -p "$target"
    cp "$directory"/{case_spec.json,command.json,exit_status.txt,outcome.json,activation_config.csv,donor_profiles.csv,boundary_faces.csv,stop_record.json,transport_metadata.json,params.input} "$target"/ 2>/dev/null || true
    [ ! -f "$directory/donor_attempts.csv" ] || { head -n 2 "$directory/donor_attempts.csv"; tail -n 1 "$directory/donor_attempts.csv"; } > "$target/donor_attempts_first_last.csv"
    [ ! -f "$directory/coupling_events.csv" ] || { head -n 2 "$directory/coupling_events.csv"; tail -n 1 "$directory/coupling_events.csv"; } > "$target/coupling_events_first_last.csv"
done
for directory in production/* production_corrected/*; do
    [ -d "$directory" ] || continue
    group="$(basename "$(dirname "$directory")")"
    target="analysis/evidence/production/${group}/$(basename "$directory")"
    mkdir -p "$target"
    cp "$directory"/{case_spec.json,command.json,exit_status.txt,outcome.json,activation_config.csv,donor_profiles.csv,boundary_faces.csv,stop_record.json,transport_metadata.json,params.input} "$target"/ 2>/dev/null || true
    [ ! -f "$directory/donor_attempts.csv" ] || { head -n 2 "$directory/donor_attempts.csv"; tail -n 1 "$directory/donor_attempts.csv"; } > "$target/donor_attempts_first_last.csv"
    [ ! -f "$directory/coupling_events.csv" ] || { head -n 2 "$directory/coupling_events.csv"; tail -n 1 "$directory/coupling_events.csv"; } > "$target/coupling_events_first_last.csv"
done
cp repo/model/p2_h_ale_coupling/{analyze_kr_av.py,plot_kr_av.py,verify_kr_av.py,run_kr_av.py,build_kr_av.slurm,pilot_kr_av.slurm,pilot_128_kr_av.slurm,matrix_kr_av.slurm,postprocess_kr_av.slurm,package_kr_av.slurm} analysis/source/
IDS="$(grep -hE '^[0-9]+$' *submission*.txt | sort -u | paste -sd, -),${SLURM_JOB_ID}"
sacct -j "$IDS" --format=JobID,JobName,Partition,State,ExitCode,Elapsed,AllocCPUS,NodeList,MaxRSS -P > analysis/logs/sacct.txt
tar -czf kr_av_delivery.tar.gz analysis
sha256sum kr_av_delivery.tar.gz > kr_av_delivery.tar.gz.sha256
hostname > package_node.txt
printf '%s\n' "$SLURM_JOB_ID" > package_job.txt
