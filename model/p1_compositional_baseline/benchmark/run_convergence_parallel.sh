#!/usr/bin/env bash
set -u

if (( $# != 4 )); then
    echo "usage: $0 BINARY INPUT EIGENPROFILE NEW_OUTPUT_ROOT" >&2
    exit 64
fi

binary=$1
input=$2
profile=$3
root=$4
if [[ -e $root ]]; then
    echo "refusing to reuse output root: $root" >&2
    exit 65
fi
mkdir -p "$root"

export OMPI_MCA_pml=ob1
export OMPI_MCA_btl=self,tcp
export OMPI_MCA_mtl=^ofi
export UCX_TLS=tcp,self
export OMP_NUM_THREADS=1
export TBB_NUM_THREADS=1

run_case() {
    local label=$1 nx=$2 nz=$3 dt=$4 output="$root/$1"
    mkdir "$output"
    cd "$output" || return
    printf 'START %s %s\n' "$label" "$(date -Is)" | tee start.txt
    /usr/bin/time -v nice -n 15 "$binary" "$input" \
        -Problem.Name "$label" \
        -Problem.EigenprofilePath "$profile" \
        -Grid.Cells0 "$nx" -Grid.Cells1 "$nz" \
        -TimeLoop.TEnd 60000 \
        -TimeLoop.DtInitial "$dt" \
        -TimeLoop.MaxTimeStepSize "$dt" \
        -TimeLoop.OutputEverySteps 50 \
        -Dune.NumThreads 1 > run.log 2> time.log
    local status=$?
    printf '%s\n' "$status" > process_exit_code.txt
    if (( status == 0 )); then
        sha256sum run_status.json component_budget.csv model_metadata.json > output_sha256.txt
        printf 'DONE %s %s\n' "$label" "$(date -Is)" | tee done.txt
    else
        printf 'FAILED %s exit=%s %s\n' "$label" "$status" "$(date -Is)" | tee failed.txt
    fi
    return "$status"
}

run_case B1_hcoarse_n140_nz100_dt227_t60000 140 100 227.3717729539139 &
pid_h=$!
run_case B1_vcoarse_n280_nz50_dt227_t60000 280 50 227.3717729539139 &
pid_v=$!
run_case B1_dt_half_n280_nz100_dt113_t60000 280 100 113.68588647695695 &
pid_t=$!

status=0
wait "$pid_h" || status=1
wait "$pid_v" || status=1
wait "$pid_t" || status=1
printf 'parallel convergence exit=%s %s\n' "$status" "$(date -Is)" | tee "$root/parallel_status.txt"
exit "$status"
