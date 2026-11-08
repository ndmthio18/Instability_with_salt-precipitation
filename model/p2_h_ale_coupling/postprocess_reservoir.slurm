#!/bin/bash
#SBATCH --partition=cnmix
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --job-name=p2h-res-post
#SBATCH --output=logs/post-%j.out
#SBATCH --error=logs/post-%j.err
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || exit 90
case "$(hostname)" in ibcln*|login*|test*) exit 91;; esac
module purge
module load compilers/gcc/v12.2.0 mpi/openmpi-gcc12/v4.1.4 tools/cmake/v3.25.2
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 DUMUX_NUM_THREADS=1
cd /WORK/hulm_work/zhixin/p2_h_kr_av_reservoir_20261108
mkdir -p analysis/source analysis/job_logs analysis/verification_evidence analysis/build
cp repo/model/p2_h_ale_coupling/{README.md,soil_methods.inc,RESERVOIR_FORMULATION.md,run_reservoir.py,verify_reservoir.py,run_old_noflux.py,verify_ale_pilot.py,analyze_reservoir.py,plot_reservoir.py,build_reservoir.slurm,verify_reservoir.slurm,pilot_reservoir.slurm,ale_pilot_reservoir.slurm,matrix_reservoir.slurm,postprocess_reservoir.slurm} analysis/source/
cp logs/* analysis/job_logs/ 2>/dev/null || true
cp verification_checks.json pilot_gate.json ale_pilot_gate.json analysis/verification_evidence/
cp unit/* analysis/build/
cp build/soil_build.log analysis/build/
cp repo/analysis/p2_h_kr_av_comparison_20261108/results/summary.csv analysis/results_old_no_flux.csv
set +e
python3 repo/model/p2_h_ale_coupling/analyze_reservoir.py > analysis/analyze_stdout.json
ANALYSIS_STATUS=$?
set -e
python3 repo/model/p2_h_ale_coupling/plot_reservoir.py
python3 - <<'PY' > analysis/figure_checks.json
import json, pathlib
from PIL import Image
root=pathlib.Path('analysis/figures')
checks=[]
for stem in ('figure1_porosity_cloud_kr','figure2_porosity_cloud_av','figure3_evaporation','figure4_height_porosity_internal_fraction','figure5_reservoir_supply'):
    files={ext:root/f'{stem}.{ext}' for ext in ('png','pdf','svg')}
    image=Image.open(files['png']).convert('RGB')
    extrema=image.getextrema()
    checks.append({'stem':stem,'all_formats':all(p.exists() and p.stat().st_size>1000 for p in files.values()),'pixels':image.size,'nonblank':any(lo<hi for lo,hi in extrema),'sizes':{ext:p.stat().st_size for ext,p in files.items()}})
print(json.dumps({'passed':all(c['all_formats'] and c['nonblank'] for c in checks),'checks':checks},indent=2))
PY
find analysis -type f ! -name CHECKSUMS.sha256 -print0 | sort -z | xargs -0 sha256sum > analysis/CHECKSUMS.sha256
printf '%s\n' "$ANALYSIS_STATUS" > analysis/scientific_gate_exit_code.txt
exit "$ANALYSIS_STATUS"
