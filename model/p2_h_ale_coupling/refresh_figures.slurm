#!/bin/bash
#SBATCH --partition=cnmix
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:05:00
#SBATCH --job-name=p2h-figures
#SBATCH --output=figures-%j.out
set -euo pipefail
[[ -n "${SLURM_JOB_ID:-}" ]] || exit 90
case "$(hostname)" in ibcln*|login*|test*) exit 91;; esac
module purge
module load soft/anaconda3/config
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
cd /WORK/hulm_work/zhixin/p2_h_ale_coupling_20261107
export MPLCONFIGDIR="$PWD/build/matplotlib-cache"
mkdir -p "$MPLCONFIGDIR"
python repo/model/p2_h_ale_coupling/plot_results.py
python - <<'PY'
import hashlib,json,os,pathlib
root=pathlib.Path.cwd()
path=root/'deliverables/remote_manifest.json'
data=json.loads(path.read_text())
for name in ('plot_results.py','refresh_figures.slurm'):
    source=root/'repo/model/p2_h_ale_coupling'/name
    data[str(source.relative_to(root))]=hashlib.sha256(source.read_bytes()).hexdigest()
path.write_text(json.dumps(data,indent=2))
(root/'deliverables/figure_revision.json').write_text(json.dumps(dict(
    job=os.environ['SLURM_JOB_ID'],changes=['salinity axes expressed as C-Ceq in mg/kg brine',
    'porosity not plotted before crust exists','explicit magnitude of evaporation difference'],
    numerical_runs_repeated=False),indent=2))
PY
tar -czf compact_deliverables_r02.tar.gz deliverables
sha256sum full_accepted_histories.tar.gz compact_deliverables_r02.tar.gz > transfer_r02.sha256
