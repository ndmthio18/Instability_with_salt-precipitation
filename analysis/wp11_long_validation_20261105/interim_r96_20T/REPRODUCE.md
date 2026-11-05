# Reproduce the interim assessment

From `paper/2`, first obtain the archive named in `ARCHIVE_POINTER.json`, verify its SHA-256, and extract it so that the snapshot root contains `manifest.json`.

```powershell
python scripts/wp11_scalar_metrics.py <snapshot-root> analysis/wp11_long_validation_20261105/interim_r96_20T/R96 R96 --snapshot-manifest <snapshot-root>/manifest.json
$env:PYTHONPATH='scripts'
python scripts/wp11_interim_field_metrics.py <snapshot-root> analysis/wp11_long_validation_20261105/interim_r96_20T/R96 R96
python scripts/wp11_interim_compare.py analysis/wp11_long_validation_20261105/interim_r96_20T/comparison --end-T 20.035 --end-step 24042 R24 analysis/wp11_long_validation_20261105/stage_a/R24 analysis/wp10_r24_long_20261104/run/transport_budget.csv R48 analysis/wp11_long_validation_20261105/stage_a/R48 analysis/wp11_long_validation_20261105/R48_run/transport_budget.csv R96 analysis/wp11_long_validation_20261105/interim_r96_20T/R96 <snapshot-root>/transport_budget.csv
python scripts/wp11_interim_figures.py analysis/wp11_long_validation_20261105/interim_r96_20T/figures analysis/wp11_long_validation_20261105/stage_a/R24 analysis/wp11_long_validation_20261105/stage_a/R48 analysis/wp11_long_validation_20261105/interim_r96_20T/R96 analysis/wp10_r24_long_20261104/run analysis/wp11_long_validation_20261105/R48_run <snapshot-root> --end-T 20.035
```

The raw archive and VTUs are intentionally not committed. The immutable manifest, archive pointer, reductions and scripts are committed; `integrity.json` records all hash, row-continuity, timestamp, conservation and field-mapping checks.
