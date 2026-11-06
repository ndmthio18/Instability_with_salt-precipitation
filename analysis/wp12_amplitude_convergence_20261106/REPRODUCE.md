# Reproduction notes

## Frozen scope

Case 36877 uses the A28 eigenprofile, the physical wavelength
`0.021428571428571432 m`, a `200`-cell vertical grid, horizontal grids R24/R48/R96,
and `dt = 113.68588647695695 s`.  Every new run was authorized only through 6T
(`7200` accepted steps; `818538.38263409003 s`) and retained the surface-solubility
stopping condition.

The two registered protocols are `protocol.json` in this directory and
`../wp13_amplitude_1e5_20261106/protocol.json`.

## Immutable provenance

- Executable SHA-256: `112076759fdff0942000307f729c2bd3c7e96963229b43c1b67c1f12bbf9df81`
- `main.cc` SHA-256: `dc7c629220fd0f3b23a6fcb6f156537d0d1d38b376550606235aeda7729dd2b4`
- A28 profile SHA-256: `d721e605d126e555f92b43200aa73b999ed4e2b61397edd1e98c0d9fddb36575`
- Persistent remote roots:
  - `A=1e-4`: `/WORK/hulm_work/zhixin/paper2_wp12_amplitude_20261106`
  - `A=1e-5`: `/WORK/hulm_work/zhixin/paper2_wp13_amplitude_1e5_20261106`

The full field and product archives remain under those roots.  Their SHA-256
values are recorded per grid in `metrics_manifest.json`, `finished.json`, and
`archive_acceptance.json`.  The local package contains the verified reductions,
not duplicate full-field archives.

## Rebuild the separated comparisons

From `paper/2`:

```powershell
python scripts/wp12_compare.py analysis/wp12_amplitude_convergence_20261106/comparison `
  A1e6_R24=analysis/wp12_amplitude_convergence_20261106/reference_metrics/R24 `
  A1e6_R48=analysis/wp12_amplitude_convergence_20261106/reference_metrics/R48 `
  A1e6_R96=analysis/wp12_amplitude_convergence_20261106/reference_metrics/R96 `
  A1e5_R24=analysis/wp13_amplitude_1e5_20261106/metrics/R24 `
  A1e5_R48=analysis/wp13_amplitude_1e5_20261106/metrics/R48 `
  A1e5_R96=analysis/wp13_amplitude_1e5_20261106/metrics/R96 `
  A1e4_R24=analysis/wp12_amplitude_convergence_20261106/enhanced_metrics/R24 `
  A1e4_R48=analysis/wp12_amplitude_convergence_20261106/enhanced_metrics/R48 `
  A1e4_R96=analysis/wp12_amplitude_convergence_20261106/enhanced_metrics/R96
```

This command writes same-grid amplitude sensitivity and identical-amplitude grid
sensitivity to different CSV files.  It never pools the two comparisons.

## Rebuild figures and validate the package

Add the Nature-figure audit scripts to `PYTHONPATH`, then run:

```powershell
$env:PYTHONPATH = 'C:\Users\12608\.agents\skills\nature-figure\scripts'
python scripts/wp12_make_figures.py
python scripts/verify_wp12_delivery.py
```

The plotting script uses physical time and the preregistered `t-t_nl` coordinate.
It contains no fitted time shift.  PDF text and collision audits are recorded
beside each figure and are part of the delivery gate.

The script also renders 600 dpi TIFF files.  They are deliberately not retained
in the lightweight Git delivery because the five deterministic TIFFs add about
216 MiB; rerunning the command above regenerates them byte-equivalently from the
tracked reductions and plotting source.

## Interpretation boundary

This package is an accelerated nonlinear convergence benchmark.  It neither
replaces the original `A=1e-6` trajectories nor validates an original or shifted
solubility-crossing time.  All comparisons stop at the actual common accepted
window, here 0--6T.
