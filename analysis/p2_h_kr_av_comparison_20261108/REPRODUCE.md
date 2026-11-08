# P2-H `k_r` / `a_v` comparison reproduction

## Frozen interpretation

- Branch: `research/p2-h-kr-av-comparison`
- Base: `63f267f9eae85bf434a2ad9ed0804a63723e9478`
- `k_r,ref = 1.0e-3 m/s`
- `a_v = alpha_av * 36 * (1-epsilon) / d_sub`, `d_sub = 25e-6 m`
- Initial dissolved NaCl mole fraction: `Problem.SaltMoleFraction = 0.05`
- Zero initial solid and zero initial crust height.
- Conservative activation when the real surface solid inventory first corresponds to
  `h_eq = 0.1 um`; this is an activation threshold, not a pre-created seed.
- Soil grid: 80 vertical cells; ALE grid: 64 cells; maximum soil time step: 0.5 s.
- Requested horizon: 180000 s (50 h).

The shared reference case is calculated once and used in both one-factor groups.
There are five unique calculations for six plotted comparison entries.

## HPCC location

`/WORK/hulm_work/zhixin/p2_h_kr_av_comparison_20261108`

All compilation, simulation, analysis and plotting commands run under Slurm on
compute nodes. The login node was used only for transfer, checksums, scheduler
queries and submission.

## Commands

From the HPCC run root:

```bash
export PATH=/rmprog/slurm/v22.05.7/bin:$PATH
sbatch repo/model/p2_h_ale_coupling/build_kr_av.slurm
sbatch repo/model/p2_h_ale_coupling/pilot_kr_av.slurm
sbatch repo/model/p2_h_ale_coupling/pilot_128_kr_av.slurm
sbatch repo/model/p2_h_ale_coupling/matrix_kr_av.slurm
sbatch repo/model/p2_h_ale_coupling/postprocess_kr_av.slurm
sbatch repo/model/p2_h_ale_coupling/package_kr_av.slurm
```

The final `matrix_kr_av.slurm` invokes the `production_corrected` stage. The
earlier `production` directory is deliberately retained as a failed provenance
record because it used the wrong parameter namespace for the requested initial
concentration.

## Acceptance

Run:

```bash
python3 repo/model/p2_h_ale_coupling/verify_kr_av.py
```

This checks the actual accepted initial field (not only the command string),
parameter matrix, zero artificial seed, grid/time settings, activation transfer,
water/salt conservation, field cadence and all PNG/PDF/SVG deliverables.

Scheduler completion is not classified as scientific completion. A case passes
the requested horizon only if its last accepted physical time reaches 180000 s.
For this matrix, all five soil prefixes reach 180000 s but
`scientific_status=FAIL_NO_ALE_ACTIVATION`; the 92/92 verification result means
the evidence package is internally valid, not that the sensitivity objective passed.

## Executed HPCC jobs

| Job | Purpose | Result |
|---|---|---|
| 16407531 | first build | failed: omitted read-only P2-D header |
| 16407548 | corrected build and unit check | completed |
| 16407572 | concurrent 64/128 pilot | orchestration race; 64 numerical prefix retained |
| 16407589 | isolated 128 pilot | completed wrapper; protected solver stop retained |
| 16407598 | first five-case matrix | invalidated: wrong initial-condition namespace |
| 16408105_[0-4] | corrected five-case 50 h matrix | all completed; no ALE activation |
| 16411448 | first final plotting | cancelled after 219 MB vector-SVG QA failure |
| 16413799 | rasterized plotting and verification | completed, 92/92 |
| 16413869 | final axis-corrected plotting and verification | completed, 92/92 |

The compiled solver SHA-256 is
`bce211633378d0fc3ef92dde787dfd96ba33ca367202634848138ff5a28f54f0`.
Raw full-step histories remain on HPCC; Git contains compact time series,
compressed field evidence, figures, logs, checksums and failure records.
Final package job `16414012` completed on a compute node.  The retrieved archive
SHA-256 is `caf205809478c6d15b18962b177210c550f2cf0b7099411ab0822ff0c1fe0eb2`.
