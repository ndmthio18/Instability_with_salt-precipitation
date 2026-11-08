# Reproduction

## Identity

- Branch: `research/p2-h-reservoir-x008-h3p5`
- Registered physical-input commit: `90417f615dd57f1995e22944a8d0f01a59f0da09`
- Final diagnostic/postprocessing commit: `3e43bbf393248dc0f5bad8acdc3ccac247c5ad55`
- HPCC root: `/WORK/hulm_work/zhixin/p2_h_reservoir_x008_h3p5_20261108`

The 3.5 µm case requires both environment variables:

```bash
export P2H_SALT_MOLE_FRACTION=0.08
export P2H_ACTIVATION_HEIGHT_M=3.5e-6
```

Without the second variable the reservoir driver intentionally retains the
historical 0.1 µm scenario for exact reproduction.

## Slurm sequence

```bash
sbatch repo/model/p2_h_ale_coupling/pilot_reservoir_x008_h3p5.slurm
sbatch repo/model/p2_h_ale_coupling/diagnose_reservoir_x008_h3p5.slurm
sbatch repo/model/p2_h_ale_coupling/postprocess_h3p5_compare.slurm
```

Observed scientific jobs were 16438585 and 16438590. Final successful
postprocessing was job 16438615; jobs 16438601 and 16438603 are retained failed
attempts documented in `FAILURES.md`.

Check every file against `CHECKSUMS.sha256`, then fully decompress the three
gzip logs. `verification.json` is fail-closed: conservation and activation pass,
but `requested_horizon_passed` and `overall_passed` remain false.
