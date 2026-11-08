# Reproduction

## Branch and commits

- Branch: `research/p2-h-kr-av-reservoir-x008`
- Scenario commit: `b9815894aac8db30b1a780a031667a9e7cb69b31`
- Bound diagnostic commit: `79b9e9980d4df2fb615e0e020ae8cb8b4d177abb`
- Postprocessing commit: `cdf9a7908299598f3d7534c8b3ebc5fff9d49b8a`

The default frozen reproduction remains `x_NaCl = 0.05`. This scenario is
selected explicitly with `P2H_SALT_MOLE_FRACTION=0.08`; both the initial field
and bottom reservoir use that value.

## HPCC execution order

From a login shell, submit jobs to Slurm; the scripts run all numerical work on
compute nodes:

```bash
sbatch model/p2_h_ale_coupling/verify_reservoir_x008.slurm
sbatch model/p2_h_ale_coupling/pilot_reservoir_x008.slurm
sbatch model/p2_h_ale_coupling/diagnose_x008_bound.slurm
sbatch model/p2_h_ale_coupling/postprocess_x008_fail.slurm
```

Observed jobs were 16437446, 16437448, 16438191, and 16438205. The pilot is
expected to exit nonzero at the ALE upper-porosity bound. The diagnostic job
passes only when that expected solver failure is reproduced and both the last
accepted and first rejected spatial states are found.

The HPCC root used here was:

`/WORK/hulm_work/zhixin/p2_h_kr_av_reservoir_x008_20261108`

## Integrity checks

Run SHA-256 checks from this analysis directory using `CHECKSUMS.sha256`. The
four `.gz` files should also be fully decompressed, not merely listed. Review
`verification.json` before interpreting plots. A completed solver run alone is
not the scientific gate; `requested_horizon_passed` must be true before the
matrix scripts are eligible.

