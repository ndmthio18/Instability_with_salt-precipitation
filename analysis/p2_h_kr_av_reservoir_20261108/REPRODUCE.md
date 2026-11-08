# P2-H reservoir-fed comparison reproduction

## Frozen scope

- Branch: `research/p2-h-kr-av-reservoir`
- Base: `1e713ae2a41e6a8ec8ba1a42f332cbb74714074c`
- HPCC root: `/WORK/hulm_work/zhixin/p2_h_kr_av_reservoir_20261108`
- Soil/ALE: 80/64 cells; `dt_max=0.5 s`; requested horizon `180000 s`.
- Bottom: `Problem.BottomBoundaryMode=reservoir`,
  `Problem.WaterLevel=0.01 m`.
- Initial and reservoir salt: `Problem.SaltMoleFraction=0.05`.
- Zero initial solid/crust; conservative activation at real surface solid
  equivalent to `0.1 µm`; `epsilon_birth=0.43`; `Sw_min=0.99`.

All compilation, simulations, numerical analysis, and plotting run under Slurm
on compute nodes. Login-node activity was limited to file transfer, checksum,
scheduler query, and submission.

## Registered commands and dependency gates

From the HPCC root:

```bash
/rmprog/slurm/v22.05.7/bin/sbatch repo/model/p2_h_ale_coupling/build_reservoir.slurm
/rmprog/slurm/v22.05.7/bin/sbatch repo/model/p2_h_ale_coupling/verify_reservoir.slurm
/rmprog/slurm/v22.05.7/bin/sbatch repo/model/p2_h_ale_coupling/pilot_reservoir.slurm
/rmprog/slurm/v22.05.7/bin/sbatch --dependency=afterok:16430419 repo/model/p2_h_ale_coupling/ale_pilot_reservoir.slurm
/rmprog/slurm/v22.05.7/bin/sbatch --dependency=afterok:16430445 repo/model/p2_h_ale_coupling/matrix_reservoir.slurm
/rmprog/slurm/v22.05.7/bin/sbatch --dependency=afterany:16430446 repo/model/p2_h_ale_coupling/postprocess_reservoir.slurm
```

The dependency chain is scientific, not merely operational: the five-case
matrix cannot start unless the 50 h reference pilot accepts ALE activation and
the ALE-grid/time diagnostic completes conservatively.

## Job register

| Job | Purpose |
|---|---|
| `16430286` | build and compiled unit check |
| `16430327` | invalid verification wrapper; retained failure |
| `16430417` | corrected three-case boundary verification |
| `16430419` | gated 50 h reference pilot |
| `16430445` | ALE128 and `dt/2` diagnostic, `afterok:16430419` |
| `16430446_[0-4]` | five unique production cases, `afterok:16430445` |
| `16430448` | analysis/figures/package, `afterany:16430446` |
| `16433911` | first FAIL diagnostic post-process; wrong frozen-history path |
| `16433914` | corrected FAIL diagnostic package |
| `16433915` | compute-node archive/checksum |
| `16435674` | final plain-axis visual-QA regeneration |

Case-level `command.json`, `params.input`, outcome, solver output, reservoir
configuration, boundary-face audit, accepted-step budget, and checksums are
retained in this delivery or the stated HPCC raw root.

The 181 MB plain-text `solver.out` is stored in Git as the lossless
`logs/solver.out.gz`; `CHECKSUMS_HPCC.sha256` records the SHA-256 of the
decompressed original. The full HPCC evidence archive is retained outside the
Git worktree as `task12_HPCC_analysis_fail_89719b1c.tar.gz` with archive
SHA-256 `89719b1c07aa686bf635df3c55c51d07a53ee65e53d5afc1ce1133fc28a8f016`.

The pilot scientific gate is `FAIL`; jobs `16430445`, `16430446_[0-4]`, and
`16430448` were cancelled without solver execution. Reproducing the delivered
result therefore stops after the pilot and invokes
`postprocess_reservoir_fail.slurm`, not the formal matrix scripts.
