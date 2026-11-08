# Reproduction record

## Frozen scope

- Base commit: `301097e3d2de3541bf63741cd729bd5cda22329a`
- Branch: `research/p2-h-av-sensitivity`
- HPCC root: `/WORK/hulm_work/zhixin/p2_h_av_sensitivity_20261108`
- DuMuX/SUNDIALS dependencies: task-local links to the previously verified P2-H dependency tree
- Primary executable SHA-256: `a4b6d59563fbc1b8caa24a6a163766215dc1995bdcff126cf8902efa720ad8b5`
- Parameter-test SHA-256: `5e9f71f1227ab0b21860cc1cdbb8f2a3c5c2d7a74c4af7e5fa573f3ec7b34456`
- Standalone-mechanism SHA-256: `d49f49252db00a61e7f45a5a15cd80fc3a9db39970b2cc8e6f72d871e7c3df0f`

The scripts intentionally refuse to overwrite an existing stage directory. Reproduction therefore
requires a fresh task root or explicitly named new stage; do not delete accepted evidence to rerun in place.

## HPCC sequence

All compilation, numerical runs, analysis and plotting were executed on Slurm compute nodes, not login nodes.

```bash
cd /WORK/hulm_work/zhixin/p2_h_av_sensitivity_20261108

sbatch repo/model/p2_h_ale_coupling/build_av_exact.slurm
sbatch repo/model/p2_h_ale_coupling/resume_av_exact.slurm

sbatch repo/model/p2_h_ale_coupling/matrix_av.slurm primary_exact
sbatch repo/model/p2_h_ale_coupling/matrix_av.slurm dt_half_exact
sbatch repo/model/p2_h_ale_coupling/matrix_av.slurm ale_refine_exact
sbatch repo/model/p2_h_ale_coupling/matrix_av.slurm ale256_exact

sbatch repo/model/p2_h_ale_coupling/audit_av.slurm
sbatch repo/model/p2_h_ale_coupling/package_av.slurm
```

The final figure-only visual correction was run as Slurm job `16406497`; it changes no numerical data.

## Principal jobs

| Purpose | Job | Slurm state |
|---|---:|---|
| corrected smoke/resume | 16405818 | COMPLETED 0:0 |
| exact build/resume and regression | 16405951 | COMPLETED 0:0 |
| primary exact matrix | 16405993 | COMPLETED 0:0 |
| dt/2 exact matrix | 16405994 | COMPLETED 0:0 |
| 128-cell ALE matrix | 16405995 | COMPLETED 0:0 |
| 256-cell selected A3 | 16406204 | COMPLETED 0:0; physical case exit 1 at bound |
| final audit | 16406416 | COMPLETED 0:0 |
| evidence packaging | 16406479 | COMPLETED 0:0 |
| final Figure 9 render | 16406497 | COMPLETED 0:0 |

An individual A3 solver exit code 1 is the registered physical-state stop and is retained in
`results/failures.csv`; Slurm completion means the matrix driver recorded the outcome, not that A3 reached 2500 s.

## Verification

```bash
cat analysis/verification_final.json
cat analysis/results/summary.csv
cat analysis/results/timestep_comparison.csv
cat analysis/results/ale_resolution_common_state.csv
sha256sum -c delivery_archives.sha256
```

Expected key checks are: exact preactivation/activation identity, historical M3 regression,
parameter unit test, global water/salt gate and timestep 1% gate all `true`; `primary_all_reached_2500_s`
and `ale_resolution_gate` are intentionally `false`.

The HPCC evidence archives are:

- `full_av_evidence.tar.gz`: SHA-256 `1cd92d3c7b6798de7c7f77f849b281dedf6a988050f32bbae2e1934f106e1cdf`
- `av_delivery.tar.gz`: SHA-256 `973394ee7c31ebf2e61b17e54696153cdd880e7c1dd557192182a227d829c831`

The 620 MB full archive remains on HPCC and is not committed to Git. The Git delivery contains compact
accepted-state tables, selected profiles, logs, checksums, failed-state records and figures.

