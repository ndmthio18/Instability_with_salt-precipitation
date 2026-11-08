# Reproduction and provenance

Base commit: `301097e3d2de3541bf63741cd729bd5cda22329a`.
Development branch: `research/p2-h-grid-precipitation-fix`.
Remote: `https://github.com/ndmthio18/Instability_with_salt-precipitation.git`.
No merge into the original research branches.

## Actual environment

HPCC root:
`/WORK/hulm_work/zhixin/p2_h_grid_precipitation_fix_20261108`.
Diagnostic build: `build/soil_solver`; corrected build:
`correction/build/soil_solver`. Production runs are serial DuMux processes,
parallelized across independent cases in Slurm. `DUMUX_NUM_THREADS=1`,
`OMP_NUM_THREADS=1`, `OPENBLAS_NUM_THREADS=1`, `MKL_NUM_THREADS=1`.

Compiler modules: `compilers/gcc/v12.2.0`, `mpi/openmpi-gcc12/v4.1.4`,
`tools/cmake/v3.25.2`. Analysis: `soft/anaconda3/config`.
SUNDIALS static libraries reuse the read-only installation under
`/WORK/hulm_work/zhixin/p2_h_ale_coupling_20261107/dependencies/install`.

Exact DuMux/DUNE source stack:
`/WORK/hulm_work/zhixin/Greal81_coupling_regime_20261106/staging/exact_stack_sources_r04_full.tar.gz`,
SHA256 `2a29d942b87c157eb0bd58dd35577c1088932a048411e069837bf36b80cad684`.
The build scripts verify this before unpacking to compute-node scratch.
Frozen P1-R2 sources are copied into an isolated application; originals are not edited.

Binary SHA256:

| Purpose | SHA256 |
|---|---|
| Diagnostic / projection control | `260119a917c9e6fddb2a6030d7ee5306ae476ba5bf06dc25c56e7df11582df73` |
| Corrected / both selectable modes | `b64d3c453c24f2e195f71f70560cf6ccc7db4b7221fc474bb81110e347a2fa72` |

Actual compiled main-file hashes, ablation binary hash, source archive hashes
and protected-source checks are in `provenance.json`. `runtime_source_sha256.csv`
records the two isolated runtime trees separately: their diagnostic headers
intentionally differ because the correction was introduced after diagnosis.

## Actual scheduled stages

| Job | Purpose | Outcome |
|---|---|---|
|16380079|Historical inventory audit, first diagnostic evidence|Completed|
|16380077|Diagnostic compile, donor unit test, two smoke cases|Completed|
|16380103|16 original / controlled precrust cases|All completed 2500 s|
|16380104|Original 160/320 activation diagnostics|160 completed 710 s; 320 reproduced physical-bound failure|
|16380113|Independent corrected compile and two smoke cases|Completed|
|16380117|12 conservative ALE interventions|All requested horizons completed|
|16380120|Three inactive/legacy regression cases|Completed|
|16380119|15 corrected precrust / mesh / timestep cases|All completed 2500 s|
|16380200|Two ablation timestep refinements|Completed|
|16380222|Four coupled grids, two dt/2, one forced rollback|All seven completed 2500 s|
|16380240|Full numerical audit and first figure rendering|Completed 0:0|
|16380253|Initial archive and additional exact regressions|Completed|
|16380264|Figure QA revision and lightweight final archive|Completed|

Intermediate audits are also retained in `accounting.txt` and `logs/`.
Compilation and calculations were performed on Slurm compute nodes, including
ibc11b11n04, ibc11b11n19 and ibc12b11n38. No solver or postprocessor ran on login nodes.

## Run sequence

The scripts use explicit working directories. To rerun, create a **new** HPCC
directory and change only the administrative working-directory paths in copies
of the scripts. Case folders refuse overwrite. Do not rerun into the archived
directories or replace a binary used by an active process.

The actual sequence used these entry points, under the recorded root:

```bash
sbatch repo/model/p2_h_ale_coupling/build_grid.slurm
sbatch repo/model/p2_h_ale_coupling/matrix_grid.slurm diagnostic
sbatch -c2 repo/model/p2_h_ale_coupling/matrix_grid.slurm ale_diagnostic
# Independent correction directory and source snapshot, after A-C evidence:
sbatch repo/model/p2_h_ale_coupling/build_fixed.slurm
# Run from correction/ after its smoke/budget audit:
sbatch repo/model/p2_h_ale_coupling/matrix_fixed.slurm fixed_pre
sbatch -c3 repo/model/p2_h_ale_coupling/matrix_fixed.slurm fixed_tests
# Only after the complete four-grid precrust audit:
sbatch repo/model/p2_h_ale_coupling/matrix_fixed.slurm fixed_coupled
# Analysis root:
sbatch repo/model/p2_h_ale_coupling/audit_grid.slurm final
sbatch repo/model/p2_h_ale_coupling/package_grid.slurm
sbatch repo/model/p2_h_ale_coupling/qa_grid.slurm
```

Standalone interventions were submitted with Slurm `--wrap` commands loading
GCC and invoking `run_ablations.py` / `refine_ablations.py`. Baseline and partial
audits loaded Anaconda and invoked `analyze_grid.py` / `analyze_ale_failure.py`
through Slurm. These small jobs are listed in `accounting.txt`; their stdout is
retained in `logs/`. Every simulation's exact argument vector and input file is
under `cases/<stage>/<case>/command.json` and `params.input`.

Corrected production options explicitly include:

```text
-P2H.PrecipitationMode implicit_surface
-P2H.DonorMode adaptive_multilayer
-P2H.MaxDonorDepth 0.005
-P2H.CrustCells 64
-P2H.GridDiagnostics true
-TimeLoop.DtInitial 0.5
-TimeLoop.MaxTimeStepSize 0.5
-TimeLoop.TEnd 2500
```

The default activation threshold is 3.5e-6 m and is deliberately omitted from
the matrix command line to exercise the default. The retained saturation is
0.99; no command lowers it. Time-halving cases use 0.25 s for both timestep
parameters. Precrust-only cases additionally select `P2H.PrecrustOnly=true`.
The legacy compatibility fallback is `projection`; corrected runners always
select `implicit_surface` explicitly.

Family A uses `Grid.Grading1=-0.97**(80/N)` evaluated before launching DuMux.
Family B uses `Grid.Positions1="0 0.198 0.2"`, `Grid.Cells1="Ndeep Ntop"`,
`Grid.Grading1="1 1"`. Actual SCV-derived dimensions are in the diagnostics;
no uniform-spacing assumption enters inventory or donor allocation.

## Archive and lightweight data policy

Full original CSV and solver logs (HPCC only):
`full_grid_evidence.tar.gz`, SHA256
`54902318559c5b8405ea991ac39f2e28ba5faf50d2ab24d6d3b42e24b81e6412`.

Final compact archive: `grid_delivery_v2.tar.gz`, SHA256
`dc89782de22fa6b4bcd05f8bd711ff1328d4035939e5f7296d3e9392034f36fb`.
The local external archive is in `../p2_h_grid_hpcc_artifacts/`, outside Git.
No binary, raw VTK or bulk raw-history archive is committed.

`raw_output_sha256.csv` identifies complete remote numerical outputs.
The `cases/` histories are explicitly sampled (every 100th accepted record,
activation, changed rejection count and final state); they are **not** passed
off as the full histories used for verification. Full accepted raw histories
were audited on HPCC before sampling. Selected soil fields retain time zero,
activation, 1200/1900/2500 s where actually available, and last saved state.
ALE fields retain every accepted step in the first ten seconds after activation,
100 s snapshots, last accepted state and all bound-rejected output states.

The v2 packaging verifies that all raw case-file hashes match the first package
before reusing the full archive. Only figure display and compact selection
changed. Whole-run raw output checksums did not change during visual QA.
Some original P1-R2 diagnostics describe **soil-only** inventories; they omit
the external solid/ALE inventory and must not be interpreted as combined
budget checks. They remain in the full archive. Combined checks are provided
by `independent_mass_budgets.csv` and `independent_residual_series.csv`.

## Failures and limitations retained

- Original320 ALE: exit1 at1417.629079 s, epsilon bound; this is the scientific
  failed control, not a missing or selectively excluded result.
- First partial audit16380204: failed because the not-yet-launched coupled
  directory was absent. The audit now labels absent stages NOT TESTED; later
  full audits passed. No physical calculation was changed or rerun to fix it.
- Explicit8G analysis submission: rejected before job creation because this
  cluster advertises1 MB Slurm MEMORY metadata per node. The supported CPU-only
  allocation was used; no physical parameter or tolerance changed.
- A transient SSH banner timeout during staging was retried; no running solver
  was restarted or canceled. Tar source timestamps produced harmless clock/old
  timestamp warnings; all numerical identity checks use content SHA256.
- Figure QA changed the failure-field color scale to expose both initial and
  subsequent concentrations/reactions, and removed an overlapping axis offset.
  No numerical result or acceptance threshold changed.

`local_delivery.sha256` is generated after adding the final local Chinese report
and reproduction notes. It covers the committed scientific delivery, excluding
itself. The remote compact archive checksum remains an independent transport
check even though local narrative/provenance annotations are added afterward.
