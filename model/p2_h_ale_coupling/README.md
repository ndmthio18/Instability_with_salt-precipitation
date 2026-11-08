# P2-H C++ ALE and conditional DuMuX coupling

## Reservoir-fed long-time revision

`RESERVOIR_FORMULATION.md` documents the accepted-step bottom component-flux
ledger used by branch `research/p2-h-kr-av-reservoir`. Its registered sequence
is `build_reservoir.slurm`, `verify_reservoir.slurm`, the gated 50 h
`pilot_reservoir.slurm`, `ale_pilot_reservoir.slurm`, the five unique cases in
`matrix_reservoir.slurm`, and `postprocess_reservoir.slurm`. Simulation and
post-processing scripts require a Slurm allocation. Every reservoir case passes
`Problem.BottomBoundaryMode=reservoir` and
`Problem.SaltMoleFraction=0.05` explicitly.

The follow-up `research/p2-h-kr-av-reservoir-x008` keeps the same compiled
physics and accepts the initial/reservoir mole fraction through the explicit
`P2H_SALT_MOLE_FRACTION` run environment. The registered x=0.08 Slurm scripts
set this value for both the accepted initial soil field and the bottom
Dirichlet reservoir; the default remains 0.05 for frozen-task reproduction.

The isolated `research/p2-h-av-sensitivity` revision adds an internal-only
`P2H.AvMultiplier`, exact accepted-state precipitation/dissolution BDF ledgers,
and registered 80/160-soil plus 64/128/256-ALE sensitivity checks at fixed
`kr=1e-3 m/s`. See the [parameter contract](AV_SENSITIVITY.md),
[Chinese scientific report](../../analysis/p2_h_av_sensitivity_20261108/REPORT_ZH.md),
[failure record](../../analysis/p2_h_av_sensitivity_20261108/FAILURES.md), and
[reproduction record](../../analysis/p2_h_av_sensitivity_20261108/REPRODUCE.md).
The `10av` local-porosity/stop-time ALE refinement gate fails; this branch is
not a calibrated physical salt-crust model.

The `research/p2-h-adaptive-donor` extension adds optional conservative,
depth-limited multi-layer pore charging. See [registered experiments](ADAPTIVE_PLAN.md),
[allocation equations and assumptions](ADAPTIVE_FORMULATION.md), and
[reproduction instructions](../../analysis/p2_h_adaptive_donor_20261108/REPRODUCE.md).
Select `P2H.DonorMode=adaptive_multilayer` and `P2H.MaxDonorDepth` in metres;
`top_only` remains the compatibility default. The ALE operator and BDF are
unchanged. Diagnostic `ControlledDonorTest` and `PrecrustOnly` settings are
explicit test modes, not production initial conditions.

User scope: `task/9.md`, isolated branch `research/p2-h-ale-coupling`.
See [equations and limitations](../../docs/P2_H_FORMULATION.md) and
[registered staged plan](PLAN.md). The final Chinese report and accepted
results belong under `analysis/p2_h_ale_coupling_20261107/`.

## Components

- `ale.hh`: exact corrected M3 B2 spatial operator, stable top root and shared
  upwind material mesh flux; standalone constant-flux reservoir regression or
  externally supplied bottom composition and top gas-film callback.
- `bdf.hh`: SUNDIALS 6.7.0 CVODE BDF/Newton integration of conservative pore
  volume E and dissolved NaCl Q, with independent top/internal/exchange ledgers.
- `regression.cc`, `operator_check.cc`, `verify_reference.py`: compiled M3
  regression, reverse-mesh operator comparison, full-history/field comparison
  against the pinned Python archive and independent budgets.
- `soil_methods.inc`, `prepare_dumux.py`: deterministic isolated P1-R2 overlay;
  zero-crust precipitation, conservative pore charging, conditional implicit
  interface exchange, and joint acceptance/rollback.
- `run_soil.py`, `verify_soil.py`: bounded B/C/E matrices, resource refusals,
  injected rejections, separate inventories and actual common-time comparisons.
- `plot_results.py`, `collect.py`: accepted-state figures and lightweight
  evidence, retaining full trajectories and raw fields outside Git.

No artificial salt seed appears in the soil experiment. The explicitly
pre-existing 1 mm crust belongs ONLY to the required standalone M3 regression.
The soil runner reuses frozen P2-F soil inputs for provenance; its inherited
`[P2F]` parameter section and `p2f_1d` filename prefix are unused. No P2-F seed
source is compiled, and tests assert zero initial external solid/domain height.

## HPCC workflow

All builds, simulation, validation, and plotting run under Slurm. The workspace
is `/WORK/hulm_work/zhixin/p2_h_ale_coupling_20261107`.

1. Obtain the official [SUNDIALS v6.7.0 source archive](https://github.com/LLNL/sundials/tree/v6.7.0).
   `prepare_dependencies.py` verifies its pinned SHA and packages only source,
   headers, CMake files and licenses; it does not modify or build upstream code.
2. Export the pinned salt_grow source plus `grid_640.json/.npz` with `git archive`
   at `8c0df6eddeea95a8a76c61f245cef2f891255c86` into `m3_reference.tar`.
3. Stage standalone files under `source/`, the required frozen model directories
   and this module under `repo/model/`, and the two archives at workspace root.
4. `sbatch source/stage_a.slurm` builds the task-local static CVODE installation,
   performs the 80-cell smoke and 640-cell/6-h regression, then applies A gates.
5. After A passes, `sbatch repo/model/p2_h_ale_coupling/stage_b.slurm` builds the
   actual DuMuX overlay and executes the activation/rollback matrix. Each run
   directory must be fresh; the runner refuses to overwrite previous cases.
6. `stage_c.slurm` first verifies B, then runs the bounded conditional coupled
   matrix. Individual physical/numerical failures remain labeled outcomes;
   scheduler completion does not imply that every scientific gate passed.
7. `postprocess.slurm` audits complete accepted histories, generates PNG/SVG
   figures and creates compact and full-history archives with transfer hashes.

Configure P2H `ActivationHeight` (3.5e-6 m default), `CrustCells` (64),
`MinimumSupplySaturation` (.99), `ActivationOnly` (false),
`InjectRejectedTrial` (false), and `RejectAfterActive` (false). The last two
are controlled lifecycle tests. Film settings retain P2G `FilmThickness` and
`RelativeHumidity`. Sources and physical settings are frozen per run archive;
failed attempts are preserved under `failures/`, with separately hashed binaries.

The reservoir experiment driver keeps its historical 0.1-micrometre scenario
when `P2H_ACTIVATION_HEIGHT_M` is absent. Set that environment variable
explicitly for controlled activation-height comparisons. The registered
`x_NaCl=0.08` crystal-scale test uses `P2H_ACTIVATION_HEIGHT_M=3.5e-6` and does
not create a seed; it waits for the actual surface solid inventory to cross the
threshold.

## Interpretation

This is a conditional saturated-crust mass coupling, not a calibrated crust
hydraulic model. Sw>=.99 is a declared supply screen, not a derived pressure
continuity law. Instantaneous local pore charging is itself grid-dependent;
refused activations are outcomes, not evidence that natural salt crusts cannot
form. Subcrystal thresholds do not validate continuum mechanics. M4's material
coordinate distinction is retained: a fixed ALE index is not a material marker.
