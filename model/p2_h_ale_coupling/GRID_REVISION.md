# Grid precipitation revision

Branch `research/p2-h-grid-precipitation-fix`, based on
`301097e3d2de3541bf63741cd729bd5cda22329a`.

The correction is explicitly selected with:

```
-P2H.PrecipitationMode implicit_surface
-P2H.DonorMode adaptive_multilayer
-P2H.MaxDonorDepth 0.005
-P2H.GridDiagnostics true
```

`projection` remains available and is the compatibility fallback. The new
`run_grid.py fixed_pre` and `fixed_coupled` suites explicitly select the
corrected mode. Existing benchmark scripts retain their historical behavior.
The activation default remains 3.5e-6 m; epsilon_birth remains 0.43. This is a
one-horizontal-cell demonstration, not a laterally heterogeneous surface solve.

See SURFACE_COMPLEMENTARITY.md for the conservation equations and the physical
assumption retained from the original surface-equilibrium projection.

## Code boundary

- `grid_diagnostics.inc`: scalar implicit surface-flux solve and nonmutating
  finite-volume/field observers, included inside the isolated Problem class.
- `soil_methods.inc`: mode dispatch, accepted surface inventory and unchanged
  adaptive donor transfer; accepted-only commit and rejected-step rollback.
- `bdf.hh`: optional callback captures a trial rejected by the unchanged M3
  physical-bound test. The state, RHS, error weights, BDF settings and physical
  validation conditions are unchanged.
- `ale.hh`, frozen P1-R2 residual/constitutive files and the gas-film source
  are unchanged. No production clipping, new internal-soil reaction or seed.
- `run_ablations.py` generates a **test-only** ALE header overlay replacing
  bottom salt input; it never enters the DuMux build. The equilibrium initial
  reconstruction changes height/initial porosity conservatively and is also
  restricted to standalone tests. Production initialization remains actual
  accumulated solid plus the original conservative adaptive donation.

## Diagnostic semantics

`grid_budgets.csv`: every accepted preactivation step, six physical cut depths;
transport endpoint inventory/flux before projection and accepted inventory
after projection/donation are separate. Activation rows include donation and
are excluded from the ordinary precrust regional projection budget. They are
instead checked against donor inventories and combined budgets.

`coupled_history.csv`: every accepted step, top C/Sw, precipitation stock,
global and interface component inventories. It does not contain rejected
states. `grid_profiles.csv`: full actual BOX geometry and primary/reconstructed
variables at time zero, multiples of 100 s and every accepted step during the
first ten seconds after activation. Near-surface top C/Sw and all cut-region
inventories are recorded every accepted step in the other two files.

`ale_diagnostics.csv`: all accepted active steps, full epsilon/C/R/D_eff/a_v
fields, physical cell bounds, growth and interface fluxes; separately marked
`bound_rejected` for output states rejected by M3's physical-bound check.
This callback does not claim to record CVODE's internal Newton/error-test
rejects. Accepted field values are not extrapolated to a failed final time.

The bottom concentration is constant during each ALE trial. Thus accepted
advective salt integral is exactly C_bottom times the liquid ledger increment;
the remaining salt increment is diffusion. Standalone replays use these step
integrals and split at every original accepted endpoint. Prescribed salt
replays retain M3's liquid-volume continuity closure; they do not prescribe an
independent, potentially incompatible liquid input.

## Execution

`build_grid.slurm`: diagnostic observer build and geometry smoke checks.
`matrix_grid.slurm diagnostic`: historical and controlled precrust grids.
`matrix_grid.slurm ale_diagnostic`: original 160/320 short postactivation fields.
`build_fixed.slurm`: independent corrected build in `correction/`.
`matrix_fixed.slurm fixed_tests|fixed_pre|fixed_coupled`: staged tests, with
coupled execution only after the precrust audit.
`run_ablations.py`, `refine_ablations.py`: conservative causal interventions.
`audit_grid.slurm final`: budgets, independent field checks, comparisons, plots.

Scripts refuse unscheduled execution. The cluster's Slurm MEMORY metadata is
configured as 1 MB per node, so an explicit 8G analysis request was rejected
before job creation; audit scripts use the supported CPU allocation pattern.
Compiler messages and all failed audit/solver attempts are retained.
