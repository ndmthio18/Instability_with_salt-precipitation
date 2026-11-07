# P2-G minimal surface-precipitation model

This is the active simplified development path requested on 2026-11-07.
P2-E/F remain frozen references. No seed transfer, crust reaction solver,
crust pore liquid, hydraulic resistance, or evolving porosity is used here.

- [Equations and assumptions](../../docs/P2_G_FORMULATION.md)
- [Registered bounded benchmark](PLAN.md)
- [HPCC results and decision](../../analysis/p2_g_surface_precipitation_20261107/REPORT_ZH.md)

`prepare.py DESTINATION` copies frozen P1-R2 into an isolated application,
then adds the surface mass projection and independent Prat diagnostic.
It extracts only the gas-film thermodynamics from P2-D. The original files
are never edited. `surface.hh` is the small stand-alone mass/growth module;
`methods.inc` owns the accepted-step inventory and diagnostics.

The supported benchmark is a horizontally uniform column, represented by
one horizontal cell and 40/80/160 vertical cells. The held evaporation flux
and solid inventory are column averages. This implementation is not a
spatially resolved 2D surface model.

`run_matrix.py` derives the soil input from frozen P2-F to retain the same
soil/density/retention settings, but replaces its boundary mode with
`surface_film`. The inherited `[P2F]` section and `p2f_1d` output prefix are
unused historical input text; no P2-F source or seed mechanism is compiled.
Actual overrides are in each run's `command.json`. P2-G parameters are
`EnablePrecipitation`, `ContinuumHeight`, `FilmThickness` and
`RelativeHumidity`, with defaults true, 1e-6 m, .0083 m and .4.

## Reproduction on HPCC

The executed source bundle is recorded by SHA in the result manifest.
Stage the frozen source directories and this module under `source/model/`
in the workspace named in `hpcc.slurm`. The fixed external stack archive
and its required hash are specified in that batch file. Submit with
`sbatch hpcc.slurm`; it builds the isolated application, executes seven
unit checks, runs the bounded seven-case matrix, and analyzes the output.
All build and numerical work runs on allocated compute nodes.

For the independent full-history audit, copy `audit.py` and `audit.slurm`
to that workspace and submit `sbatch audit.slurm` after the first job
completes. The second job does not run the solver. It verifies the exact
growth formula, projection budgets, time continuity, threshold independence,
and paired pre-precipitation histories, and writes grid/time comparisons.
Scientific gate failures are retained in JSON; passing code checks does
not turn a failed scientific gate into a pass.

Full accepted histories and raw fields stay on HPCC. A hash-verified local
archive of the accepted histories is outside Git; Git holds compact sampled
records, all source/output hashes, formulas, and the decision report.
