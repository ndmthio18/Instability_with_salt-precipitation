# Failure and correction register

## Verification job 16430327

The first short verification wrapper exited before numerical integration. Its
input builder replaced both the problem name and
`EvaporationBoundaryMode=p2f_analytical` with the diagnostic label
`p2h_ale_reservoir`, which is not a valid DuMuX boundary mode. All three cases
therefore reported `Unknown Problem.EvaporationBoundaryMode`.

The failed directory is preserved on HPCC as
`verification_failed_16430327/`. No accepted state or scientific result from
that attempt is used. The builder was corrected to change the problem name and
boundary mode separately. Job `16430417` then completed the hydrostatic,
reservoir-flux, and frozen-no-flux regression tests.

## Evidence classification

- A Slurm `COMPLETED` state is only execution evidence.
- A run reaching 50 h without accepted ALE activation is not a successful
  coupled sensitivity calculation.
- Early minimum-porosity or solver-protection stops remain censored and are not
  extrapolated.
- Global conservation does not by itself validate local ALE clogging or the
  continuum interpretation of a `0.1 µm` activation layer.

## Pilot gate 16430419

The solver completed the requested 180000 s, but the scientific verifier
returned 1 because no surface precipitation or ALE activation occurred. Final
surface NaCl mass fraction was 0.262856472 versus equilibrium 0.264940870.
This is the registered `FAIL_NO_SURFACE_PRECIPITATION_NO_ALE_ACTIVATION`, not a
solver crash. Dependent ALE-pilot, matrix, and standard post-processing jobs
were cancelled without running a solver.

## Diagnostic post-processing

Job `16433911` failed before plotting because it searched for the frozen
no_flux compact history inside the copied source repository rather than the
frozen HPCC analysis root. The failed log is retained. The source path was
corrected without changing numerical results; job `16433914` completed.
Artificial scientific data were not generated to replace the missing matrix.
Manual visual QA then identified Matplotlib offset notation on near-unity `Sw`
axes; job `16435674` regenerated the same data with plain tick labels.
