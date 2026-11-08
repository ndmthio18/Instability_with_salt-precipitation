# Failure and correction register

## Build 16407531

The first build stopped after 6 s because the source transfer omitted the
read-only `model/p2_d_coupling/cpp/pratcrust.hh` dependency used by
`prepare_dumux.py`. The missing dependency was added without changing equations
or numerical settings. Build 16407548 then completed and passed the `a_v`
scaling/top-law-independence unit test.

## Pilot 16407572

The two concurrent pilot processes initially used the same temporary manifest
name. The 128-cell process failed its atomic rename; the independent 64-cell
calculation remained valid and stopped at its registered model screen. Temporary
manifests were changed to include the process ID. The 128-cell pilot was rerun as
job 16407589. No accepted numerical state was overwritten.

## Wrong initial-concentration namespace: job 16407598

The first five-case matrix passed `Transport.SaltMoleFraction=0.05`, but the
actual initial-condition parameter belongs to the `Problem` group. The accepted
soil field therefore remained at the historical near-equilibrium value
(`x_NaCl=0.09999`). This was detected during mandatory visual QA because the soil
mass-fraction colour scale remained near 0.26494.

These runs are retained under `production/` as failed provenance and are excluded
from every final scientific table and figure. The corrected matrix uses
`Problem.SaltMoleFraction=0.05` under `production_corrected/`. Verification was
strengthened to reconstruct the initial mole fraction from the accepted mass
fraction and require agreement with 0.05.

The earlier 64/128 ALE pilots (jobs 16407572 and 16407589) predate this
correction and therefore also used the historical near-equilibrium accepted
field.  They remain only a bounded ALE cost/early-scalar check and are not
evidence for the requested `x_NaCl=0.05` production initial condition.

## Scientific stopping records

The final report and `results/summary.csv` retain every early model stop. A Slurm
wrapper exit of zero only means that the solver exit code and accepted prefix were
recorded; it does not convert a protected model stop into a completed 50 h run.

## Post-processing job 16411448

The first final post-processing attempt was deliberately cancelled after visual
delivery QA found that vectorizing every `pcolormesh` cell produced a 219 MB SVG.
No numerical calculation was cancelled or changed.  The cloud-map artists were
marked `rasterized=True`; job 16413799 then passed all 92 checks.  A second visual
QA found that empty no-activation panels inherited Matplotlib's default 0--1
vertical range, which could be misread as a physical 1 um crust.  The axes were
fixed explicitly to 0--0.1 um and final job 16413869 passed all 92 checks.

## Corrected production scientific result

All five corrected cases reached 180000 s but never accepted an ALE activation.
The first surface precipitation was at 108399 s and the first activation attempt
at 108996 s.  Within 5 mm, every donor had `Sw < 0.99` (surface `Sw=0.864725`;
largest recorded value about 0.877492), so admissible liquid was exactly zero.
This is retained as `FAIL_NO_ALE_ACTIVATION`; reaching the requested horizon is
not classified as successful coupled sensitivity analysis.
