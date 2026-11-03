# WP8 constrained background implementation

Parent: frozen `model/transport_wp7_localflux_v1/source` at reviewed Paper 2
commit `323fd639c81e47b201ecb0ae8b52ba8449892a59`.

Only `appl/simpleRichardsnc/main.cc` is modified: include the new header, replace
the full Newton solve with the equal-height reduced solve, retain fixed dt.
`wp8_equal_height.hh` is added. Constitutive laws, normalized storage, SCVF
transport, geometry, boundary data and threshold stopping are byte-identical.
The old Newton object remains constructed but is not used for time stepping;
its printed default options are not the constrained Newton tolerances.

The new solver is restricted to one rank and a one-cell horizontal extrusion.
Each height has exactly two identical lifted copies of one pressure/concentration
pair. Residual and Jacobian are summed before solving. Failed Newton steps abort
without being accepted. The dense pivoted reduced solve is limited to this
402-unknown validation prototype; no general parallel production claim is made.

All three isolated compile attempts succeeded. r01/r02 failed before accepting a
step because the common 1e-20 algebraic cutoff was below the measured water
residual rounding floor. r02 added trace diagnostics only. r03 separates water
(2e-18) and salt (2e-20) algebraic cutoffs and retains relative shift <=1e-10.
The independent preregistered 1e-14 uniformity, salt/water closure and bridge
limits were unchanged. r03 passed ten accepted steps and 16 physical/field checks;
T/4T validation is not yet authorized or performed.

Evidence: `analysis/wp8_zero_mode_20261103/build`, `background_smoke10*` and
`BACKGROUND_VALIDATION_PLAN_ZH.md`. Prior header versions are saved by SHA in
`build/source_history`; remote build roots and failed-run archives are preserved.
The final exact source file manifest is `source_manifest.json`.
