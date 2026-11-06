# P1 numerical verification

## Predeclared gates

Acceptance criteria were stored and hashed before the accepted run. The absolute component-budget limit was \(10^{-9}\) mol m\(^{-2}\); nonzero signals additionally used a relative check. V3 required pressure drift below 0.1 Pa and spurious Darcy velocity below one percent of the Case 36877 evaporation velocity, \(3.0655444\times10^{-11}\) m s\(^{-1}\). V4 required less than two percent coarse-grid analytical error and an error-reduction ratio of at least 1.5. V5 required a maximum 1D/2D relative difference below \(10^{-8}\). V6 required a phase-state change and transition-time disagreement below five percent under time-step halving.

## Accepted run

The accepted suite is `/temp-hd/zhixin/p1_compositional_20261106/verification_run_06`; the compact committed record is `model/p1_compositional_baseline/results/p1_verification_20261106/gate_summary.json`.

| Gate | Result | Evidence |
|---|---:|---|
| V1 closed box | PASS | maximum accepted-step residual \(1.95\times10^{-12}\) mol m\(^{-2}\) |
| V2 prescribed H2O/air/NaCl flux | PASS | worst residual \(2.66\times10^{-12}\) mol m\(^{-2}\) |
| V3 two-phase hydrostatic equilibrium | PASS | pressure drift 0.00801 Pa; maximum phase Darcy velocity \(2.154\times10^{-11}\) m s\(^{-1}\) |
| V4 passive diffusion | PASS | relative amplitude error \(1.78\times10^{-4}\) at 40 cells and \(7.31\times10^{-5}\) at 80 cells; error ratio 2.44 |
| V5 1D/2D consistency | PASS | maximum relative difference \(9.17\times10^{-9}\) |
| V6 phase switching | PASS | appearance at 110 s for dt 10/5 s; disappearance at 670 s for dt 10/5 s; physical bounds retained |
| V7 independent budget | PASS | worst independently reconstructed residual \(1.82\times10^{-12}\) mol m\(^{-2}\) |

V6 uses controlled air-component injection into a liquid-only state and controlled air removal from a two-phase state. These are verification tests, not Case 36877 physics claims.

## Failed attempts retained

- Build attempt 01 lacked the module pkg-config template.
- Build attempt 02 exposed the SuperLU header mismatch and a missing local-residual constructor argument.
- Verification run 02 exposed inactive-phase diagnostic handling, a non-hydrostatic gas initial pressure and an ill-posed V5 phase label.
- Verification runs 04 and 05 failed the unchanged V3/V4/V6 gates and motivated physically traceable corrections; their failures were not reclassified as passes.

Passing V1--V7 demonstrates numerical accounting and controlled state transitions. It does not validate density-instability onset, the `BrineAir` constitutive choice against experiment, or any precipitation mechanism.
