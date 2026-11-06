# P1-R2 conservation verification

## Registered gates

`model/p1_r2_conservative_richards/verification/acceptance_criteria.json` was written before the solver runs. Controlled tests use an absolute budget tolerance of \(10^{-10}\) kg plus \(10^{-8}\) of the larger of initial inventory or integrated absolute boundary transfer. Case 36877 uses a relative budget tolerance of \(10^{-6}\). Conversion identities use \(5\times10^{-13}\), the analytical/finite-difference diffusion derivative uses \(10^{-9}\), and prescribed/independently reconstructed fluxes use \(10^{-8}\) relative tolerance.

These gates test discrete numerical closure. They do not validate the retained constitutive physics or the mixture-velocity approximation.

## Unit conversion and V4

The corrected audit reads the exact frozen Case 36877 input rather than hardcoding approximate constants. At \(x_s=0.0036\), it obtains

| Quantity | Corrected value |
|---|---:|
| \(\bar M\) | 0.018160808992 kg mol^-1 |
| \(w_s\) | 0.0115845059596561 |
| \(\rho_l\) | 1006.02410443354 kg m^-3 |
| \(c_l\) | 55395.3353551984 mol m^-3 |

Regression tests reject the rounded substitutions used by the earlier audit. The component identities \(m_i=M_in_i\) and \(m_w+m_s=m_{total}\) pass. The analytical derivative

\[
dw_s/dx_s=M_sM_w/\bar M^2
\]

agrees with centered finite differences to better than the registered \(10^{-9}\) relative tolerance over \(x_s=0\) to 0.1.

## V1-V3 and V6 solver tests

| Test | Accepted steps | Final time (s) | max total residual (kg) | max salt residual (kg) | max water residual (kg) | independent face error (kg s^-1) | Result |
|---|---:|---:|---:|---:|---:|---:|---|
| V1 closed box | 10 | 1000 | 3.61e-16 | 2.08e-17 | 3.59e-16 | 0 | PASS |
| V2 prescribed components | 10 | 1000 | 5.38e-16 | 2.00e-17 | 5.18e-16 | 0 | PASS |
| V3 variable density | 20 | 2000 | 6.94e-16 | 2.60e-16 | 4.44e-16 | 0 | PASS |

For V2, the expected and observed total flux are both \(6.0\times10^{-7}\) kg s^-1, and the expected and observed salt flux are both \(6.0\times10^{-9}\) kg s^-1 at recorded precision.

V3 is paired with a two-control-volume algebraic discriminator. Redistribution can preserve the old normalized quantities \(\sum\theta\) and \(\sum\theta x_s\) while changing both \(\sum\theta\rho\) and \(\sum\theta\rho w_s\) by more than \(10^{-4}\) relative. This proves non-equivalence of the residuals; it is not presented as an experimental validation of the new closure.

## V5-V7

All three full 1D Case 36877 conservative runs reached 136423.0637723484 s.

| Run | max total residual (kg) | max salt residual (kg) | max water residual (kg) |
|---|---:|---:|---:|
| Mode A, base \(\Delta t\) | 1.38e-13 | 1.60e-15 | 1.36e-13 |
| Mode B, base \(\Delta t\) | 1.15e-13 | 1.33e-15 | 1.13e-13 |
| Mode A, \(\Delta t/2\) | 5.60e-14 | 6.79e-16 | 5.54e-14 |

The controlled retry test requested 100000 s, rejected that step, and accepted 50000 s after one retry. Its second 50000 s step completed the 100000 s horizon. Budgets are evaluated only at accepted states, so the rejected nonlinear iterate is not counted as physical evolution.

The bounded 2D test completed 600 accepted steps with no retries. Its maximum total, salt and reconstructed-water residuals were \(9.05\times10^{-14}\), \(8.06\times10^{-16}\), and \(8.98\times10^{-14}\) kg, respectively.

## Evidence locations

- Unit tests: `analysis/p1_r1_unit_audit_20261106/test_conversion_audit.py` and `model/p1_r2_conservative_richards/verification/test_physics.py`.
- V1-V3 summary: `model/p1_r2_conservative_richards/results/stage_b_20261106/gate_summary.json`.
- V5-V7 accepted-step and component budgets: `model/p1_r2_conservative_richards/results/stage_c_20261106/runs/`.
- 2D accepted-step and component budgets: `model/p1_r2_conservative_richards/results/stage_d_20261106/run/`.

Raw field archives are intentionally excluded from Git; their remote paths and SHA-256 values are recorded in the provenance manifest.
