# P1 GO/NO-GO decision

## Decision

**NO-GO for starting surface precipitation development.**

The compositional implementation and V1--V7 numerical foundation are a **GO**: the model compiles reproducibly, conserves physical molar H2O/air/NaCl components, respects nonvolatile salt, supports controlled phase appearance/disappearance and passes an independent saved-state budget audit.

The full P1 scientific success criterion additionally requires an instability benchmark under matched assumptions and a one-factor-at-a-time spatial/temporal convergence assessment. The bounded benchmark did not reproduce instability: the seeded mode-28 amplitude decreased monotonically through the registered 60,000 s window. Consequently, an onset time and post-minimum growth rate do not exist in the observed interval, so the onset-sensitive convergence gate cannot be satisfied even if individual terminal observables are insensitive.

## Evidence that passed

1. Both 1D and 2D DuMuX executables were built and executed from an isolated module without modifying shared upstream sources.
2. V1--V7 passed their predeclared gates. The worst solver-side component residual was \(2.66\times10^{-12}\) mol m\(^{-2}\); the worst independent V7 residual was \(1.82\times10^{-12}\) mol m\(^{-2}\).
3. Hydrostatic spurious velocity was \(2.154\times10^{-11}\) m s\(^{-1}\), below the predeclared \(3.066\times10^{-11}\) m s\(^{-1}\) limit.
4. The V4 diffusion error decreased from \(1.78\times10^{-4}\) to \(7.31\times10^{-5}\) under vertical refinement.
5. Controlled phase appearance and disappearance times were unchanged under time-step halving.
6. No precipitation, solid NaCl, mineral source or crust geometry is present.
7. B0 and B1 both completed 264 accepted steps to 60,000 s on the full physical domain with the exploratory 280 × 100 grid. The final top salinity, liquid saturation and horizontally averaged vertical fields agree closely, as expected for a small zero-mean perturbation.

## Evidence that did not pass or remains incomplete

1. B1 mode-28 amplitude decreased monotonically from \(7.61\times10^{-7}\) to \(4.70\times10^{-7}\). Its minimum occurred at the final saved state, so 60,000 s is a right-censored non-detection, not an onset time.
2. The amplitude was still decreasing at the frozen Richards global LSA crossing, mode-28 crossing and numerical-detection reference times. Thus P1 did not reproduce the pre-precipitation instability branch under the matched primary inputs.
3. A post-minimum log-amplitude slope is undefined. This makes the predeclared onset-sensitive convergence check unavailable and therefore failed by construction; it is not imputed from the overall decay slope.
4. The executed benchmark grid was 280 × 100, not the formal 600 × 200 production grid. Its role is bounded mechanism screening, not a production-grid validation.
5. The registered one-factor comparisons use one coarser horizontal mesh, one coarser vertical mesh and a halved time step. They are sensitivity checks around the exploratory base case, not a demonstrated asymptotic convergence order.
6. The constitutive models are not identical: P1 retains upstream `BrineAir` phase equilibrium and transport properties, whereas the frozen Richards solver advances a normalized salt variable with its own density and viscosity closures. The present evidence does not identify which constitutive difference suppresses growth.
7. No experimental comparison validates either instability prediction.

## Registered sensitivity outcome

| Pair relative to 280 × 100, \(\Delta t=227.372\) s | Final top salinity difference | Final mode-amplitude difference | Budget | Post-minimum growth rate |
|---|---:|---:|---:|---:|
| 140 × 100, same \(\Delta t\) | \(1.34\times10^{-11}\) (pass) | 1.53% (pass) | pass | unavailable (fail) |
| 280 × 50, same \(\Delta t\) | 26.96% (fail) | 48.61% (fail) | pass | unavailable (fail) |
| 280 × 100, \(\Delta t/2\) | 0.00773% (pass) | 0.0282% (pass) | pass | unavailable (fail) |

All four trajectories remained monotonically nonincreasing in mode amplitude and placed their minimum at 60,000 s. The maximum absolute component-budget residual among the sensitivity cases was \(4.50\times10^{-10}\) mol m\(^{-2}\), below the registered \(10^{-8}\) mol m\(^{-2}\) limit. Conservation therefore passes while onset-sensitive numerical qualification does not.

## Scientific interpretation

P1 verification establishes that later differences from Richards are not caused by an obvious dimensional budget defect. It does not establish that upstream `BrineAir` properties are the uniquely correct constitutive choice, that the numerical onset is grid-converged, or that post-precipitation instability follows the same branch.

The Richards global LSA crossing, observable mode-28 crossing, Richards numerical detection, P1 amplitude minimum and any nonlinear transport time remain distinct quantities.

## Required next decision

Treat this P1 assessment as complete with a scientific NO-GO. Before P2, define and separately authorize a new P1 revision that isolates the constitutive cause of the missing growth (for example, controlled closure substitutions and, only after a positive growth branch is recovered, the formal 600 × 200 run and a longer registered window). Do not interpret a longer run or finer grid alone as a justified fix without first demonstrating the relevant unstable branch.

P2 precipitation, P3 Prat growth and P4 crust feedback were not implemented or started.
