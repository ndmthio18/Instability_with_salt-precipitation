# P1-R2 GO/NO-GO decision

## Decision

**GO for accepting P1-R2 as the isolated numerically conservative Richards-salt implementation.**

**NO-GO for starting Prat/crust coupling or making converged physical onset claims from the present evidence.**

This split decision is necessary because discrete component conservation and physical validation are different propositions.

## What is established

1. The existing custom Richards framework now advances physical total-liquid and dissolved-NaCl mass balances in kg and kg s^-1 while retaining pressure and NaCl mole fraction as primary variables.
2. H2O is reconstructed exactly as total minus NaCl within the declared mass-average diffusion closure.
3. V1-V4 pass registered conversion, component-identity, prescribed-flux, variable-density and independent-face checks.
4. V5-V7 reach the requested Case 36877 horizon, close component budgets, show negligible base/half-step differences, and correctly account for a forced rejected/retried step.
5. The single authorized 2D run reaches 136423 s and closes the total and salt budgets to \(9.05\times10^{-14}\) kg and \(8.06\times10^{-16}\) kg.
6. At fixed R24 discretization, the conservative solver preserves the qualitative decay-to-growth mode trajectory. Its late-window projected growth rate is only 0.65% above the frozen normalized reference.
7. No frozen source was modified, no precipitation or Prat model was introduced, and raw fields remain outside Git with archive provenance.

## What is not established

1. The retained Darcy velocity is treated as a mixture mass-average advective velocity. This is an explicit approximation, not exact multicomponent mixture mechanics.
2. Mode A and Mode B produce a 2.30% final mass-fraction-profile difference, far exceeding the \(5.87\times10^{-6}\) time-step sensitivity. The atmospheric boundary semantics therefore remain a material physical-model choice.
3. Stage D contains one R24 x 200 run. It does not provide a new horizontal/vertical refinement sequence or a continuum-extrapolated growth rate.
4. The profile-error threshold for V5 was not separately pre-registered; those comparisons are quantitative diagnostics rather than a retroactively declared convergence PASS.
5. No experiment validates the mixture closure, evaporation boundary mode, onset time, surface salinity or downward salt-flux profile.
6. Numerical conservation does not by itself validate a future crust-growth law, precipitation kinetics or a no-lateral-transport interface approximation.

## Evidence-conclusion alignment

The strongest defensible conclusion is:

> P1-R2 supplies a reproducible, discretely conservative total-liquid/NaCl solver and numerically audited component fluxes on the tested grids. The solver is suitable as the numerical foundation for a later, separately authorized Prat-interface study, but the present work does not validate that coupled physics or its production-grid predictions.

Calling the current model "physically validated" would be an evidence-conclusion mismatch. Calling it "numerically conservative under the declared closure" is supported.

## Gates required before P2

Before adding a crust sink, source or geometry update, a separately authorized task should:

1. justify or calibrate Mode A versus Mode B from the intended atmospheric control variable;
2. pre-register horizontal, vertical and temporal refinement tolerances for surface salinity, component flux, onset/minimum time and late growth rate;
3. demonstrate near-surface vertical resolution and mode-growth stability on at least one finer grid;
4. specify the exact top-patch sign convention and time discretization for soil/crust exchange;
5. preserve the identities \(F_{s,soil}=dM_{s,crust}/dt\) and \(F_{w,soil}=j_{atm}+dM_{w,crust}/dt\) without an arbitrary salt sink;
6. separate numerical closure tests from calibration and experimental validation.

Until these gates are met, P2/P3/P4 remain held. No additional production simulation is justified by P1-R2 alone.
