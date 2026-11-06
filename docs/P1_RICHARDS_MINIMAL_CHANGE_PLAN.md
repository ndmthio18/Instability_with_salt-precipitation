# P1-R1 Minimal Change and Verification Plan

## Decision

Adopt a **total-liquid-mass plus salt-mass residual** in an isolated copy of the existing custom Richards solver, while retaining pressure and NaCl mole fraction as primary variables.

Do not:

- continue the rejected 2pnc route;
- create a separate RichardsNC framework;
- edit the frozen Paper-1 snapshot;
- change permeability, porosity, retention/relative-permeability laws, density law, viscosity, gravity, mesh, time stepping, initial perturbation or Case 36877 demand to improve agreement;
- run production cases before the unit tests close;
- merge or push without explicit authorization.

No solver prototype was created in P1-R1 because the diffusion reference system and the physical meaning of the evaporation demand remain decision gates. Guessing either would change the scientific model.

## Minimal implementation boundary

After the two gates are approved, create a new, versioned source tree rather than modifying `model/transport_wp8_background_v2`.

The smallest defensible code change is:

1. Keep the current geometry, primary variables, Darcy velocity, capillary/relative-permeability closures, density function and initial/boundary topology.
2. Add one audited conversion helper for \(x_s\leftrightarrow w_s\), \(\bar M\), \(\rho\), and \(c\).
3. Replace only the custom storage kernel with \(\theta\rho\) and \(\theta\rho w_s\).
4. Replace only the custom component flux kernel with \(\rho\mathbf q\) and \(\rho w_s\mathbf q-\rho D_{\mathrm{eff}}\nabla w_s\).
5. Convert the top boundary to explicitly declared mass units; retain zero salt escape before precipitation.
6. Rewrite diagnostics to report independent total-mass, salt-mass and derived water-mass inventories and boundary fluxes.
7. Record all source hashes and parameters before any run.

This reuses the existing framework and the scientifically relevant upstream conventions without importing the full upstream RichardsNC residual, which would require a broader fluid-system and volume-variable redesign.

## Required pre-implementation gates

### G1 — conservation basis

Approve formulation C: equation 0 is total liquid mass; equation 1 is dissolved salt mass. Water mass is equation 0 minus equation 1.

### G2 — diffusion reference system

Approve mass-averaged Fick diffusion in mass fraction:

\[
\mathbf j_s^m=-\rho D_{\mathrm{eff}}\nabla w_s.
\]

If the research intent is instead to preserve \(-D_{\mathrm{eff}}\nabla x_s\), define its physical coefficient and units as a new empirical closure and repeat the derivation before implementation.

### G3 — evaporation basis

Document whether `EvaporationRate` is:

- an equivalent pure-water volume flux with a specified reference density;
- an H2O mass flux;
- an H2O molar flux; or
- an atmospheric demand requiring a constitutive surface-limitation law.

No silent conversion is allowed.

### G4 — registered numerical tolerances

Record nonlinear/linear solver tolerances and the conservation thresholds below before observing test results. Tighten the proposed thresholds only for demonstrated floating-point or quadrature reasons; do not loosen them to accept a failed implementation.

## Verification ladder

Define the normalized closure error for component/balance \(i\) over \([t_0,t_1]\):

\[
\epsilon_i=
\frac{|M_i(t_1)-M_i(t_0)+\int_{t_0}^{t_1}F_{i,\partial\Omega}\,dt|}
{\max(|M_i(t_0)|,\int|F_{i,\partial\Omega}|dt,M_{i,\mathrm{floor}})}.
\]

For controlled verification cases, the proposed pass threshold is \(\epsilon_i\le10^{-8}\); for the full Case 36877 pilot, \(\epsilon_i\le10^{-6}\). Absolute floors must be registered in the input manifest to avoid declaring a near-zero component accurate through an unstable relative denominator.

### V1 — closed-domain conservation

- Close every boundary.
- Use nonuniform \(S_w\) and \(x_s\) so storage is nontrivial.
- Reconstruct total, salt and water mass independently from saved cell states.
- Pass: each component closes to \(10^{-8}\), and water mass equals total minus salt within the same absolute budget.

### V2 — constant component-flux test

- Prescribe a known constant total/water flux and, separately, a known salt component flux.
- Pass: inventory slope matches the imposed integrated boundary flux in magnitude and sign to \(10^{-8}\).
- Include a salt-impermeable evaporation case to prove that water can leave while salt cannot.

### V3 — variable-density discriminator

- Construct a transient composition field for which \(\rho(x_s)\) and \(c(x_s)\) vary measurably.
- Evaluate both the old diagnostic-only conversion and the new physical residual.
- Pass: the physical residual closes to \(10^{-8}\); the test must be strong enough that the old diagnostic mismatch exceeds ten times the numerical closure floor. This prevents an accidentally constant-property test from appearing diagnostic.

### V4 — mass/molar conversion identities

- For every cell and test state, verify \(m_i=M_in_i\), \(m_w+m_s=\theta\rho\), and \(w_s=x_sM_s/\bar M\).
- Pass: maximum relative identity error \(\le5\times10^{-13}\), subject to an explicit absolute floor near zero.

### V5 — zero-precipitation Case 36877 comparison

- Run the frozen normalized model and the conservative prototype over identical accepted physical times.
- Compare \(x_s\), \(S_w\), \(\rho\), total/salt inventories, perturbation amplitude and fitted growth rate.
- Do not tune material properties to recover the old curve.
- Pass as a **comparison package** if both runs are numerically accepted and the conservative run closes to \(10^{-6}\). Report profile \(L_2/L_\infty\) differences and growth-rate shift. A growth-rate shift above 5% triggers reassessment of prior instability claims; it is not grounds to reject the physical balance law.

### V6 — independent boundary-flux reconstruction

- Sum physical fluxes directly over boundary sub-control-volume faces using saved states and face geometry.
- Reconstruct the bottom Dirichlet reaction independently rather than defining it as the remainder required to close the budget.
- Pass: face-sum and residual-derived boundary totals agree to \(10^{-8}\) in unit tests and \(10^{-6}\) in the Case 36877 pilot.

### V7 — accepted-step and time-refinement consistency

- For every accepted step, compare the backward-Euler storage increment with the time-discrete face flux used by the residual.
- Repeat V2/V3 with \(\Delta t\), \(\Delta t/2\), and \(\Delta t/4\).
- Pass: accepted-step closure meets the relevant threshold; independently quadrature-based budget error decreases consistently under refinement and the finest result meets the registered tolerance.

## Evidence package for any future prototype

The prototype is not eligible for a production run until it contains:

- a source and parameter hash manifest;
- exact equation ordering and sign convention;
- declared mass/molar and diffusion reference systems;
- V1–V4 machine-readable results and logs;
- independent inventory and face-flux reconstruction scripts;
- a comparison manifest pairing runs by actual common accepted time;
- a statement separating numerical closure from physical validation.

Passing V1–V7 establishes implementation consistency and conservation. It does **not** validate the density law, diffusion coefficient, evaporation closure, or precipitation physics against experiment.

## Answers to the six P1-R1 questions

1. **Can we use the original Richards–salt solver directly?** Yes as the retained numerical architecture and as a normalized baseline; no as an already physical component-conservation solver.
2. **Can its units be converted without modifying the governing equations?** No. A reference-state post-processing conversion is approximate and fails exact closure when density or mean molar mass varies.
3. **What is the smallest necessary source modification?** In an isolated copy, replace only storage, component flux, top-boundary conversion and diagnostics; retain pressure/mole-fraction primary variables and all hydraulic/material parameter functions.
4. **Should we use molar conservation or mass conservation?** Use total liquid mass plus NaCl mass with a declared mass-average diffusion frame. It connects directly to crust mass storage and avoids an unnecessary molar-reference ambiguity.
5. **Can the resulting model provide physical water and NaCl fluxes to Prat?** Yes, conditionally: after V1–V7, total and salt mass flux yield water by subtraction and enter the two interface balances without an approximate unit factor. This is numerical consistency, not experimental validation.
6. **What must pass before P2?** G1–G4 and V1–V7, including independent saved-state inventories, independent boundary-face fluxes, variable-density discrimination, accepted-step closure, time refinement and the zero-precipitation Case 36877 comparison. P2 remains blocked until these pass and the user separately authorizes it.
