# Research Task: Conservative Two-Moving-Boundary ALE Validation of the Corrected Porous Salt Crust Growth Law

## 1. Repository and baseline

Repository:

`ndmthio18/Instability_with_salt-precipitation`

Remote:

`instability`

Current branch:

`research/p2-h-corrected-law-longterm`

Baseline commit:

`c7da956041b7e42f85278603a068c3915394ec20`

Create a new branch:

`research/p2-h-double-moving-boundary-validation`

Review the previous studies before implementation:

- `analysis/p2_h_corrected_vs_ale_x008_20261112/`
- `analysis/p2_h_corrected_law_validation/`
- `analysis/p2_h_corrected_law_longterm/`
- `model/p2_h_ale_coupling/ale.hh`
- `model/p2_h_ale_coupling/bdf.hh`
- `docs/P2_H_FORMULATION.md`

Preserve all existing validated models and results.

Do not modify the original ALE implementation in ways that break previous results. Implement the double-moving-boundary formulation as a separate model or an explicitly enabled alternative.

## 2. Main scientific objective

We proposed a corrected porous salt growth law incorporating downward diffusion of NaCl into the underlying soil.

Previous results established that:

1. The corrected finite-Da analytical top-growth law agrees with the short-time ALE solution under matched boundary conditions.
2. At large Da, downward diffusion has a weak effect on top precipitation because of reaction–diffusion screening.
3. Nevertheless, significant dissolution develops near the bottom of the crust.
4. With u > 0, previous long-term simulations stopped when the local solid fraction became nearly zero.
5. The current fixed-bottom ALE cannot independently distinguish top growth, bottom recession, and net crust thickening.

The objective is therefore:

**Develop a conservative two-moving-boundary ALE model and use it to independently evaluate the corrected analytical predictions for top growth, bottom recession, and net crust thickness.**

This task is primarily a validation study, not an unrestricted effort to develop a new salt-crust model.

The two-moving-boundary ALE must not use the corrected analytical growth law to determine either moving boundary.

The analytical law must remain a prediction that can be falsified.

## 3. Clearly distinguish three observables

Define:

\[
z_t(t)=\text{top position}
\]

\[
z_b(t)=\text{bottom position of the solid-bearing crust}
\]

\[
h_c(t)=z_t(t)-z_b(t)
\]

The geometric net-thickening rate is:

\[
\boxed{
v_{\mathrm{net}}=\frac{dh_c}{dt}=v_t-v_b
}
\]

where:

\[
v_t=\frac{dz_t}{dt},\qquad v_b=\frac{dz_b}{dt}.
\]

Also track total crystalline salt mass:

\[
M_{\mathrm{cr}}(t)
=
\int_{z_b}^{z_t}
\rho_{\mathrm{cr}}[1-\epsilon(z,t)]\,dz.
\]

Do not confuse:

- Actual top position.
- Actual bottom position.
- Geometric crust thickness.
- Total crystal mass.
- Inventory-equivalent thickness.

These are distinct observables.

## 4. Analytical predictions to validate

### 4.1 Corrected finite-Da top-growth law

Use the previously validated formula:

\[
\boxed{
v_t^{\mathrm{corr}}=
\frac{C_{\mathrm{eq}}j}
{\rho_{\mathrm{cr}}(1-\epsilon_t)(1-C_{\mathrm{eq}})}
\frac{
X[1-u\,\mathrm{sech}(\sqrt{Da})]
}{
X+\tanh(\sqrt{Da})
}
}
\]

where:

\[
X=\sqrt{Da^*},
\qquad
u=\frac{B_{\mathrm{soil}}}{C_{\mathrm{eq}}j}.
\]

Retain the original definitions of Da and Da* and check that the correct active-crust thickness is used.

### 4.2 Corrected net-growth relation

The previous study also proposed a crystal-inventory-based net-growth relation:

\[
\boxed{
\frac{dh_{\mathrm{net}}^{\mathrm{corr}}}{dt}
=
\frac{P(t)}
{\rho_{\mathrm{cr}}\phi_{\mathrm{core}}}
}
\]

with:

\[
P(t)=
\frac{C_{\mathrm{eq}}j(t)-B_{\mathrm{soil}}(t)}
{1-C_{\mathrm{eq}}}
\]

and:

\[
\phi_{\mathrm{core}}
=
(1-\epsilon_t)
\frac{1+\sqrt{Da^*}}{\sqrt{Da^*}}.
\]

CRITICAL:

This relation was previously validated only as a crystal-inventory-equivalent height.

It must not automatically be treated as a prediction of actual geometric net growth.

First investigate the additional assumptions needed for:

\[
h_{\mathrm{net}}^{\mathrm{corr}}
\approx z_t-z_b.
\]

In particular, analyze the role of:

- Spatially nonuniform solid fraction.
- Internal precipitation and dissolution.
- Evolving crust structure.
- Stable core solid fraction.
- Crust translation and material deformation.
- Dissolved-salt storage.

Do not assume that the geometric net-growth rate equals P divided by a theoretical solid fraction unless the conditions required for this interpretation are demonstrated.

This distinction is central to the validation.

## 5. Stage A — Derive the moving-boundary conservation equations

Before coding, derive the equations for a system containing:

1. A solid-bearing porous NaCl crust.
2. A liquid region beneath the crust.
3. A moving top precipitation boundary.
4. A moving bottom dissolution boundary.

Initially assume isothermal conditions and continuous liquid contact.

Do not introduce air gaps, mechanical collapse, or unsaturated flow in this first implementation.

### 5.1 Moving coordinate

For the porous crust:

\[
z=z_b(t)+\xi h_c(t),
\qquad 0\le\xi\le1.
\]

The grid velocity is:

\[
\boxed{
w(\xi,t)=v_b+\xi(v_t-v_b)
}
\]

Use this velocity consistently in the ALE conservative transport equations.

### 5.2 Internal salt transport

Retain the existing ALE physics:

- Dissolved-salt storage.
- Advective transport.
- Molecular/effective diffusion.
- Internal precipitation.
- Internal dissolution.
- Porosity evolution.
- Moving-mesh transport.
- Top surface precipitation.

### 5.3 Bottom interface conservation

Define the bottom dissolution flux \(D_b\) as a positive mass flux of crystalline NaCl dissolved at the moving crust–liquid interface.

If the solid framework is stationary in the laboratory frame, the sharp-front solid mass jump gives:

\[
\boxed{
D_b=\rho_{\mathrm{cr}}\phi_bv_b
}
\]

where \(\phi_b\) is the solid fraction of the intact material being consumed immediately ahead of the recession front.

This formula is valid only for the explicitly stated stationary-skeleton, sharp-interface interpretation.

If the solid skeleton moves, include the skeleton velocity in the interface kinematics.

Do not estimate \(\phi_b\) using a nearly depleted cell where \(\phi_b\rightarrow0\), as this would produce a singular or physically misleading recession velocity.

The dissolved-salt jump condition must use fluxes relative to the moving interface.

For upward-positive fluxes, derive the Rankine–Hugoniot relation consistently for both regions, including the interfacial dissolution source.

Also derive the corresponding liquid/water and total-volume balances.

Do not assume that the liquid mass flux must be identical on both sides of a moving, phase-changing interface without accounting for interfacial storage and displacement.

### 5.4 Overall solid conservation

Verify:

\[
\boxed{
\frac{dM_{\mathrm{cr}}}{dt}
=
G_{\mathrm{top}}
-D_b
+
\int_{z_b}^{z_t}R(z,t)\,dz
}
\]

for the selected sharp-interface representation.

Here:

- \(G_{\mathrm{top}}\) is solid formation at the top boundary.
- \(D_b\) is dissolution at the receding bottom boundary.
- \(R\) is internal precipitation/dissolution in the remaining crust.

Do not double-count bottom dissolution as both a volume reaction and an interfacial reaction.

## 6. Stage B — Determine a physically defensible bottom-motion closure

This is the most important design requirement.

The conservation equations alone do not necessarily determine the bottom recession velocity uniquely.

Before implementation, establish how the model independently determines \(D_b\) and \(v_b\).

Evaluate the following candidate approaches.

### Option 1: Crystal-depletion-front tracking

Retain the finite-rate internal precipitation/dissolution model.

As the bottom solid-bearing material becomes depleted, identify the connected depleted region starting from the soil-facing boundary.

Convert fully depleted material into the adjacent liquid region through a conservative front-tracking or finite-volume remapping procedure.

The moving bottom boundary represents the boundary of the remaining solid-bearing material.

The front position should be determined by actual solid mass depletion, not by the corrected net-growth law.

Use a conservative subcell or cut-cell treatment where necessary to reduce artificial grid dependence.

Do not treat an isolated liquid-only cell inside the crust as automatic evidence that the entire bottom boundary has moved.

### Option 2: Explicit kinetic sharp-interface dissolution

If Option 1 is not sufficiently well-defined, formulate a bottom reaction-front law using local undersaturation and a physically justified reactive surface measure.

The bottom dissolution flux must depend on local interface conditions and available crystalline material.

State all additional assumptions and parameters explicitly.

Do not introduce an arbitrary bottom recession coefficient solely to force agreement with the corrected growth law.

### Required decision

Compare these approaches theoretically before choosing an implementation.

Select the simplest approach that is mass-conservative, physically interpretable, and numerically testable.

If neither approach can uniquely determine bottom motion without additional physical assumptions, document the underdetermination explicitly.

Do not conceal it by prescribing a fitted \(v_b(t)\).

## 7. Stage C — Two-region liquid/crust transport

If the bottom solid-bearing front moves upward, a liquid-only region may remain beneath the crust.

This region must retain its water and dissolved NaCl inventory.

Do not simply delete depleted ALE cells.

Use:

\[
\text{soil boundary}
\rightarrow
\text{liquid region}
\rightarrow
\text{porous salt crust}
\rightarrow
\text{atmosphere}.
\]

The liquid region has no crystalline salt skeleton, whereas the porous crust contains both liquid and solid NaCl.

Solve conservative dissolved-salt transport in both regions.

At their moving interface, enforce:

- Appropriate concentration/contact conditions.
- Relative dissolved-salt flux balance.
- Crystal dissolution and solid mass balance.
- Compatible liquid and volume conservation.

For the first version, prescribed bottom forcing is permitted.

Do not claim independent hydraulic supply prediction unless pressure/flow equations are actually solved.

Make clear that the liquid region is a conditional model assumption, not evidence of an experimentally observed liquid gap.

## 8. Stage D — Verification tests before production calculations

Implement and pass the following tests.

### Test 1: No bottom dissolution

Set conditions with \(D_b=0\).

Require:

\[
z_b(t)=z_b(0).
\]

The new model must recover the original validated fixed-bottom ALE solution to numerical tolerance.

### Test 2: Bottom dissolution only

Disable top precipitation for this test and prescribe a well-defined bottom dissolution problem.

Verify bottom recession and total crystalline mass loss independently.

The recession solution must be obtained from the physical interface balance, not prescribed from the desired outcome.

### Test 3: Top growth and bottom recession

Allow both processes.

Verify:

\[
\frac{dh_c}{dt}=v_t-v_b.
\]

Also verify global and interface-local water and salt budgets.

### Test 4: Zero-solid limit

Test conservative conversion of depleted porous cells into liquid-only cells.

Check that:

- Crystal mass remains nonnegative.
- Dissolved salt is not lost.
- Pore liquid is not duplicated or deleted.
- Material transfer does not introduce discontinuous mass jumps.
- The result converges with mesh refinement.

### Test 5: No artificial validation

Introduce a controlled perturbation to the predicted bottom recession or kinetic parameter.

Confirm that the corrected analytical growth-law comparison responds independently.

This test is intended to detect circular validation or a hard-coded analytical answer.

## 9. Stage E — Controlled u > 0 simulations

After the verification tests pass, run:

\[
u=0,\ 0.3,\ 0.6,\ 0.9.
\]

Use the previous baseline parameters:

- \(h_0=3.5\,\mu\mathrm m\).
- Initial porosity 0.43.
- \(C_{\mathrm{eq}}=0.264940869768\).
- \(k_r=10^{-3}\,\mathrm{m/s}\).
- \(\beta=2\).
- \(D_s=1.3\times10^{-9}\,\mathrm{m^2/s}\).
- Constant initial reference evaporation \(j=1.78\times10^{-5}\,\mathrm{kg\,m^{-2}\,s^{-1}}\).

The target physical duration is 50 h.

Unlike the previous model, do not automatically stop when a bottom cell loses all its crystal material if the new conservative moving-front formulation can validly transfer that region into the liquid domain.

However, continue only while:

- The governing equations remain applicable.
- The interface conditions remain valid.
- The remaining solid-bearing crust is physically admissible.
- Concentrations and inventories satisfy bounds.
- The numerical solution remains convergent.

If a simulation cannot validly reach 50 h, record the actual physical or numerical limitation.

Do not bypass conservation failures or invent mechanical support to force a complete trajectory.

Then repeat the study using the available time-dependent x0=0.08 evaporation and downward-diffusion histories.

Do not mix this one-way forcing with a fully coupled DuMuX calculation.

## 10. Stage F — Independent validation of our corrected growth law

The central comparison is:

\[
\boxed{
v_{\mathrm{net}}^{\mathrm{ALE}}
=v_t^{\mathrm{ALE}}-v_b^{\mathrm{ALE}}
}
\]

against the corrected analytical net-growth prediction, when its assumptions allow a geometric interpretation.

Also independently compare:

\[
v_t^{\mathrm{ALE}}
\quad\text{and}\quad
v_t^{\mathrm{corr}}.
\]

Integrate the ALE boundary velocities to obtain:

\[
z_t(t),\qquad z_b(t),\qquad h_c(t).
\]

Do not convert ALE crystalline mass into height using the analytical \(\phi_{\mathrm{core}}\) as the primary comparison.

The analytical model must not be used to generate the ALE recession speed.

### Separate two hypotheses

**H1 — Corrected top-growth law**

The finite-Da corrected law predicts the actual top precipitation velocity under its stated transport assumptions.

**H2 — Corrected net-growth law**

The proposed crystal-inventory net-growth relation also predicts actual geometric crust thickening under additional assumptions about the internal solid fraction and geometry.

Evaluate H1 and H2 separately.

It is possible for H1 to be valid while H2 fails.

Such a result would be scientifically important.

### Additional mass diagnostics

Compute:

\[
M_{\mathrm{cr}}(t)
\]

\[
M_{\mathrm{dissolved}}(t)
\]

\[
G_{\mathrm{top}}(t)
\]

\[
D_b(t)
\]

\[
\int R\,dz
\]

and compare the exact conservation budget to the approximations used in the corrected analytical model.

Test whether the predicted \(\phi_{\mathrm{core}}\) actually corresponds to a stable spatial region of the computed crust.

If no stable core develops, do not interpret inventory-equivalent growth as geometric growth.

## 11. Stage G — Grid convergence and model independence

Perform spatial refinement using at least three resolutions for representative u values.

Convergence checks must include:

- Top velocity.
- Bottom recession velocity.
- Net-thickening velocity.
- Top and bottom positions.
- Geometric crust thickness.
- Total crystalline mass.
- Dissolved-salt inventory.
- Interface-local conservation residuals.

A converged total crystal inventory alone is not sufficient.

Use time-tolerance sensitivity tests and conservative remapping tests.

Pay particular attention to the onset of bottom recession and whether the front position depends on an arbitrary porosity threshold.

Report remaining grid dependence rather than hiding it.

## 12. Required figures

Generate publication-quality SVG, PDF and PNG figures.

### Figure 1 — Two moving interfaces

Plot:

\[
z_t(t),\qquad z_b(t),\qquad h_c(t)
\]

for each u.

### Figure 2 — Corrected net-growth validation

Compare actual ALE net-thickening velocity against the corrected analytical net-growth prediction.

This is the primary new scientific figure.

### Figure 3 — Top-growth validation

Compare ALE top precipitation velocity with the corrected finite-Da top-growth law.

### Figure 4 — Solid distribution

Plot the evolving solid fraction across the liquid and porous-crust regions, clearly marking both moving boundaries.

### Figure 5 — Salt redistribution

Plot top precipitation, bottom interfacial dissolution, internal precipitation/dissolution and dissolved-salt storage.

### Figure 6 — Crust geometry and inventories

Compare geometric crust thickness, crystal mass inventory and inventory-equivalent thickness, clearly explaining why these differ.

### Figure 7 — Grid convergence

Show top/bottom/net velocity and interface-position convergence.

### Figure 8 — x0=0.08 soil forcing

Compare analytical predictions and two-moving-boundary ALE under the time-dependent soil forcing, using only the valid simulation interval.

All figures must use actual computed data.

## 13. Scientific deliverables

Create:

`analysis/p2_h_double_moving_boundary_validation/`

Include:

- `README.md`
- `GOVERNING_EQUATIONS.md`
- `BOTTOM_FRONT_CLOSURE.md`
- `CONSERVATION_AND_REMAP_TESTS.md`
- `VALIDATION_REPORT.md`
- `LIMITATIONS.md`
- Source code and test drivers.
- Parameters and input histories.
- Results in CSV/JSON.
- Interface-position and mass-inventory histories.
- Numerical convergence tables.
- Eight scientific figures.
- HPCC execution records and source hashes.

The report must explicitly answer:

1. Can the moving-front ALE recover the original fixed-bottom ALE when bottom dissolution is absent?
2. Can a bottom recession velocity be independently determined from the physical equations?
3. Does the model conservatively track crystal depletion and transfer the depleted region into the liquid domain?
4. Can u > 0 cases be calculated beyond the previous epsilon=0.999 termination?
5. Does top growth continue while the bottom recedes?
6. Does the corrected top-growth law remain accurate?
7. Does the corrected net-growth law predict the actual geometric net-thickening rate?
8. When and why do analytical and ALE predictions diverge?
9. Is the theoretical core solid fraction independently supported by the ALE spatial structure?
10. Which additional physical closures are needed before interpreting the model as a realistic soil-supported salt crust?

## 14. Execution priorities

Proceed in this order:

**P0 — Derive and verify the moving-interface physics.**

**P1 — Implement conservative bottom recession and liquid-region treatment.**

**P2 — Validate recovery, interface motion, and mass conservation.**

**P3 — Run the controlled u > 0 long-term comparisons.**

**P4 — Compare the corrected analytical top and net-growth laws.**

**P5 — Apply the time-dependent x0=0.08 soil forcing.**

Do not begin expensive 50-hour production runs before P0–P2 pass.

Use HPCC compute nodes and preserve all failed cases.

Commit the resulting code, reports, tests, figures and reproducibility evidence to:

`research/p2-h-double-moving-boundary-validation`

## 15. Critical scientific constraints

This study is primarily designed to validate our corrected growth law.

Do not introduce complex mechanics, air-gap physics or additional empirical fitting unless the simpler problem is mathematically underdetermined.

If additional interface kinetics or structural assumptions are necessary, state them explicitly and test their sensitivity.

**Do not define bottom recession using the corrected analytical net-growth law.**

**Do not use analytical \(\phi_{\mathrm{core}}\) to convert ALE crystal inventory into the geometric thickness being used for validation.**

**Do not count dissolved crystalline mass twice when moving cells from the porous crust into the liquid region.**

**Do not equate a mathematically continuous liquid region with an experimentally established mechanically stable suspended crust.**

A discrepancy between the corrected law and an independently computed two-moving-boundary ALE is a valid research result.

The primary goal is to establish whether the corrected model predicts actual net crust thickening, and under precisely which physical conditions.