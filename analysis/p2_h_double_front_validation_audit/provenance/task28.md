# Task 28 — Resolve Numerical Reliability and Complete the Independent Validation of the Corrected Salt-Crust Net-Growth Law

## 1. Repository and baseline

Repository: `ndmthio18/Instability_with_salt-precipitation`

Remote: `instability`

Baseline branch:

`research/p2-h-double-moving-boundary-validation`

Baseline commit:

`04f9bb2d856512531ebcb6265e4b46ef35cc43b0`

Create a new research branch:

`research/p2-h-double-front-validation-audit`

Preserve the previous implementation, numerical results, failed runs and validation reports.

Before changing code, review:

- `analysis/p2_h_double_moving_boundary_validation/VALIDATION_REPORT.md`
- `BOTTOM_FRONT_CLOSURE.md`
- `GOVERNING_EQUATIONS.md`
- `CONSERVATION_AND_REMAP_TESTS.md`
- `LIMITATIONS.md`
- `two_front.hh`
- `two_bdf.hh`
- `run_case.cc`
- `results/convergence.csv`
- `results/flux_resolution.csv`
- `results/summary.csv`

The overall scientific objective remains **independent validation of our corrected salt-crust growth law**, not development of an increasingly complicated ALE model.

## 2. Current status and unresolved problems

The previous double-moving-boundary ALE introduced two evolving interfaces:

\[
h_c=z_t-z_b
\]

\[
v_{\mathrm{net}}=v_t-v_b.
\]

The model completed 50 h simulations for \(u=0,0.3,0.6,0.9\).

At 50 h, the finest-grid results were approximately:

| u | ALE net velocity (nm/s) | Corrected law (nm/s) | Relative discrepancy |
|---|---:|---:|---:|
| 0 | 3.5494 | 3.3301 | 6.18% |
| 0.3 | 2.4810 | 2.3310 | 6.05% |
| 0.6 | 1.3945 | 1.3320 | 4.48% |
| 0.9 | 0.3291 | 0.3330 | 1.20% |

However, the previous review identified three important limitations.

**Problem A — Time integration robustness**

For \(u=0.6\):

- `rtol=1e-8`: completed 50 h.
- `rtol=1e-9`: CVODE failed near 2.5201 s.
- `rtol=1e-10`: CVODE failed near 7.26963 s.

A complete strict-tolerance convergence test has not been established.

**Problem B — Unconverged early interface fluxes**

For the x008 case at 1 s, the crust-side diffusion flux \(F_+\) varied between approximately:

\[
+1.58\times10^{-5}
\]

and

\[
-1.06\times10^{-4}\ \mathrm{kg\,m^{-2}s^{-1}}
\]

across spatial resolutions.

The early flux is not trustworthy without additional verification.

**Problem C — Bottom kinetic closure**

The model currently uses:

\[
D_b=\rho_l k_r\kappa_b\phi_b
\max(C_{\mathrm{eq}}-C_b,0)
\]

and

\[
v_b=\frac{\rho_l}{\rho_{\mathrm{cr}}}
k_r\kappa_b\max(C_{\mathrm{eq}}-C_b,0).
\]

This is an additional interfacial kinetic assumption.

The parameter \(\kappa_b\) must not be calibrated to reproduce the corrected analytical growth law.

The previous sensitivity experiments already demonstrated that this assumption can substantially affect net growth.

The current scientific status is NO-GO for final quantitative validation.

The purpose of Task 28 is to determine whether this status can be resolved with independently verified numerics and defensible assumptions.

---

## 3. P0 — Diagnose and resolve the CVODE failures

**Highest priority. Complete this before running another full 50 h simulation matrix.**

First reproduce the existing failures without changing material parameters.

Instrument the CVODE solver to report:

- Exact internal failure time.
- Current integration step size.
- Previous accepted step size.
- Number of nonlinear iterations.
- Nonlinear convergence failures.
- Error-test failures.
- Jacobian evaluations.
- Most restrictive error-weighted state variables.
- Minimum and maximum concentrations.
- Minimum pore volume and solid volume.
- Minimum liquid-region cell size.
- Top and bottom interface velocities.

Examine whether the failures arise from:

1. Poor scaling of conservative state variables.
2. Absolute error weights for very small liquid cells.
3. Finite-difference Jacobian perturbations.
4. Strong stiffness caused by the initial 1 nm liquid region.
5. Incorrect Jacobian derivatives near the nonlinear interface root.
6. Constraint handling near zero solid volume.
7. Time integration across rapid initial boundary-layer relaxation.

Specifically audit `two_bdf.hh`.

Check the physical scaling of Jacobian perturbations and the use of centered versus one-sided differences near state constraints.

Consider a nondimensional formulation, block-scaled variables, improved Jacobian evaluation or analytically derived interface derivatives if these are justified.

**Do not simply relax the tolerance to make the model run.**

Do not clip accepted states or modify physical material parameters to hide failures.

### P0 acceptance gate

For representative \(u=0.6\), \(u=0.9\), and x008 cases:

- Complete the intended physical simulation period under at least two independently tested tolerance settings.
- Verify conservation separately from convergence.
- Compare \(v_t\), \(v_b\), \(v_{\mathrm{net}}\), \(z_t\), \(z_b\), and \(h_c\).
- Use an appropriate nonzero reference velocity to normalize errors near \(v_{\mathrm{net}}=0\).
- Aim for less than 1% difference in the relevant velocity observables on the common time interval.
- Report any unresolved failure explicitly.

If the strict-tolerance tests still cannot pass, preserve the failed evidence and do not claim final quantitative validation.

---

## 4. P1 — Investigate the early interface diffusion spikes

The current model starts with a very thin liquid interval:

\[
b_0=1\ \mathrm{nm}.
\]

Such a thin region can introduce extremely short diffusion timescales and large boundary gradients.

Investigate whether the early flux spikes are caused by physical initial relaxation, numerical stiffness, inconsistent initial boundary conditions, or the interface discretization.

### Required tests

**Test A — Initial liquid thickness**

Repeat controlled early-time runs using:

\[
b_0=0.1,\ 1,\ 10\ \mathrm{nm}.
\]

Report initial-condition sensitivity separately from numerical convergence.

**Test B — Spatial refinement**

For \(u=0.6\) and x008, compare at least 128, 256, and 512 crust cells with correspondingly refined liquid cells.

Where necessary, also test a less stretched liquid-region mesh.

**Test C — Temporal refinement**

Use sufficiently small adaptive steps during the initial transient.

Do not assume that a small requested output interval guarantees accurate internal integration.

**Test D — Consistent initialization**

Investigate constructing an initial concentration profile that satisfies the instantaneous liquid–crust interface conditions.

Compare this with the current uniform \(C=C_{\mathrm{eq}}\) initial state.

If the initial state is not compatible with the imposed forcing, clearly distinguish its genuine transient response from numerical errors.

Do not remove the physical initial transient simply to improve agreement.

### Required outputs

Plot the following quantities over the early transient:

\[
F_+(t),\quad F_-(t),\quad C_b(t),
\]

\[
v_t(t),\quad v_b(t),\quad v_{\mathrm{net}}(t).
\]

Include the first seconds explicitly, using separate panels when needed.

Evaluate convergence in absolute physical units and normalized by the reference NaCl flux.

Report the earliest time after which both spatial and temporal flux convergence criteria are satisfied.

Do not arbitrarily discard the first hour.

The appropriate exclusion interval, if any, must be established through numerical evidence.

---

## 5. P2 — Independently audit the bottom interface physics

Once P0 and P1 are adequately controlled, examine whether the bottom recession closure is physically consistent.

Current model:

\[
D_b=\rho_l k_r\kappa_b\phi_b
\max(C_{\mathrm{eq}}-C_b,0).
\]

Do not assume \(\kappa_b=1\) is necessarily realistic.

### 5.1 Verify the interfacial salt and water balances

Independently check the Rankine–Hugoniot jump relations for moving interfaces using the accepted-step flux integrals and reconstructed inventories.

Confirm that:

- Interfacial dissolution is counted once.
- Volume reaction is not duplicated at the bottom.
- Liquid-region salt and water inventories remain conservative.
- Interface velocity is computed independently of the analytical net-growth law.
- Solid fraction immediately ahead of the recession front is used consistently.

### 5.2 Evaluate bottom kinetic sensitivity

Use:

\[
\kappa_b=0.5,\ 1,\ 2.
\]

Compare the actual:

\[
v_b(t),\quad z_b(t),\quad h_c(t)
\]

at matching physical times and numerical resolutions.

Calculate sensitivities of the net growth law error to \(\kappa_b\).

Determine whether the apparent agreement between the corrected law and ALE is robust to reasonable changes in bottom interfacial kinetics.

Do not choose the value of \(\kappa_b\) that makes the growth law agree best.

### 5.3 Clarify the scope of bottom recession

Distinguish:

- Dissolution in the crust interior.
- Dissolution directly at the exposed bottom front.
- Transport through the liquid-only layer.
- Crystal mass loss.
- Geometric movement of the active solid-bearing boundary.

A surface reaction closure is an additional model assumption, not a mathematically necessary consequence of the previous internal area law.

If bottom motion is not uniquely determined without this closure, state that explicitly.

---

## 6. P3 — Re-evaluate the corrected analytical model

After numerical verification, evaluate two hypotheses separately.

### H1 — Corrected top-growth law

\[
v_t^{\mathrm{corr}}=
\frac{C_{\mathrm{eq}}j}
{\rho_{\mathrm{cr}}(1-\epsilon_t)(1-C_{\mathrm{eq}})}
\frac{X[1-u\,\mathrm{sech}(\sqrt{Da})]}
{X+\tanh(\sqrt{Da})}.
\]

Compare directly with:

\[
v_t^{\mathrm{ALE}}.
\]

Report finite-Da effects, spatially evolving coefficients and the difference between prescribed soil diffusion and the actual crust-side interface diffusion.

Use the analytical model only within its stated boundary conditions.

### H2 — Corrected net-growth law

\[
v_{\mathrm{net}}^{\mathrm{corr}}
=
\frac{P}{\rho_{\mathrm{cr}}\phi_{\mathrm{core}}}
\]

where:

\[
P=\frac{C_{\mathrm{eq}}j-B_{\mathrm{soil}}}
{1-C_{\mathrm{eq}}}.
\]

Compare directly with:

\[
\boxed{
v_{\mathrm{net}}^{\mathrm{ALE}}
=
v_t^{\mathrm{ALE}}-v_b^{\mathrm{ALE}}.
}
\]

Do not derive the ALE net growth from the theoretical solid fraction.

Independently calculate:

\[
\bar\phi(t)=
\frac{M_{\mathrm{cr}}(t)}
{\rho_{\mathrm{cr}}h_c(t)}
\]

and evaluate the exact geometric identity:

\[
\boxed{
\dot h_c
=
\frac{\dot M_{\mathrm{cr}}}
{\rho_{\mathrm{cr}}\bar\phi}
-
\frac{h_c}{\bar\phi}
\frac{d\bar\phi}{dt}.
}
\]

Use this identity to explain the difference between the corrected net-growth law and the actual ALE net growth.

### Critical validation check — Error cancellation

The previous \(u=0.9\) result showed that top and inferred bottom-velocity discrepancies nearly cancel.

Explicitly decompose:

\[
\Delta v_{\mathrm{net}}
=
\Delta v_t-\Delta v_b.
\]

Be careful: the analytical bottom velocity inferred from corrected top minus corrected net velocity is only a diagnostic, not an independently validated bottom growth law.

Report whether a small net error arises from accurate prediction or cancellation between larger errors.

Also test whether the agreement persists over the full time history rather than only at 50 h.

---

## 7. P4 — Reassess the theoretical core solid fraction

The corrected analytical model predicts:

\[
\phi_{\mathrm{core}}=0.88988992.
\]

The previous double-front ALE predicted values around 0.84–0.85.

Investigate whether this discrepancy results from:

- Long-term spatially variable internal reaction.
- Unresolved boundary layers.
- The new bottom interface kinetics.
- Finite-thickness effects.
- Material properties and specific surface area evolution.
- Assumptions used in the stable-core derivation.

Report both the actual volume-averaged solid fraction and the core-region solid fraction.

Do not use the theoretical core solid fraction as an input to the ALE.

If the analytical value cannot be independently recovered, describe its limited validity and determine the consequences for the net-growth prediction.

Do not fit a new constant \(\phi_{\mathrm{core}}\) to the ALE data and call the model validated.

---

## 8. P5 — Production validation only after numerical gates pass

Once P0–P2 are addressed, run the necessary controlled cases:

\[
u=0,\ 0.3,\ 0.6,\ 0.9
\]

for 50 h.

Use the original reference material parameters.

Then evaluate the x008 time-dependent soil forcing over its supported 44.70875 h interval.

Use the improved implementation for every compared case.

Do not compare a newly converged numerical model with an old coarse-grid run and interpret the resulting difference as a physical mechanism.

For each run, report:

- Actual top position and speed.
- Actual bottom position and speed.
- Net geometric thickness and growth rate.
- Total crystal mass.
- Dissolved NaCl inventory.
- Crust and liquid-region concentrations.
- Top precipitation.
- Bottom dissolution.
- Internal precipitation and dissolution.
- Global and local conservation errors.
- Spatial/time convergence estimates.

Preserve failed and censored cases.

If P0 remains unresolved, do not present the full-duration simulations as final quantitative validation.

---

## 9. Required figures

Create publication-quality PNG, PDF and editable SVG figures.

**Figure 1 — Numerical stability**

CVODE failure diagnostics, step sizes and error-weighted state variables.

**Figure 2 — Early interface flux convergence**

Compare \(F_+\), \(F_-\) and \(C_b\) across mesh and temporal settings.

**Figure 3 — Double-front evolution**

Plot \(z_t(t)\), \(z_b(t)\) and \(h_c(t)\) for each u.

**Figure 4 — Independent growth-law validation**

Compare ALE and analytical top/net growth velocities.

**Figure 5 — Error decomposition**

Separate top prediction error, bottom diagnostic discrepancy and net-growth error cancellation.

**Figure 6 — Bottom kinetics sensitivity**

Compare \(\kappa_b=0.5,1,2\).

**Figure 7 — Actual versus theoretical core solid fraction**

Plot spatial solid fraction and relevant averages over time.

**Figure 8 — Final numerical convergence**

Show grid/time convergence of top, bottom and net velocities and actual geometric thickness.

All figures must use actual numerical results. No schematic or fabricated data may be presented as calculated evidence.

---

## 10. Deliverables

Create:

`analysis/p2_h_double_front_validation_audit/`

Include:

- `README.md`
- `NUMERICAL_FAILURE_DIAGNOSIS.md`
- `INTERFACE_FLUX_CONVERGENCE.md`
- `BOTTOM_KINETICS_AUDIT.md`
- `TOP_AND_NET_GROWTH_VALIDATION.md`
- `CORE_SOLID_FRACTION_AUDIT.md`
- `FINAL_REVIEW.md`
- C++ and Python source files.
- Automated tests.
- Numerical results and configuration files.
- Spatial and temporal convergence tables.
- Eight scientific figures.
- HPCC execution records, logs and source hashes.

Provide a machine-readable validation status with separate gates:

1. Mathematical conservation.
2. Moving-interface consistency.
3. Spatial convergence.
4. Temporal convergence.
5. Interface-flux convergence.
6. Bottom kinetic sensitivity.
7. Corrected top-growth law.
8. Corrected geometric net-growth law.

The overall result must remain NO-GO if any essential numerical gate is unresolved.

---

## 11. Execution requirements

Use HPCC compute nodes for production calculations.

Do not run expensive simulations on login nodes.

Keep successful and failed evidence separately.

Do not suppress failures by relaxing physical bounds or solver tolerances.

Do not modify the corrected analytical law to force agreement with ALE.

Do not calibrate the bottom kinetic coefficient using the target net-growth velocity.

Prioritize:

**P0: Numerical robustness and time convergence.**

**P1: Early interface-flux convergence.**

**P2: Bottom kinetic closure audit.**

**P3/P4: Independent validation and error decomposition.**

**P5: Final production only if numerical gates pass.**

If P0 or P1 cannot be resolved, stop the dependent validation claim and deliver a complete diagnostic report with the unresolved numerical issue.

Commit to:

`research/p2-h-double-front-validation-audit`

Push to remote:

`instability`

Report the final commit SHA, validated cases, numerical failures, uncertainties and scientific conclusions.

## Critical scientific requirement

Our objective is to determine whether the corrected analytical net-growth law is a physically meaningful approximation for actual geometric salt-crust thickening.

The comparison must not be circular.

A successful mathematical conservation test does not establish spatial or temporal convergence.

A successful 50 h simulation does not by itself validate bottom kinetic assumptions.

A small net-growth error does not establish accurate top and bottom motion when errors cancel.

**The final conclusion must state what the corrected law genuinely predicts, what requires additional bottom-interface physics, and which claims remain unsupported.**