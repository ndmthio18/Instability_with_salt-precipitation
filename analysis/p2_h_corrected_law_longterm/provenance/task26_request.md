# Task: Long-Term Physical Validation and Applicability Limits of the Corrected Porous Salt Growth Law

## 1. Repository and branch

Repository: `ndmthio18/Instability_with_salt-precipitation`

Remote: `instability`

Current validated branch:

`research/p2-h-corrected-law-validation`

Baseline commit:

`c6fb6277c4cdf1689d80f4bbdda7a421060d874a`

Create a new branch:

`research/p2-h-corrected-law-longterm`

Do not overwrite or modify previous validated results.

Review these existing studies before implementing anything:

- `analysis/p2_h_corrected_law_validation/`
- `analysis/p2_h_corrected_vs_ale_x008_20261112/`
- `model/p2_h_ale_coupling/ale.hh`
- `model/p2_h_ale_coupling/bdf.hh`
- `docs/P2_H_FORMULATION.md`

## 2. Scientific objective

Our previous study demonstrated excellent agreement between the corrected analytical top-growth law and the full ALE model in short-time controlled benchmarks.

However, those benchmarks covered only approximately 0.014–1.425 seconds and kept crust thickness, porosity and transport properties nearly unchanged.

The next objective is to determine:

**Can the corrected analytical growth law predict the top growth of a porous NaCl crust over long periods when internal precipitation, dissolution, porosity evolution and salt-storage changes become significant?**

We also want to identify which assumptions fail first.

Do not treat agreement with ALE as guaranteed. The purpose is to critically test the analytical model.

## 3. Analytical model

Use the previously verified corrected law:

\[
v_{\mathrm{corr}}=
\frac{C_{\mathrm{eq}}j}
{\rho_{\mathrm{cr}}(1-\epsilon_t)(1-C_{\mathrm{eq}})}
\frac{
X[1-u\,\mathrm{sech}(\sqrt{Da})]
}{
X+\tanh(\sqrt{Da})
}
\]

where

\[
X=\sqrt{Da^*},\qquad
u=\frac{B_{\mathrm{soil}}}{C_{\mathrm{eq}}j}
\]

\[
Da=
\frac{h^2(1-C_{\mathrm{eq}})a_vk_r}
{\epsilon D_{\mathrm{eff}}}
\]

\[
Da^*=
\frac{\beta^2(1-\epsilon_t)^2(1-C_{\mathrm{eq}})k_r}
{\epsilon D_{\mathrm{eff}}a_v}
\]

Check the precise material-property definitions against the previous derivation.

Distinguish the original derivation's constant properties from any later use of instantaneous effective properties.

Do not silently replace the analytical assumptions with fitted parameters.

The primary variable for validation is always the actual ALE top velocity \(v_t\), NOT an inventory-equivalent height.

## 4. Stage P0 — Long-term ALE study (highest priority)

Start from an existing porous NaCl crust.

Reference parameters:

- Initial height: \(h_0=3.5\,\mu\mathrm m\)
- Initial porosity: 0.43
- \(C_{\mathrm{eq}}=0.264940869768\)
- \(\rho_{\mathrm{cr}}=2165\,\mathrm{kg/m^3}\)
- \(k_r=10^{-3}\,\mathrm{m/s}\)
- \(\beta=2\)
- \(D_s=1.3\times10^{-9}\,\mathrm{m^2/s}\)
- Grain parameter: \(25\,\mu\mathrm m\)

Use the previously validated ALE conservative formulation.

Retain:

- Transient dissolved-salt storage.
- Liquid advection.
- Salt diffusion.
- Moving ALE coordinates.
- Internal precipitation and dissolution.
- Spatial porosity evolution.
- Finite-rate top precipitation.
- Complete local crystal dissolution only when permitted by the conservative state formulation.

### P0-A: Constant forcing tests

First run constant-evaporation and constant-bottom-diffusion cases to isolate long-term structural evolution.

Use:

\[
j=1.78\times10^{-5}\ \mathrm{kg\,m^{-2}s^{-1}}
\]

\[
B_{\mathrm{soil}}=uC_{\mathrm{eq}}j
\]

with:

\[
u=0,\ 0.3,\ 0.6,\ 0.9
\]

Run over physically meaningful durations, targeting 1, 5, 10, 20 and 50 hours when the model remains physically and numerically admissible.

Do not force simulations to reach 50 hours if the crust loses its solid framework or the model assumptions become invalid.

Preserve the existing baseline and record exactly when and why any simulation terminates.

### P0-B: Time-dependent forcing from x0=0.08

Use the existing \(x_0=0.08\) soil-derived forcing:

- Time-dependent evaporation \(j(t)\).
- Time-dependent downward salt diffusion \(B_{\mathrm{soil}}(t)\).

Run the ALE over the available 44.70875-hour history.

Important: this remains a one-way prescribed-forcing experiment, not a fully coupled soil–crust model.

Do not silently equate the following:

1. Soil original SCVF net salt flux.
2. Analytical supply reconstructed using \(j\) and \(B_{\mathrm{soil}}\).
3. ALE-computed bottom net salt flux.

They must be reported separately.

## 5. Compare different analytical approximations

For each ALE run, evaluate three analytical predictions without fitting to the ALE velocity.

### Model A — Original finite-Da law

Use the initial constant material parameters and the evolving ALE height.

This tests the influence of thickness evolution.

### Model B — Instantaneous-property analytical law

Evaluate Da and Da* using explicitly defined evolving material properties.

Because ALE porosity varies spatially, test at least:

- Initial/reference porosity.
- Crust volume-averaged porosity.
- Porosity representative of the upper reaction–diffusion zone.

Clearly identify each as an additional approximation rather than an exact consequence of the original derivation.

### Model C — Quasi-steady numerical BVP with frozen ALE properties

At selected times, freeze the actual ALE porosity and transport-property profiles.

Solve a quasi-steady reaction–diffusion boundary-value problem using those spatially varying coefficients and consistent boundary conditions.

This separates errors caused by assuming constant material properties from errors caused by neglecting transient and advective terms.

Where possible, also compare a quasi-steady model retaining advection.

**Important:** Do not use the ALE top velocity or top concentration as prescribed inputs to the BVP being evaluated against those observables.

The prediction must be generated independently from the boundary data and frozen material fields.

## 6. Decompose the governing-equation terms

This is essential.

The analytical law was derived by simplifying Marc Prat's transient salt-transport equation.

We need to know quantitatively whether those simplifications remain valid over long times.

For representative ALE states, calculate and compare the magnitude of:

\[
T_{\mathrm{storage}}
\]

\[
T_{\mathrm{ALE-motion}}
\]

\[
T_{\mathrm{advection}}
\]

\[
T_{\mathrm{diffusion}}
\]

\[
T_{\mathrm{reaction}}
\]

Use the actual ALE conservative spatial discretization and distinguish physical flux terms from mesh-motion terms.

Normalize each term using a consistent reaction–diffusion scale.

Calculate:

\[
Pe_\ell=
\frac{j\ell}{\rho_l\epsilon D_{\mathrm{eff}}}
\]

using a clearly documented definition of the effective transport coefficient.

Also report:

\[
Da(t),\qquad Da^*(t),\qquad
\frac{\ell(t)}{h(t)}
\]

and suitable transient-to-reaction or transient-to-diffusion measures.

Determine when the assumptions of quasi-steady transport, negligible advection, negligible mesh-motion effects and approximately constant porosity cease to be justified.

Do not infer the validity of the quasi-steady approximation solely from a slowly changing top velocity.

## 7. Investigate salt redistribution and finite crystal inventory

For every run, output:

\[
G_{\mathrm{top}}(t)
\]

\[
R_+(t)=\int_0^h\max(R,0)\,dz
\]

\[
D_-(t)=\int_0^h\max(-R,0)\,dz
\]

\[
M_{\mathrm{solid}}(t)
\]

\[
M_{\mathrm{dissolved}}(t)
\]

\[
J_{s,b}(t)
\]

and actual geometric height \(h(t)\).

Verify the integrated conservation identity:

\[
\frac{dM_{\mathrm{dissolved}}}{dt}
=
J_{s,b}-G_{\mathrm{top}}-R_++D_-
\]

Use independently accumulated face fluxes, reaction integrals and inventory changes to assess conservation.

Do not present an algebraic rearrangement of the same RHS as an independent numerical conservation test.

Track regions where:

\[
\epsilon\rightarrow1
\]

and regions where:

\[
\epsilon\rightarrow0
\]

Detect solid depletion and potential pore clogging before states become nonphysical.

The key scientific question is:

**When the soil delivers less dissolved salt, how much of the top growth is maintained through dissolution of the pre-existing crust, and for how long?**

Do not assume that such buffering can persist indefinitely.

## 8. Stage P1 — Parameter study of Da* and Pe

After completing P0, investigate the dependence on the independent reaction-distribution parameter Da* and the importance of advection.

### Da* study

Investigate at least three regimes:

\[
Da^*<1,\quad Da^*\approx1,\quad Da^*>1
\]

Change parameters systematically, recording their actual dimensional values.

In particular, distinguish changing the top kinetic parameter \(\beta\) from changing the internal reactive surface area \(a_v\).

Do not claim these changes are physically equivalent.

### Pe study

Vary evaporation and/or transport parameters to obtain different values of \(Pe_\ell\).

Determine when the corrected reaction–diffusion law becomes inaccurate.

Construct an applicability diagram involving:

\[
Da,\qquad Da^*,\qquad Pe_\ell,\qquad u
\]

Classify cases using measured analytical-vs-ALE error, not arbitrary theoretical labels.

## 9. Stage P2 — Actual soil–crust flux consistency

Only after P0 and P1, return to the DuMuX soil problem.

The prior x0=0.08 study reported cumulative salt quantities of:

- Soil SCVF net salt transport: 0.594990 kg/m².
- Analytical reconstructed supply: 0.606399 kg/m².
- ALE bottom salt transport: 0.698607 kg/m².

Investigate the reasons for these differences.

Reconstruct the actual soil–crust interface fluxes rather than treating the first internal soil SCVF as the interface without geometric justification.

Record:

- Soil liquid mass flux.
- ALE-required liquid mass flux.
- Soil advective salt flux.
- Soil diffusive salt flux.
- Soil total salt flux.
- ALE bottom advective and diffusive salt fluxes.
- Actual interface concentration.
- Salt dissolved-storage and solid-storage contributions.

Determine whether the prescribed saturated ALE state is compatible with the water supply available from the soil.

Do not impose two incompatible water fluxes simultaneously.

If full hydraulic coupling is unavailable, stop at a clearly identified interface-consistency diagnostic rather than introducing an arbitrary correction factor.

## 10. Required figures

Produce publication-quality PNG, PDF and editable SVG figures.

Figure 1: Long-term ALE top growth velocity versus corrected analytical predictions.

Figure 2: Relative growth-law error versus time, together with Da and Da* evolution.

Figure 3: Concentration and porosity depth profiles at selected physical times.

Figure 4: Top precipitation, internal precipitation, internal dissolution, and dissolved-salt storage versus time.

Figure 5: Normalized magnitude of storage, advection, diffusion, reaction and moving-grid terms.

Figure 6: Applicability map of the corrected analytical law across Da, Da*, Pe and u.

Figure 7: Original soil salt supply versus ALE bottom salt supply under x0=0.08 forcing.

Figure 8: Numerical grid and time convergence for representative long-term cases.

Use real simulated data and state which cases satisfy each comparison's boundary assumptions.

## 11. Numerical verification requirements

For the long-term cases, perform:

- At least three ALE spatial resolutions.
- At least two BDF tolerance settings.
- Sensitivity to mesh stretching.
- Independent water and NaCl budget checks.
- Checks for local crystal depletion and excessive pore filling.

Use appropriate stopping criteria based on physical states.

Do not continue computing through a depleted salt layer and interpret its liquid-only cells as a structurally intact porous crust.

Preserve failed cases and report their causes.

## 12. Scientific deliverables

Create:

`analysis/p2_h_corrected_law_longterm/`

Deliver:

- `README.md`
- `GOVERNING_EQUATIONS.md`
- `LONGTERM_VALIDATION_REPORT.md`
- `ASSUMPTION_BREAKDOWN.md`
- `SOIL_INTERFACE_AUDIT.md`
- Reproducible Python/C++ scripts.
- Input parameter tables.
- Numerical results.
- Validation and convergence tables.
- Publication-quality figures.
- Automated numerical tests.
- HPCC execution logs and commit SHA.

The scientific report must answer:

1. Does the corrected law predict actual ALE top velocity beyond the initial quasi-steady period?
2. When do the analytical and ALE predictions diverge?
3. Which approximation is primarily responsible for the divergence?
4. Does an evolving porosity correction improve the prediction?
5. Does the theoretical reaction–diffusion screening survive substantial internal dissolution?
6. Can internal crystal dissolution maintain top precipitation when external net salt supply decreases?
7. Are soil and ALE bottom fluxes physically consistent?
8. What is the demonstrated applicability range of the corrected law?

## 13. Execution priorities and gates

Priority order:

**P0: Long-term ALE + term decomposition + assumption audit.**

**P1: Da* and Pe parameter sensitivity.**

**P2: Real soil–crust interface consistency.**

Do not proceed to a more complicated model solely to obtain closer agreement.

If a stage reveals a fundamental mathematical or physical inconsistency, document it and address the cause before proceeding.

Use HPCC compute nodes for numerical production runs. Do not run expensive jobs on login nodes.

Commit all completed source code, tests, reports, figures and relevant results to:

`research/p2-h-corrected-law-longterm`

At completion, report the branch name, commit SHA, tested cases, quantitative findings, failed tests, unresolved assumptions and recommended next steps.

## Critical instruction

The previous short-time study already demonstrated that the corrected law can reproduce the ALE top velocity under tightly controlled near-steady conditions.

**The objective now is not to obtain another near-perfect agreement. The objective is to determine when and why this agreement breaks down.**

Do not calibrate the analytical law to the ALE results.

A well-documented failure of the analytical approximation is scientifically as valuable as a successful prediction.