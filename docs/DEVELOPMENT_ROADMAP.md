# Development roadmap and validation gates

## 1. Governing principle

The project may use one volumetric DuMuX domain for soil and a lower-dimensional surface state for external salt. Every stage must close water and NaCl components in physical units before adding the next closure. Numerical closure, grid/time-step sensitivity, physical validation, and nucleation validation are separate claims.

The stages are intentionally nested. A failure at one stage is not repaired by adding the next stage's complexity.

## 2. P0 -- theoretical coupling and source audit

### Governing equations and unknowns

No numerical unknowns are added. P0 defines:

- soil component fluxes (F_{w,\mathrm{soil}}) and (F_{s,\mathrm{soil}});
- crust inventories (M_{w,c},M_{s,d},M_{s,\mathrm{solid}});
- growth state (h);
- Reynolds-transport and moving-interface balances;
- the validity limit under which the reduced system recovers Prat equation (28)/(29).

### Boundary conditions and simplifications

- fixed soil--crust base;
- no lateral crust transport for the derivation;
- no second volumetric crust domain;
- explicit distinction between relative and laboratory-frame evaporation flux;
- (R_{\mathrm{precip,soil}}=0).

### Acceptance and decision

P0 passes if component flux signs/units close, moving-volume storage is explicit, the base is not over-constrained, and the two code bases are classified by actual source behavior. This audit passes that theoretical gate, with implementation caveats recorded in `P0_GO_NO_GO.md`.

## 3. P1 -- homogeneous no-precipitation compositional baseline

### Scientific purpose

Establish a reproducible dimensional soil component model before introducing surface salt. This stage asks whether the target compositional architecture reproduces the accepted no-precipitation base state and instability observables under matched assumptions.

### Governing equations

In the soil domain (Omega_s), solve component conservation in molar units:

\[
\frac{\partial N_k}{\partial t}+\nabla\cdot\mathbf J_k=0,
\qquad k\in\{\mathrm{H_2O,air,NaCl}\},
\]

with the selected two-phase Darcy laws, capillary law, phase equilibrium, density, viscosity, and multicomponent diffusion. Set all mineral reactions to zero.

### Unknowns

- phase pressure/formulation variable;
- phase-state switch (saturation or missing-phase composition);
- dissolved NaCl composition;
- no solid salt state.

### Boundary conditions

- bottom pressure/composition matched to the accepted homogeneous benchmark;
- sides no flow;
- top water/air transfer matched to the declared evaporation model;
- top NaCl flux exactly zero.

### Simplifications

- homogeneous two-dimensional soil;
- isothermal initially;
- no crust, no surface precipitation, no interior precipitation;
- fixed porosity/permeability.

### Numerical conservation tests

1. closed-box component conservation;
2. constant prescribed component flux with analytic inventory change;
3. hydrostatic equilibrium with zero spurious flux;
4. homogeneous 1-D versus laterally uniform 2-D equivalence;
5. phase appearance/disappearance and liquid-connectivity recovery;
6. independent recomputation of volume-integrated inventories and boundary fluxes;
7. space and time refinement for the selected onset/transport observables.

### Limiting analytical/independent solutions

- hydrostatic pressure profile;
- one-dimensional diffusion of a passive salt pulse in constant properties;
- constant-flux inventory slope;
- saturated/no-gas limit against the frozen Strategy A model;
- LSA (t_{\mathrm{crit}}) and seeded growth rate under matched constitutive assumptions.

### GO/NO-GO

**GO** if physical molar balances close, the matched base state is reproduced, phase switching is robust, and onset differences are understood within a predeclared tolerance. **NO-GO** if agreement is obtained only by changing density, (k_{rw}), diffusion, perturbation, or onset definitions; that would not be a solver-only comparison.

## 4. P2 -- conservative surface precipitation without crust geometry

### Scientific purpose

Test thermodynamic activation and salt transfer to an external surface inventory without invoking Prat growth or changing evaporation geometry.

### Governing equations

Retain P1 soil equations and add a surface solid inventory (M_{s,0}(x,t)):

\[
F_{s,\mathrm{soil}}=\frac{dM_{s,0}}{dt}
\]

for a zero-volume surface store. Water evaporation remains the P1 top boundary; (h=0), (M_{w,c}=M_{s,d}=0).

For an equilibrium prototype, use a surface complementarity relation such as

\[
M_{s,0}\ge0,
\quad C_{\mathrm{eq}}-C_b\ge0,
\quad M_{s,0}(C_{\mathrm{eq}}-C_b)=0.
\]

This is a thermodynamic equilibrium closure, not nucleation kinetics.

### Unknowns

- all P1 soil unknowns;
- surface solid mass per projected area (M_{s,0});
- surface transfer flux solved by the coupled residual.

### Boundary conditions

- water/air as in P1;
- salt boundary flux is the negative of the surface inventory rate;
- do not independently prescribe both (C_b=C_{\mathrm{eq}}) and salt flux.

### Simplifications

- no crust thickness or pore liquid;
- no internal crust precipitation;
- no hydraulic/evaporation feedback;
- no interior soil precipitation.

### Tests and analytical limits

1. closed soil-plus-surface total NaCl conservation;
2. no precipitation below equilibrium;
3. exact inventory increase under imposed supersaturated flux;
4. dissolution returns mass to soil and stops at zero surface solid;
5. time-step convergence of event activation and complementarity residual;
6. zero-solid limit recovers P1 exactly.

### GO/NO-GO

**GO** if activation, precipitation, and dissolution are conservative and nonnegative without an arbitrary simultaneous Dirichlet/Neumann constraint. **NO-GO** for interpreting activation time as experimental nucleation without separate evidence. Stop if the soil interior becomes supersaturated under the imposed (R_{\mathrm{precip,soil}}=0) assumption.

## 5. P3 -- Prat surface growth with independent crust inventories

### Scientific purpose

Add geometric growth and crust liquid/solid bookkeeping while preserving one volumetric soil domain.

### Governing equations

For each surface patch:

\[
F_{w,\mathrm{soil}}=j_{\mathrm{atm}}+\frac{dM_{w,c}}{dt},
\]

\[
F_{s,\mathrm{soil}}=
\frac{d}{dt}\left(M_{s,d}+M_{s,\mathrm{solid}}\right),
\]

plus the reduced crust transport/kinetic problem that returns

\[
(F_w,F_s,C_b,C_{\mathrm{top}},v_s).
\]

In the saturated equilibrium limit it must reproduce Prat equation (28), and for (Da>10) equation (29).

### Unknowns

- P2 soil unknowns;
- (h,M_{w,c},M_{s,d},M_{s,\mathrm{solid}});
- crust saturation/connectivity state;
- algebraic (C_b,C_{\mathrm{top}},U,j,v_s) as required by the closure;
- optional internal/top solid split only if a reaction integral or material history closes it.

### Boundary conditions

- continuity of component fluxes at the soil--crust base;
- either equilibrium/chemical-potential continuity or finite-rate base transfer, not both a prescribed concentration and flux;
- atmospheric evaporation relation at the moving top;
- moving-top water/total-mass and salt balances;
- top crystal-growth kinetics.

### Simplifications

- independent local columns, no lateral crust transport;
- initially fixed \(\varepsilon_c,a_v,D_{\mathrm{eff}},k_r\);
- no crust feedback on gas-side evaporation beyond the prescribed (j);
- do not use the model at (h=0); activate from a P2 surface inventory with a declared finite seed thickness or an asymptotically derived zero-thickness limit.

### Numerical conservation tests

1. soil-plus-crust water and NaCl closure each time step;
2. explicit audit of new pore-liquid and solid-skeleton terms;
3. independence of results from splitting total solid into top/internal bookkeeping;
4. nonnegative inventories and bounded composition;
5. zero-growth and zero-evaporation limits;
6. dissolution/undersaturation limit;
7. mesh-independent surface quadrature and projected/actual-area conversion;
8. time-step refinement across activation and growth.

### Analytical limits

- reproduce the Prat finite-\(Da\) \(\tanh(\sqrt{Da})\) law;
- reproduce the large-(Da) law only when (Da>10);
- recover (v_s\to0) as (j\to0);
- recover the ideal upper bound as (Da^*\to\infty);
- match a direct one-dimensional reaction--diffusion solve for representative (Pe,Da,Da_{\mathrm{top}});
- demonstrate that a flux-limited case departs from the closed Prat law while retaining component closure.

### GO/NO-GO

**GO** if the coupled subgrid reproduces Prat in its own limit and remains conservative under finite supply. **NO-GO** if supply limitation is implemented only by clipping velocity, if (C_b) and (F_s) are both arbitrary, or if (1-f(Da^*)) is used as an unverified internal precipitation fraction.

## 6. P4 -- crust feedback on evaporation and capillary supply

### Scientific purpose

Only after P3, allow crust geometry and saturation to modify evaporation, pressure loss, and hydraulic connectivity. This stage is required to ask whether density instability can persist *after* precipitation under feedback, rather than under a one-way imposed evaporation rate.

### Governing additions

- gas-side vapor transfer resistance depending on (h), surface geometry, and possibly dry-crust diffusion;
- capillary/viscous supply relation within the crust;
- evolution of crust saturation/connectivity;
- optional evolution of \(\varepsilon_c,a_v,D_{\mathrm{eff}}\) only with independent evidence.

### Unknowns

P3 variables plus crust hydraulic pressure/supply potential or a validated reduced conductance, and a gas-side transfer state.

### Boundary conditions

The atmospheric flux is no longer prescribed independently:

\[
j_{\mathrm{atm}}=
\mathcal E(h,S_c,a_w,\mathrm{RH},T,\text{geometry}),
\]

and liquid supply follows a capillary/viscous relation. Loss of connectivity must transition to a declared dry-crust regime.

### Tests and limits

1. (h\to0) recovers P2/P1 evaporation;
2. zero added resistance recovers P3;
3. large dry resistance suppresses evaporation without creating/losing mass;
4. capillary conductance limits recover prescribed-supply and no-supply bounds;
5. independent grid/time/patch refinement;
6. comparison with measured evaporation and crust thickness/saturation data, not only a solubility crossing.

### GO/NO-GO

**GO** for post-precipitation instability claims only if the no-crust and crust-coupled solutions are numerically qualified over the claimed time/wavelength window and the feedback parameters are independently constrained. **NO-GO** if a non-converged short-wave mode, an enhanced artificial initial perturbation, or an unvalidated nucleation time determines the conclusion.

## 7. Dependencies and scientific risks

| Dependency/risk | Earliest stage | Required control |
|---|---:|---|
| Reproducible compositional DuMuX build | P1 | Pin DUNE/DuMuX/compiler and archive source manifests. |
| Exact density/viscosity/activity formulation | P1 | Unit and limiting tests; no silent switch between repositories. |
| LSA compatibility | P1 | Matched constitutive branches and separate (t_{\mathrm{crit}},t_{\mathrm{seed}},t_{\mathrm{detect}}). |
| Surface equilibrium versus nucleation | P2 | Label equilibrium activation; obtain independent nucleation evidence before validation claims. |
| Interior supersaturation under (R_{\mathrm{precip,soil}}=0) | P2 onward | Stop/flag at first violation or revise the model assumption. |
| Initial crust thickness | P3 | Derive an activation/seed rule; do not insert an arbitrary geometric parameter. |
| Crust microstructure \(\varepsilon_c,a_v,D_{\mathrm{eff}}\) | P3 | Sensitivity plus independent microscopy/CT constraints. |
| Lateral crust transport | P3/P4 | Patch-size refinement and comparison with a lateral reduced model if gradients are large. |
| Liquid disconnection/dry crust | P4 | Explicit regime switch and mass-conservative dry transport. |
| Late nonlinear numerical convergence | all interpretation stages | Carry WP9 limits; use WP14 only as a branch-specific numerical stress test. |

## 8. Recommended order of work

1. Implement and qualify P1 only.
2. Freeze the physical-unit component contract and independent budget checks.
3. Add P2 surface inventory with equilibrium activation/dissolution tests.
4. Derive and test the flux-driven P3 subgrid against the original Prat solution.
5. Add P4 feedback only after crust properties and connectivity closures are constrained.

No production campaign or post-precipitation instability interpretation is justified before P3 is conservative and the relevant nonlinear no-precipitation window is numerically qualified.
