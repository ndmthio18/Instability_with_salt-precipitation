# P0 GO/NO-GO decision

## 1. Decision

### Feasibility of one volumetric soil domain plus a surface subgrid

**GO, conditionally.** A conservative model does not require a second volumetric DuMuX crust domain. The soil can remain one volumetric compositional domain while each top-boundary patch owns crust thickness and water/dissolved/solid salt inventories. Conservation is feasible if the surface subgrid returns water flux, NaCl flux, and growth speed from one coupled solve.

### Direct insertion of the simplified Prat law into the current Richards boundary

**NO-GO.** The current Richards--salt residual conserves normalized (phi S_w) and (phi S_wx_s), not physical water and NaCl components. Prat equation (29)/(32) assumes an existing saturated connected crust, (C_b=C_{\mathrm{eq}}), quasi-steady transport, (Pe\ll Da), and large (Da). It cannot be activated at the first solubility crossing or supplied with an independently prescribed soil salt flux.

### Pointwise supply clipping

**NO-GO as a coupled closure.** Taking the smaller of Prat growth demand and a salt-flux-derived velocity omits liquid storage in new pores, internal reaction, dissolved storage, water supply, and the concentration changes that caused the supply limitation. It is permissible only as an explicitly labeled offline diagnostic bound.

## 2. Main findings

1. The exact patch balances are

   \[
   F_{w,\mathrm{soil}}=j_{\mathrm{atm}}+\frac{dM_{w,c}}{dt}
   \]

   and

   \[
   F_{s,\mathrm{soil}}=
   \frac{d}{dt}(M_{s,d}+M_{s,\mathrm{solid}}),
   \]

   when (j_{\mathrm{atm}}) is measured relative to the moving interface and lateral fluxes are zero. Laboratory-frame vapor flux requires the swept-vapor correction (-\rho_vv_s).

2. New crust volume requires both liquid-filled pore mass and solid skeleton mass. For a saturated constant-porosity layer, the new-volume rates include

   \[
   \varepsilon_c\rho_l(1-C_{\mathrm{new}})v_s,
   \quad
   \varepsilon_c\rho_lC_{\mathrm{new}}v_s,
   \quad
   (1-\varepsilon_c)\rho_{cr}v_s.
   \]

3. Prat's closed growth law is one limit of a coupled reaction--diffusion/moving-interface problem. The implied soil salt flux is

   \[
   F_{s,\mathrm{soil}}
   =\rho C_bU-\rho\varepsilon_cD_{s,\mathrm{eff}}C_z(0^+),
   \]

   and the water flux is (F_w=\rho U-F_s) for the binary mass-averaged system. These fluxes cannot be chosen independently of (j), (v_s), and (C(z)).

4. The Prat reaction--diffusion solution already influences the top-growth efficiency. The complement (1-f(Da^*)) is not proved to be the internal-deposition fraction and must not be applied as a second salt sink.

5. The current Richards solver is the stronger regression/LSA platform but the weaker physical conservation platform. The archived `2pncmin` solver is the stronger dimensional component/phase platform, but its existing mineral source precipitates inside original pores and its delivery is not a complete build environment.

6. The recommended target is a homogeneous compositional `2pnc` soil solver in the reproducible Paper 2 environment, with volumetric soil precipitation disabled and independent surface states. The frozen Richards model remains the no-precipitation/LSA regression reference.

7. Existing WP9/WP14 results do not validate post-precipitation instability. WP9 failed late full-spectrum/transport gates for the canonical small perturbation; WP14 qualified a larger-amplitude numerical stress-test branch but explicitly found it non-transferable to the canonical branch.

## 3. What can be implemented now

The following next-stage work is adequately specified:

- a homogeneous two-dimensional, no-mineral, physically dimensional compositional P1 baseline;
- independent molar water/air/NaCl inventory and boundary-flux diagnostics;
- matched no-precipitation comparisons with the frozen Richards/LSA baseline;
- phase appearance/disappearance and liquid-connectivity limiting tests;
- a design-only interface API that will later return ((F_w,F_s,j,v_s)), without yet implementing growth.

No new empirical kinetic coefficient or crust geometry parameter is required for P1.

## 4. What must be derived or validated first

Before P2/P3 coding or scientific interpretation, resolve:

1. the exact DUNE/DuMuX version and upstream `BrineAir` density/viscosity formulas needed to rebuild the compositional model;
2. a matched constitutive limit between the compositional solver and the existing LSA/Richards equations;
3. the surface equilibrium activation/dissolution complementarity in physical units;
4. the rule for transition from zero-volume surface solid to finite crust thickness;
5. whether the crust base is equilibrium-buffered or governed by finite-rate transfer when soil is undersaturated;
6. a flux-driven or mixed-boundary crust closure for insufficient water/salt supply;
7. the behavior after loss of liquid connectivity;
8. the model-domain response when soil interior supersaturation occurs despite (R_{\mathrm{precip,soil}}=0);
9. independent constraints on \(\varepsilon_c,a_v,D_{\mathrm{eff}},k_r\) before feedback studies;
10. numerical qualification of the actual perturbation/wavelength/time window used for any post-precipitation instability claim.

## 5. Exact recommended scope for the next Codex task

### Task title

**P1 -- dimensional compositional baseline and conservation qualification**

### Authorized implementation scope

1. Create a new isolated model directory on `research/prat-surface-coupling`; do not modify frozen model/result directories.
2. Port the minimum homogeneous 2-D isothermal `2pnc` component model into the Paper 2 build environment.
3. Disable all volumetric mineral reactions and omit a solid primary variable unless the framework requires a zero, inactive row.
4. Pin and record DUNE/DuMuX/compiler versions and the exact upstream fluid-property sources.
5. Implement accepted-step component budgets in mol m(^{-2}) and mol m(^{-2}) s(^{-1}).
6. Add closed-box, constant-flux, hydrostatic, 1-D/2-D, passive-diffusion, and phase-switch tests.
7. Construct one matched no-precipitation case against the frozen Richards model and LSA; report differences without tuning.
8. Stop after P1 GO/NO-GO. Do not implement surface precipitation, Prat growth, nucleation, crust feedback, or long production simulations.

### P1 acceptance criteria

- total water, air, and NaCl molar residuals are below predeclared absolute and signal-relative tolerances;
- top NaCl flux is exactly zero in the baseline;
- hydrostatic spurious flux is below the evaporation/transport signal floor;
- 1-D and laterally uniform 2-D solutions agree within discretization tolerance;
- phase switching is nonnegative and repeatable under time-step refinement;
- differences from Strategy A/LSA are attributable to declared equation/constitutive differences;
- no claim of physical validation or post-precipitation instability is made.

## 6. Claim boundary

P0 establishes mathematical feasibility and identifies a viable architecture. It does **not** establish that:

- external crust actually stays saturated over the target conditions;
- (C_b=C_{\mathrm{eq}}) holds at an unsaturated soil surface;
- the equilibrium activation time is a nucleation time;
- Prat's microstructure parameters apply to this soil/salt system;
- density instability persists after precipitation;
- any existing late-stage DuMuX mode is converged or experimentally validated.

The appropriate publication-level statement is: a conservation-consistent single-soil-domain/surface-subgrid formulation is possible, but the existing solvers do not yet implement that formulation and the original Prat law is recovered only as a restricted saturated-supply limit.
