# DuMuX model compatibility and architecture audit

## 1. Source provenance

### 1.1 Unsaturated density-instability source

- repository: `ndmthio18/surrogate_model`, standalone Paper 2 worktree;
- audited commit: `a9a042c72dc258e71174705bb597a39560f3ef11`;
- protected source branch: `codex/paper2-main`;
- audited copies: `model/transport_wp7_localflux_v1/` and `model/transport_wp8_background_v2/`;
- relevant scientific limits: `analysis/wp9_r96_convergence_20261104/DECISION_ZH.md`, `R24_R48_R96_COMPARISON_ZH.md`, and the later `analysis/wp14_nonlinear_timestep_20261106/phase_b/REPORT_FINAL_ZH.md`.

The original Paper 2 worktree had one unrelated untracked `analysis/wp11_long_validation_20261105/stage_b/` directory. It was not copied, staged, or modified.

### 1.2 Two-phase mineralization source

- repository: `ndmthio18/heterogeniety`, inspected read-only at local HEAD `c817504312be0b9ab1b0dbdd59b38070b82b7227`;
- archived delivery: `deliveries/RH0_fixed_174h_20261103/model_revisions/20261102_cylinder_fixed_hydraulics_v12/`;
- every audited v12 source file was last changed by commit `87dd21fea29d5689f0beb52f2334ca4d84399eb3`;
- the delivery explicitly excludes the complete DuMuX/DUNE dependency tree and executable.

References below use `HET-v12:<path>:<lines>` for that archived delivery.

## 2. Strategy A: existing Richards--salt solver

### 2.1 Primary variables and equation count

The solver has two primary variables:

\[
\mathbf u=(p_w,x_s),
\]

where `pressureIdx=0`, `moleFractionIdx=1`, and there are two equations ([`model.hh:238-244`](../model/transport_wp8_background_v2/source/appl/simpleRichardsnc/model.hh#L238)). Saturation and liquid density are derived from pressure and mole fraction ([`model.hh:49-53`](../model/transport_wp8_background_v2/source/appl/simpleRichardsnc/model.hh#L49)).

### 2.2 Exact implemented storage and fluxes

The implemented storage is explicitly described in source as normalized:

\[
S_0=\phi S_w,
\qquad
S_1=\phi S_w x_s,
\]

at [`model.hh:117-132`](../model/transport_wp8_background_v2/source/appl/simpleRichardsnc/model.hh#L117).

The water/volume flux through an SCVF is

\[
q_n=-\mathbf n^T\mathbf K(\nabla p_w-\rho_l\mathbf g)
\frac{k_{rw}(S_w)}{\mu},
\]

and the normalized salt flux is

\[
J_{s,n}=q_n x_{s,\mathrm{up}}
-\mathbf n^TD_{\mathrm{eff}}\nabla x_s.
\]

The exact shared residual kernel, including arithmetic face averages for \(D_{\mathrm{eff}}\) and \(\rho_l\), upstream mobility/composition, and SCVF area, is at [`model.hh:155-210`](../model/transport_wp8_background_v2/source/appl/simpleRichardsnc/model.hh#L155).

This is conservative in the code's **normalized length-equivalent variables**, not in physical moles or kilograms. There is no molar-density factor in storage or advective/diffusive salt flux. The diagnostic metadata confirms `inventory="phi*Sw*x"`, inventory unit `m`, and flux unit `m/s` ([`transportdiagnostics.hh:157-164`](../model/transport_wp8_background_v2/source/appl/simpleRichardsnc/transportdiagnostics.hh#L157)).

Consequences:

- (phi S_wx_s) is not dissolved NaCl molar or mass inventory when liquid molar density varies with salinity;
- (phi S_w) is liquid volume, not water-component mass;
- the model cannot directly close Prat's water and NaCl mass balances;
- numerical closure of these normalized equations is not evidence of conservative crystallization.

The code contains a molar-density accessor ([`main.cc:512-515`](../model/transport_wp8_background_v2/source/appl/simpleRichardsnc/main.cc#L512)) but does not use it in the residual.

### 2.3 Constitutive laws

- (S_w(p_w)): van Genuchten form in [`main.cc:408-427`](../model/transport_wp8_background_v2/source/appl/simpleRichardsnc/main.cc#L408). The intended saturation cap for (p_w>p_{nw}) is unreachable after an earlier unconditional return; this is a latent compatibility issue if positive capillary-pressure assumptions are violated.
- (k_{rw}): Mualem expression, with the corrected exponent-form option selected for Paper 2, in [`main.cc:444-467`](../model/transport_wp8_background_v2/source/appl/simpleRichardsnc/main.cc#L444).
- viscosity: constant input, [`main.cc:469-475`](../model/transport_wp8_background_v2/source/appl/simpleRichardsnc/main.cc#L469).
- effective diffusion: Millington--Quirk-type (phi^{1.5}S_w^{3.5}D_m), [`main.cc:477-485`](../model/transport_wp8_background_v2/source/appl/simpleRichardsnc/main.cc#L477).
- density: mole fraction is converted to mass fraction and passed to a quadratic brine-density correlation, [`main.cc:499-510`](../model/transport_wp8_background_v2/source/appl/simpleRichardsnc/main.cc#L499).

### 2.4 Boundary and initial conditions

- bottom: Dirichlet pressure from a hydrostatic profile and constant salt mole fraction ([`main.cc:186-207`](../model/transport_wp8_background_v2/source/appl/simpleRichardsnc/main.cc#L186));
- top and sides: Neumann; top imposes volume evaporation `Problem.EvaporationRate` in the first equation and zero salt flux ([`main.cc:209-227`](../model/transport_wp8_background_v2/source/appl/simpleRichardsnc/main.cc#L209));
- initial pressure: hydrostatic; salt is constant, periodic/random, or an LSA eigenprofile ([`main.cc:231-370`](../model/transport_wp8_background_v2/source/appl/simpleRichardsnc/main.cc#L231)).

The WP7/WP8 diagnostics reuse the residual's conservative SCVF kernel rather than VTK velocity ([`transportdiagnostics.hh:214-260`](../model/transport_wp8_background_v2/source/appl/simpleRichardsnc/transportdiagnostics.hh#L214)). They additionally enforce zero top/side salt flux in the current schema ([`transportdiagnostics.hh:296-301`](../model/transport_wp8_background_v2/source/appl/simpleRichardsnc/transportdiagnostics.hh#L296)). A precipitation boundary therefore requires a new diagnostic schema as well as a new physical boundary flux.

### 2.5 What can and cannot be reused

Can be reused directly for a frozen P1 regression reference:

- homogeneous unsaturated hydraulic setup and corrected (k_{rw}) branch;
- density-dependent Darcy flow and existing LSA initialization/detection machinery;
- reproducible build/run environment;
- accepted-step SCVF diagnostics and conservation-test patterns.

Cannot be reused unchanged for P2/P3:

- normalized storage and salt flux;
- first equation as a water-component mass balance;
- zero-salt top boundary;
- solubility threshold as a precipitation law;
- any claim that `moleFraction >= 0.1` predicts nucleation.

## 3. Strategy B: archived `2pncmin` implementation

### 3.1 Model class, phases, and primary variables

The v12 application inherits DuMuX `TwoPNCMin` and uses cell-centered TPFA (`CCTpfa`) (`HET-v12:source/appl/properties.hh:61-65,85`). It uses a brine--air fluid system and a compositional solid phase with NaCl plus inert granite (`properties.hh:92-109`). The formulation is `p0s1`, i.e. phase-0 pressure plus the switching variable for phase state (`properties.hh:127-128`).

The locally modified volume variables require molar balances (`HET-v12:source/dumux/porousmediumflow/2pnc/volumevariables.hh:100-108`). The phase-state switch maps the switch primary variable to saturation or a missing-phase composition, while the remaining component primary variable is liquid NaCl mole fraction (`volumevariables.hh:185-309`). The additional mineral primary variable is solid NaCl volume fraction; this is visible in the initial state (`HET-v12:source/appl/problem.hh:146-152`).

Thus the effective isothermal unknown set is state-dependent but structurally comprises:

\[
(p_0,\;s_1\ \text{or missing-phase composition},\;x_{\mathrm{NaCl}},\;\phi_{\mathrm{NaCl,solid}}).
\]

### 3.2 Component and mineral conservation

The local residual adds solid NaCl storage to the NaCl component equation and makes the mineral row algebraic when equilibrium precipitation is enabled (`HET-v12:source/appl/equilibrium_localresidual.hh:14-20`). The equilibrium constraint is

\[
n_s\ge0,\qquad g_{\mathrm{cap}}\ge0,\qquad
\min(n_s,g_{\mathrm{cap}})=0,
\]

with both (n_s) and the capacity gap in mol m(^{-3}) bulk (`HET-v12:source/appl/precipitation_equilibrium.hh:7-14`). The source evaluates this condition in **every soil SCV** (`HET-v12:source/appl/problem.hh:208-226`). Therefore this implementation describes precipitation inside the original porous medium.

The v12 fixed-hydraulics option retains solid salt in conservative storage but returns reference porosity and reference intrinsic permeability instead of allowing the mineral volume to modify them (`HET-v12:source/dumux/porousmediumflow/2pnc/volumevariables.hh:133-160,420-427`; `source/appl/spatialparams.hh:148-162`). The archived report notes a maximum equivalent solid volume variable of 10.976736 and explicitly states that it is an inventory diagnostic, not geometric pore occupancy or an external salt crust.

### 3.3 Boundary conditions and evaporation

All boundaries are Neumann in the application (`HET-v12:source/appl/problem.hh:145`). At the top, the code imposes water-component molar loss and separately vents air (`problem.hh:158-179`). At the bottom, a hydraulic-resistance/reservoir law supplies all three fluid components with donor upwinding (`problem.hh:182-200`).

Evaporation has several registered modes:

- constant or interval water mass flux;
- activity-weighted redistribution over top faces;
- a gas-film transfer mode using water activity, ambient RH, vapor density, and diffusivity (`problem.hh:50-109,117-143,158-180`).

The Pitzer-type NaCl water-activity approximation is in `HET-v12:source/appl/water_activity.hh:15-40`. The fluid-system overlay applies activity to liquid-water fugacity/vapor pressure and preserves nonvolatile salt (`transport_fluidsystem.hh:15-45`; `phase_equilibrium.hh:19-49`).

### 3.4 Density, viscosity, diffusion, and hydraulic laws

The application delegates brine/air density and viscosity to DuMuX `BrineAir`; the archived delivery does not vendor those complete upstream formula files. It locally replaces phase-composition closure and permits a scalar multiplier on liquid water--salt binary diffusion (`HET-v12:source/appl/transport_fluidsystem.hh:34-45`). Effective diffusion, phase mobility, capillary pressure, and phase composition are updated by the standard `2pnc` volume-variable machinery (`volumevariables.hh:127-160,175-318,407-443`).

The absence of the complete dependency tree means the exact upstream density/viscosity implementation and a clean rebuild cannot be verified from this archive alone. The 174 h report documents numerical conservation and an executable hash, not a self-contained build environment.

### 3.5 Reuse classification

| v12 capability | Reuse status | Reason |
|---|---|---|
| Two-phase brine/air phase switching | Reuse concept and, after version audit, code | Needed for liquid-connectivity loss and vapor feedback. |
| Molar water/air/NaCl component accounting | Reuse | Correct dimensional basis for soil-to-crust fluxes. |
| NaCl solubility conversion and complementarity pattern | Reuse after unit/temperature audit | Suitable for a **surface** equilibrium inventory, not unchanged as an SCV source. |
| Solid NaCl inventory | Reuse the accounting idea | Existing variable lives in soil SCVs; external crust must be a boundary state. |
| Volumetric precipitation source | **Do not reuse for the stated model.** | It violates (R_{\mathrm{precip,soil}}=0) and creates salt inside original pores. |
| Solid-dependent porosity/permeability | **Do not reuse.** | The intended soil is fixed; external crust geometry is separate. |
| Fixed-hydraulics solid volume as geometry | **Do not reuse.** | The archived report shows it can exceed pore capacity and is inventory-only. |
| Water activity and nonvolatile-salt phase equilibrium | Reuse with tests | Relevant to evaporation and phase behavior. |
| 3-D heterogeneous cylinder geometry and boundary supply | **Do not reuse.** | Different scientific system. |
| 174 h comparison as validation of surface crust | **Do not reuse as validation.** | It overpredicts solution loss and contains no external crust geometry. |

## 4. Architecture comparison

| Criterion | Strategy A: extend Richards--salt | Strategy B: adapt compositional `2pnc(min)` |
|---|---|---|
| Conservation consistency now | Normalized volume/salt equations only | Dimensional molar component balances plus mineral inventory |
| Primary variables | (p_w,x_s) | Phase-state-dependent (p_0,s/x), (x_{NaCl}), optional mineral state |
| Surface solid | Requires a new boundary inventory | Also requires a new boundary inventory; existing mineral is volumetric |
| Phase transition/connectivity | No explicit gas phase switch | Explicit liquid/gas phase switch |
| Density-driven convection | Present and tied to existing LSA | Present in compositional flow, but must be benchmarked against LSA |
| Solubility equilibrium | Threshold diagnostic only | Implemented complementarity pattern, currently in soil cells |
| LSA compatibility | High for the frozen model | Not established; requires a limiting-case benchmark |
| Numerical complexity | Lower | Higher: phase switching, multiple components, optional algebraic mineral row |
| Source changes for P2/P3 | Fundamental dimensional residual and boundary changes | Geometry simplification, disable soil mineralization, add boundary state |
| Reproducible build | Available in Paper 2 | Archived v12 environment incomplete; port/version audit required |
| P4 crust drying/vapor feedback | Structurally weak | Structurally suitable |

## 5. Recommendation

### 5.1 Target architecture

**Recommend Strategy B as the target physical architecture, but not a direct copy of the v12 application.** Build a homogeneous two-dimensional compositional soil model in the reproducible Paper 2 DuMuX environment, using the `2pnc` component/phase machinery and independently audited v12 constitutive overlays. Set

\[
R_{\mathrm{precip,soil}}=0
\]

everywhere in the volumetric domain. Represent precipitated salt only in surface-owned state variables. Do not map the v12 mineral volume fraction to crust thickness.

This recommendation is based on conservation and phase behavior, not code length. The post-precipitation question includes variable liquid supply and possible hydraulic disconnection, which the single-phase normalized Richards model cannot represent without a governing-equation rewrite.

### 5.2 Role of Strategy A

Keep Strategy A frozen as:

- the P1 regression oracle for the no-precipitation instability problem;
- the link to existing LSA initialization and onset definitions;
- a source of verified SCVF diagnostic/test patterns.

It is not an acceptable P2/P3 physical residual without conversion to dimensional component conservation. If Strategy B cannot be made reproducibly buildable, the fallback is a new compositional Richards/`RichardsNC` implementation with physical molar storage and flux, not insertion of a surface sink into the existing normalized salt row.

### 5.3 Required architecture gate

Strategy B is only admitted after all of the following:

1. exact DuMuX/DUNE versions and upstream fluid-property formulas are pinned;
2. a homogeneous 2-D no-mineral build is reproducible;
3. closed-box and boundary-flux water/NaCl molar balances close to solver tolerance;
4. its no-precipitation base state and density-instability onset are compared with Strategy A/LSA under matched constitutive assumptions;
5. phase switching and loss/recovery of liquid connectivity pass manufactured or limiting tests.

Until those gates pass, the architecture decision is a conservation-based recommendation, not a validated replacement solver.

## 6. Numerical evidence limits carried forward

WP9 established that the canonical (10^{-6}) perturbation had stable early seeded-mode growth but failed late full-spectrum and shallow-transport agreement: G1/G2 passed and G3--G5 failed. R96 also changed both (dx) and (dt), vertical resolution was fixed, and the restricted box excluded the global critical mode. Those results prohibit post-precipitation claims from late small-amplitude short-wave behavior.

WP14 later showed that a (10^{-4}) perturbation is a useful nonlinear numerical stress test and that, for that branch, observed time-step differences were generally smaller than horizontal-grid differences. It also changed the late modal branch and cannot replace the (10^{-6}) physical protocol. Neither work includes precipitation or a crust. Therefore no existing run validates post-precipitation instability, crust feedback, or a physical (t_{\mathrm{eff}}).
