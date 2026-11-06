# Prat (2026) physics audit for a DuMuX surface subgrid model

## 1. Audit scope and evidence labels

This audit concerns Marc Prat, *The Evaporation-Driven Growth Law of Porous Salt Structures*, **Langmuir** 42 (2026) 21269--21282, DOI [10.1021/acs.langmuir.6c02211](https://doi.org/10.1021/acs.langmuir.6c02211). The checked local publisher PDF has SHA-256
`82E4F450481F7185AFAF3F1B1C2A7A60C25C16FF04B0A570FCB2C989EB23AB3A` (14 pages). Crossref metadata were queried by DOI on 2026-11-06.

Statements below are marked as:

- **[Prat]** explicitly stated or derived in the paper;
- **[DuMuX]** observed in one of the audited implementations;
- **[Proposed]** a new coupling relation or modeling decision;
- **[Unverified]** a physical closure that still needs derivation or validation.

No formula below is reconstructed from memory. Equation numbers refer to the publisher PDF.

## 2. What the paper actually models

**[Prat]** The modeled object is an already existing porous salt structure of height (h(t)), with halite as solid skeleton and a connected pore space saturated by NaCl solution. Capillarity is assumed to keep solution supplied from the porous substrate. At the fixed base (z=0), the liquid salt mass fraction is prescribed as (C=C_{\mathrm{eq}}). Evaporation occurs at the moving outer surface (z=h(t)). Dissolved salt is advected and diffused through the crust, can precipitate or dissolve on internal pore walls, and precipitates at the outer surface to create new porous layers (paper pp. 21270--21271).

This is not a nucleation model and not a model of the original soil pores. Its independent domain is the porous salt structure. The paper does not derive how an unsaturated soil produces the base flux needed to maintain (C=C_{\mathrm{eq}}).

### 2.1 Conservation equations in the salt structure

Using the paper's notation, \(C\) is dissolved NaCl **mass fraction**, \(\rho\) is solution density, \(\varepsilon_s\) is crust porosity, \(\mathbf U\) is Darcy filtration velocity, \(D_{s,\mathrm{eff}}\) is the effective salt diffusivity, \(a_v\) is internal crystal surface area per bulk volume, \(k_r\) is the precipitation/dissolution coefficient, and \(\rho_{cr}\) is halite density.

**[Prat]** Dissolved-salt conservation is equation (1):

\[
\frac{\partial(\rho\varepsilon_s C)}{\partial t}
+\nabla\!\cdot(\rho C\mathbf U)
=\nabla\!\cdot(\rho\varepsilon_sD_{s,\mathrm{eff}}\nabla C)
-a_v\rho k_r(C-C_{\mathrm{eq}}).
\]

The local porosity change is equation (2):

\[
-\rho_{cr}\frac{\partial\varepsilon_s}{\partial t}
=a_v\rho k_r(C-C_{\mathrm{eq}}).
\]

The solution-mass balance is equation (3):

\[
\frac{\partial(\rho\varepsilon_s)}{\partial t}
+\nabla\!\cdot(\rho\mathbf U)
=-a_v\rho k_r(C-C_{\mathrm{eq}}).
\]

For nearly constant \(\rho=\rho(C_{\mathrm{eq}})\), equations (1)--(3) give the conservative mass-fraction equation (5). The factor \(1-C\) in its reaction term is a consequence of combining salt and solution balances; it must not be dropped when reconstructing the derivation.

### 2.2 Moving outer-interface balances

**[Prat]** For the one-dimensional moving top, the total solution/solid/vapor mass balance is equation (6):

\[
\rho U=
\left[\rho\varepsilon_s+\rho_{cr}(1-\varepsilon_s)-\rho_v\right]\frac{dh}{dt}+j,
\qquad z=h(t),
\]

where \(j\) is the evaporation mass flux and \(\rho_v\) is external water-vapor density. The salt balance is equation (7):

\[
\rho C U-\rho\varepsilon_sD_{s,\mathrm{eff}}\frac{\partial C}{\partial z}
=\left[\rho C\varepsilon_s+\rho_{cr}(1-\varepsilon_s)\right]\frac{dh}{dt},
\qquad z=h(t).
\]

These equations explicitly account for the solution filling the pore volume of newly created crust and the solid skeleton in that same new volume. The (-\rho_v dh/dt) term is a swept-volume term associated with the moving interface. It is negligible in the subsequent paper derivation relative to liquid and crystal densities, but it identifies which definition of evaporation flux is being used.

Multiplying equation (6) by (C) and subtracting equation (7), followed by (C\approx C_{\mathrm{eq}}), gives equation (8), which links top growth to evaporation and the salt diffusive gradient. It is not an independent kinetic law.

## 3. Reaction--diffusion reduction and dimensionless groups

**[Prat]** The moving coordinate (z^*=z/h(t)) converts the growing interval to a fixed interval. The transformed equation (9) contains the Reynolds/moving-coordinate term

\[
-\rho\varepsilon_s\frac{z^*}{h}\frac{dh}{dt}\frac{\partial C}{\partial z^*}.
\]

The paper then makes four reductions:

1. solution density is approximated by \(\rho(C_{\mathrm{eq}})\);
2. dissolved-salt transport in the crust is quasi-steady because the diffusion time is small relative to growth time;
3. spatial variation of \(\varepsilon_sD_{s,\mathrm{eff}}\) is neglected;
4. (C-C_{\mathrm{eq}}) is small enough to linearize the reaction source.

With (U_{\mathrm{evap}}=j/\rho), the nondimensional groups in equation (13) are

\[
Pe=\frac{U_{\mathrm{evap}}h}{\varepsilon_sD_{s,\mathrm{eff}}},
\qquad
Da=\frac{h^2(1-C_{\mathrm{eq}})a_vk_r}
{\varepsilon_sD_{s,\mathrm{eff}}}.
\]

The boundary conditions are equations (14)--(15):

\[
C^*(0)=0,
\qquad
C^*(1)=C^*_{\mathrm{top}},
\qquad C^*=C-C_{\mathrm{eq}}.
\]

The outer-surface kinetic postulate introduces the geometrical/upscaled factor \(\beta\) in equation (20). The paper explicitly notes that this macroscopic interface relation is postulated rather than formally derived by homogenization. For the ideal cubic construction in the Appendix, \(\beta=2\).

The top Damköhler number and the height-independent Damköhler number are

\[
Da_{\mathrm{top}}
=\frac{\beta(1-\varepsilon_s)(1-C_{\mathrm{eq}})h k_r}
{\varepsilon_sD_{s,\mathrm{eff}}},
\]

\[
Da^*=
\frac{\beta^2(1-\varepsilon_s)^2(1-C_{\mathrm{eq}})k_r}
{\varepsilon_sD_{s,\mathrm{eff}}a_v}.
\]

### 3.1 General and large-(Da) growth laws

**[Prat]** The derivation of equations (18)--(27) also assumes (Pe\ll Da). The general growth law retaining the finite-(Da) correction is equation (28):

\[
v_s\equiv\frac{dh}{dt}
=\frac{C_{\mathrm{eq}}j}
{\rho_{cr}(1-\varepsilon_s)(1-C_{\mathrm{eq}})}
\frac{\sqrt{Da^*}\tanh(\sqrt{Da})}
{1+\sqrt{Da^*}\tanh(\sqrt{Da})}.
\]

For \(Da>10\), the paper takes \(\tanh(\sqrt{Da})\approx1\), giving equation (29), equivalently equation (32):

\[
v_s=
\frac{C_{\mathrm{eq}}j}
{\rho_{cr}(1-\varepsilon_s)(1-C_{\mathrm{eq}})}
\frac{\sqrt{Da^*}}{1+\sqrt{Da^*}}.
\]

The last factor is not generally one. The paper explicitly shows that it remains below one over relevant parameter ranges.

### 3.2 Local surface generalization

**[Prat]** Equation (36) generalizes the large-(Da) result to a moving surface:

\[
\frac{d\mathbf X_i}{dt}\cdot\mathbf n
=\frac{C_{\mathrm{eq}}}{\rho_{cr}(1-\varepsilon_s)(1-C_{\mathrm{eq}})}
\frac{\sqrt{Da^*}}{1+\sqrt{Da^*}}
\,\mathbf j\cdot\mathbf n.
\]

This is a local normal-speed law only after the local evaporation flux is known. It does not eliminate the need to solve the external gas-side evaporation problem, nor does it supply the soil-side water and salt fluxes.

## 4. Validity envelope

The closed growth law is applicable only when all of the following are defensible:

- an external porous salt structure already exists and has a meaningful continuum thickness;
- its pore liquid is connected and effectively saturated by capillary supply;
- the base is in equilibrium with halite, so (C_b=C_{\mathrm{eq}});
- the crust transport time is short enough for a quasi-steady concentration profile;
- \(\rho\), \(\varepsilon_s\), \(a_v\), \(D_{s,\mathrm{eff}}\), and \(k_r\) are sufficiently uniform over the local column;
- the linearized reaction term is valid;
- (Pe\ll Da) for equation (28), and additionally (Da\gtrsim10) for equation (29)/(32);
- lateral transport inside the crust is negligible for a local one-dimensional column;
- growth occurs by a self-replicating new layer of similar porosity, rather than pore blockage or a changing multiscale microstructure;
- the salt structure remains hydraulically connected to the substrate.

The paper itself identifies changing crystal size, porosity reduction, gradual crust drying, vapor resistance, viscous pressure drop, substrate drying, and soil--crust gap formation as unresolved or outside the closed-law treatment (pp. 21278--21280). It also notes that the concentration boundary layer may be thinner than a representative elementary volume, which limits a continuum interpretation.

## 5. Compatibility with the proposed DuMuX conditions

| Proposed condition | Direct use of equation (29)/(32)? | Reason |
|---|---|---|
| Unsaturated soil surface | **No, unless a connected saturated crust/base meniscus is separately established.** | Soil unsaturation does not imply a saturated crust, and the law assumes continuous capillary supply. |
| Time- and space-varying soil surface concentration | **No as a prescribed-(C_{\mathrm{eq}}) law.** | The paper fixes the crust base concentration. A varying trace requires a new base boundary condition and generally a new solution of the crust subproblem. |
| Salt supply determined by the soil solution | **Not from the growth law alone.** | Equation (29) maps (j\) to growth under an implicit supply state; it is not a map from arbitrary (F_{s,\mathrm{soil}}) to growth. |
| Insufficient liquid flow to keep the crust saturated | **No.** | Quasi-steady saturated reaction--diffusion and the capillary-supply assumption fail. |
| Evaporation and soil liquid influx differ | **Only if the difference exactly supplies moving pore-liquid and solid inventories as in equations (6)--(7).** | Setting (F_w=j) discards new liquid-filled pore volume. |
| Crust initially absent | **No at activation.** | (h=0) removes the continuum domain and invalidates the diffusion-time and boundary-layer arguments. A nucleation/activation rule is required. |
| Horizontally variable thickness | **Conditionally.** | Equation (36) permits local normal growth, but surface-area mapping, the gas-side flux, and any lateral crust transport must be treated consistently. |

## 6. Internal precipitation and double-counting risk

**[Prat]** Equation (1) contains distributed precipitation/dissolution on internal crystal area. Equation (2) converts that reaction to local porosity change. Equations (7), (20), and the Appendix describe creation of solid at the moving outer surface. Thus the paper contains both a volumetric reaction and top-layer growth.

The analytical efficiency factor

\[
f(Da,Da^*)=\frac{\sqrt{Da^*}\tanh(\sqrt{Da})}
{1+\sqrt{Da^*}\tanh(\sqrt{Da})}
\]

emerges from the coupled reaction--diffusion solution and outer-surface kinetic condition. It already reflects the influence of distributed reaction on the concentration gradient and top growth. It is therefore inconsistent to use (f) for top growth and independently prescribe a second internal-precipitation fraction from the same incoming salt without re-solving the mass balance.

Conversely, (1-f) is **not established by the paper as an exact mass fraction deposited inside existing crust pores**. Such an interpretation would require integration of the volumetric reaction, inclusion of dissolved storage and moving pore volume, and a precise definition of the reference salt supply. The paper does not provide that partition identity.

**[Proposed]** A subgrid model should either:

1. solve the reduced crust boundary-value problem and integrate its volumetric reaction to obtain (P_{\mathrm{internal}}), while the moving boundary condition gives (P_{\mathrm{top}}); or
2. use the closed growth law for top growth and track only **total** crust salt unless a separately derived partition closure is available.

It must not apply the closed efficiency and an independent internal sink to the same salt flux.

## 7. Direct-law decision

The Prat law is a defensible **effective boundary model** only within its validity envelope and only when coupled with the implied component balances. It is not directly applicable at first precipitation, during undersupply, after hydraulic disconnection, or under a freely evolving soil-side concentration.

A pointwise `min(Prat growth demand, available soil salt flux)` is at most a diagnostic cap. It is not conservation-consistent because it changes (v_s) without simultaneously updating (C_b), (C_{\mathrm{top}}), (U), dissolved storage, pore-water creation, and internal reaction. A conservative supply-limited model must instead solve a flux-driven or mixed-boundary crust subproblem, as derived in `INTERFACE_CONSERVATION.md`.
