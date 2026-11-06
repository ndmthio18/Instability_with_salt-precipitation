# Soil--crust interface conservation

## 1. Purpose and control-volume convention

This document derives a conservative coupling between a fixed DuMuX soil domain and an unmeshed external porous salt structure. The derivation is independent of a particular kinetic closure. It distinguishes identities required by conservation from constitutive assumptions imported from Prat (2026).

Consider one surface patch. The soil--crust base is fixed at (z=0); (z) is positive upward. The crust occupies (0<z<h(t)). Unless stated otherwise, quantities are per unit **projected soil area** and lateral crust transport is zero. Positive soil flux is upward from soil to crust. Positive atmospheric flux is outward from crust to gas.

### 1.1 Fluxes, inventories, and units

| Symbol | Definition | SI unit |
|---|---|---|
| (F_{w,\mathrm{soil}}) | water-component mass flux from soil to crust | kg H2O m(^{-2}) s(^{-1}) |
| (F_{s,\mathrm{soil}}) | total dissolved-NaCl mass flux from soil to crust, advection plus diffusion | kg NaCl m(^{-2}) s(^{-1}) |
| (j_{\mathrm{atm}}) | water-vapor flux through the moving crust surface, relative to that surface | kg H2O m(^{-2}) s(^{-1}) |
| (M_{w,c}) | water stored in crust pore liquid | kg H2O m(^{-2}) |
| (M_{s,d}) | dissolved NaCl stored in crust pore liquid | kg NaCl m(^{-2}) |
| (M_{s,t}) | solid NaCl assigned to creation of new top layers | kg NaCl m(^{-2}) |
| (M_{s,i}) | net solid NaCl assigned to precipitation within already existing crust | kg NaCl m(^{-2}) |
| (h) | local vertical thickness; use normal position for a general surface | m |
| \(\varepsilon_c\) | crust porosity | 1 |
| (S_c) | liquid saturation in crust pore space | 1 |
| (v_s=dh/dt) | outward growth speed for a flat patch | m s(^{-1}) |

Mass units are used because Prat's (C) is a mass fraction. A molar implementation may divide every salt mass by (M_{\mathrm{NaCl}}) and every water mass by (M_{\mathrm{H2O}}), but must not mix the two bases.

Let (Y_s=C) and (Y_w=1-C) for a binary liquid. Define bulk component densities

\[
m_w=\varepsilon_cS_c\rho_lY_w,
\qquad
m_{s,d}=\varepsilon_cS_c\rho_lY_s.
\]

For a one-dimensional column,

\[
M_{w,c}=\int_0^{h(t)}m_w\,dz,
\qquad
M_{s,d}=\int_0^{h(t)}m_{s,d}\,dz.
\]

## 2. Reynolds-transport derivation

For any crust component with bulk density (m_i), Reynolds transport gives

\[
\frac{d}{dt}\int_0^{h(t)}m_i\,dz
=\int_0^{h(t)}\frac{\partial m_i}{\partial t}\,dz
+m_i(h^-,t)v_s.
\]

The last term is not optional: it is the component stored in the newly swept volume. If the base also moves at speed (v_b), an additional term (-m_i(0^+,t)v_b) appears.

### 2.1 Water component

With no internal water source,

\[
\frac{\partial m_w}{\partial t}+\frac{\partial J_w}{\partial z}=0.
\]

Integration over the moving crust and use of the moving-top jump condition give

\[
\boxed{F_{w,\mathrm{soil}}=j_{\mathrm{atm}}+\frac{dM_{w,c}}{dt}}.
\]

This is the proposed balance in the task, provided \(j_{\mathrm{atm}}\) is the water flux **relative to the moving interface**. If \(j_{\mathrm{lab}}\) is instead the Eulerian gas-side flux through a fixed spatial surface and the swept gas contains water vapor density \(\rho_v\), then

\[
j_{\mathrm{atm}}=j_{\mathrm{lab}}-\rho_vv_s,
\]

so the exact balance is

\[
F_{w,\mathrm{soil}}
=j_{\mathrm{lab}}-\rho_vv_s+\frac{dM_{w,c}}{dt}.
\]

This is the component-level counterpart of the (-\rho_vdh/dt) term in Prat equation (6). Neglecting it is normally harmless numerically but must be stated, not hidden in a sign convention.

For constant \(\varepsilon_c,S_c,\rho_l,C\), new pore volume alone contributes

\[
\left.\frac{dM_{w,c}}{dt}\right|_{\mathrm{new\ volume}}
=\varepsilon_cS_c\rho_l(1-C_{\mathrm{new}})v_s.
\]

Thus (F_{w,\mathrm{soil}}=j_{\mathrm{atm}}) is generally false during growth even for a saturated crust.

### 2.2 Total NaCl

Let (r_i\) be the positive volumetric rate of dissolved NaCl converted to internal solid, in kg m(^{-3}) s(^{-1}). Then

\[
\frac{\partial m_{s,d}}{\partial t}
+\frac{\partial J_{s,d}}{\partial z}=-r_i,
\qquad
\frac{\partial m_{s,i}}{\partial t}=r_i
\]

in existing material. With no volatile salt and no lateral salt flux,

\[
\boxed{
F_{s,\mathrm{soil}}
=\frac{d}{dt}\left(M_{s,t}+M_{s,i}+M_{s,d}\right).
}
\]

The moving-volume contribution to dissolved storage is

\[
\left.\frac{dM_{s,d}}{dt}\right|_{\mathrm{new\ volume}}
=\varepsilon_cS_c\rho_l C_{\mathrm{new}}v_s.
\]

For a constant-porosity new layer whose skeleton is entirely halite,

\[
\left.\frac{dM_{s,t}}{dt}\right|_{\mathrm{new\ skeleton}}
=(1-\varepsilon_c)\rho_{cr}v_s.
\]

These two terms reproduce the right-hand side of Prat equation (7). They show why converting all incoming salt directly to solid overpredicts growth: some salt remains dissolved in the new pore liquid.

### 2.3 Lateral transport and a non-flat surface

For a local projected patch with height (h(x,y,t)), actual surface area is multiplied by

\[
J_\Gamma=\sqrt{1+|\nabla_{\parallel}h|^2}.
\]

If \(j_n\) is specified per actual crust-surface area and \(\mathbf Q_{i,\parallel}\) is laterally integrated crust transport, the projected balance becomes

\[
F_{i,\mathrm{soil}}
=\frac{dM_{i,c}}{dt}
+J_\Gamma j_{i,n}
+\nabla_{\parallel}\cdot\mathbf Q_{i,\parallel}.
\]

The task's independent-patch approximation sets the last term to zero. Equation (36) in Prat is local in the normal direction, but it does not justify neglecting lateral transport when neighboring columns have strongly different (h), pressure, or liquid connectivity.

## 3. Exact relation to the Prat subproblem

For the saturated one-dimensional Prat crust, use upward-positive salt flux

\[
J_s(z)=\rho C U-\rho\varepsilon_cD_{s,\mathrm{eff}}\frac{\partial C}{\partial z}.
\]

At the crust base,

\[
\boxed{
F_{s,\mathrm{soil}}=J_s(0^+)
=\rho C_bU-\rho\varepsilon_cD_{s,\mathrm{eff}}
\left.\frac{\partial C}{\partial z}\right|_{0^+}.
}
\]

For a binary liquid with mass-averaged diffusion, total solution mass flux is \(\rho U\), so

\[
\boxed{F_{w,\mathrm{soil}}=\rho U-F_{s,\mathrm{soil}}.}
\]

Under Prat's base equilibrium, (C_b=C_{\mathrm{eq}}). The top balances are equations (6)--(7), and (C(z)) follows the reaction--diffusion problem. Therefore the closed growth law is recovered only when the subgrid solution simultaneously satisfies:

1. (C_b=C_{\mathrm{eq}});
2. the moving-top total-mass balance;
3. the moving-top salt balance;
4. the outer-surface kinetic relation with \(\beta\);
5. distributed reaction/dissolution in the crust;
6. the quasi-steady, uniform-property, (Pe\ll Da) reductions.

The soil fluxes are consequences of that same solution. Equation (29) by itself supplies (v_s(j)), but not a complete soil boundary condition. To impose it consistently, reconstruct (U) from equation (6), reconstruct (C(z)) from equations (18)--(24), and evaluate the two boxed base component fluxes. The resulting flux demand must be compatible with the soil state; otherwise the closed-law assumptions fail.

## 4. Solid partition without double counting

At any instant, the physically invariant solid inventory is total solid salt. The decomposition

\[
M_{s,\mathrm{solid}}=M_{s,t}+M_{s,i}
\]

requires a material-history convention:

- (M_{s,t}): skeleton carried into the domain when the moving top creates a new material layer;
- (M_{s,i}): subsequent net reaction in material that already existed.

If porosity changes internally, the same solid cannot also be assigned to top-layer creation. A spatially unresolved subgrid can preserve the split only by carrying a deposition-age/material coordinate or two independently closed source terms. Otherwise it should track total solid only.

**[Proposed]** The safest P3 state is (M_{s,\mathrm{solid}}) plus (h). Report a top/internal partition only if the reduced reaction--diffusion problem is explicitly integrated. Do not define

\[
\dot M_{s,i}=(1-f)F_{s,\mathrm{soil}}
\]

from the analytical growth efficiency (f); that identity has not been derived.

## 5. Base boundary condition and well-posedness

### 5.1 When (C_b=C_{\mathrm{eq}}) is defensible

The base can be held at equilibrium only if:

- wetted halite exists at or immediately above the interface;
- liquid contact is continuous;
- dissolution/precipitation is fast relative to transport;
- sufficient solid is available to dissolve during undersaturation;
- the coupled flux does not require negative dissolved or solid inventory.

It is not permissible to prescribe both (C_b=C_{\mathrm{eq}}) and an arbitrary (F_{s,\mathrm{soil}}). One is a Dirichlet condition and the other a Neumann condition for the same dissolved-salt subproblem. A coupled model should instead enforce state continuity/equilibrium and flux continuity, allowing the flux to be solved, or use a mixed kinetic condition.

### 5.2 Before solid exists

For (h=0) and zero surface solid, there is no crust domain and no base halite buffer. The soil keeps its ordinary evaporation/no-salt-loss boundary. A separate nucleation or activation rule is needed. A solubility threshold crossing is only a thermodynamic eligibility condition, not a validated nucleation time.

A minimal equilibrium activation can be expressed as a complementarity constraint for a surface solid inventory (M_{s,0}):

\[
M_{s,0}\ge0,\qquad C_{\mathrm{eq}}-C_b\ge0,\qquad
M_{s,0}(C_{\mathrm{eq}}-C_b)=0,
\]

but this represents instantaneous equilibrium, not nucleation kinetics. It must first be tested in P2 without geometric growth.

### 5.3 Undersaturation and dissolution

If (C_b<C_{\mathrm{eq}}) while connected solid exists, dissolution may buffer the base until available solid is exhausted. A conservative model must reduce (M_{s,\mathrm{solid}}) and can allow (v_s<0) only under a separately validated dissolution/morphology law. Prat's equations contain dissolution through the sign of (C-C_{\mathrm{eq}}), but the simplified positive-growth law was derived for evaporation-driven growth and does not by itself define recession, fragmentation, or loss of connectivity.

### 5.4 Loss of liquid connectivity

When (S_c) or a connectivity variable falls below its hydraulic threshold, the saturated Prat closure must be disabled. Possible remaining transport through films or vapor diffusion is not specified by the paper. The model must then switch to a dry-crust resistance/drying closure or stop with a declared model-domain exit. Continuing equation (29) with (j>0) would create salt without a connected salt/water supply.

### 5.5 Supersaturation inside the soil

The imposed modeling choice

\[
R_{\mathrm{precip,soil}}=0
\]

does not prevent the numerical soil solution from becoming supersaturated. If that occurs away from the external interface, the state violates the assumed physical regime. The model may be used only up to the first interior supersaturation event, or the scientific assumption must be revised. Silently allowing arbitrarily supersaturated soil while claiming equilibrium surface crystallization is not evidence-consistent.

## 6. Minimum coupled state and boundary data

For each surface patch, a minimally well-posed conservative subgrid needs:

- soil trace variables from DuMuX: liquid pressure/saturation, temperature, dissolved composition, density, and component chemical potentials or equivalent state;
- crust geometry (h);
- water and dissolved-salt inventories (M_{w,c},M_{s,d});
- total crust solid inventory (M_{s,\mathrm{solid}}), with optional top/internal history only if closed;
- crust liquid saturation or a hydraulic-connectivity state;
- \(\varepsilon_c,a_v,D_{s,\mathrm{eff}},k_r,\rho_{cr}\), either fixed parameters or explicit evolving states;
- the atmosphere-side relation that determines (j_{\mathrm{atm}});
- one base coupling relation: equilibrium/chemical-potential continuity **or** finite-rate transfer, together with flux continuity;
- one top kinetic relation and the moving-interface balances;
- nonnegativity and phase-existence constraints.

The DuMuX boundary must receive (F_{w,\mathrm{soil}}) and (F_{s,\mathrm{soil}}) returned by this same solve. It must not independently impose (C_b), (F_s), (j), and (v_s).

## 7. Supply-limited closure: acceptable and unacceptable forms

### 7.1 Unacceptable coupled rule

\[
v_s=\min\left(v_{\mathrm{Prat}}(j),
\frac{F_{s,\mathrm{soil}}}{(1-\varepsilon_c)\rho_{cr}}\right)
\]

is not a closed physical model. It ignores dissolved storage, salt in new pore liquid, internal precipitation, water required to fill new pore volume, and the concentration/gradient changes caused by limiting salt supply.

### 7.2 Diagnostic use

The same minimum may be reported offline as a **diagnostic upper bound** if both terms and omitted inventories are explicitly labeled. It must not feed back into DuMuX or be called conservative.

### 7.3 Conservation-consistent alternatives

Two defensible alternatives are:

1. **Flux-driven reduced crust problem.** Use the DuMuX component fluxes as base Neumann data, solve transient or quasi-steady crust transport with moving-top balances, and obtain (C_b,C_{\mathrm{top}},j,v_s). The base no longer has prescribed (C_{\mathrm{eq}}) unless equilibrium emerges through a kinetic/complementarity constraint.
2. **Mixed interface/inventory problem.** Solve simultaneously for (F_w,F_s,C_b,C_{\mathrm{top}},v_s) and inventory increments subject to component balances, nonnegativity, base transfer kinetics, atmospheric evaporation, and top kinetics. In the saturated, equilibrium, quasi-steady limit it must reduce to Prat equation (28)/(29).

The second form is the recommended abstraction for a surface subgrid implementation because it can degrade gracefully when supply is insufficient. Its limiting behavior, rather than a pointwise minimum, is the validation target.
