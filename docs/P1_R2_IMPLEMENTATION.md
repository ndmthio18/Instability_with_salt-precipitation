# P1-R2 conservative Richards-salt implementation

## Scope and source boundary

P1-R2 changes only the balance-law storage, component fluxes, top-boundary units and diagnostics of the existing custom Richards-salt model. Pressure and NaCl mole fraction remain the primary variables. The Box/CVFE geometry, Van Genuchten retention, Mualem relative permeability, density and viscosity closures, gravity, eigenmode initialization and bottom reservoir condition are retained. No precipitation, solid NaCl, Prat growth or crust feedback is present.

The audit reference is commit `eb79b9882a83a11359c47b25498ccb38d6dad102`. `transport_wp8_background_v2` contains a specialized one-horizontal-cell loop and was therefore not used as the 2D execution path. The isolated source was mechanically copied from the frozen `transport_wp7_localflux_v1` standard Newton path. At the reference commit, WP7 and WP8 `model.hh` and `transportdiagnostics.hh` are byte-identical; only the specialized WP8 driver differs materially. Exact hashes and the build stack are recorded in `model/p1_r2_conservative_richards/provenance/source_manifest.json`.

The complete patch against the selected frozen standard-2D source is `model/p1_r2_conservative_richards/source_diff_from_wp7_standard_2d.patch`.

Composition molar masses and density-law coefficients are mandatory runtime parameters; rounded fallback defaults were removed. The final hardened executable was rebuilt and reproduced V1. Earlier Stage B-D runs used explicit exact parameters, so this hardening does not change their residual evaluation.

## Implemented equations and units

For porosity-saturation content \(\theta=\phi S_l\), liquid density \(\rho_l\), salt mole fraction \(x_s\), mixture molar mass \(\bar M=(1-x_s)M_w+x_sM_s\), and salt mass fraction

\[
w_s=\frac{x_sM_s}{\bar M},
\]

the two residuals are

\[
\partial_t(\theta\rho_l)+\nabla\cdot(\rho_l\mathbf q_l)=0,
\]

\[
\partial_t(\theta\rho_lw_s)+
\nabla\cdot\left(\rho_lw_s\mathbf q_l-\rho_lD_{\mathrm{eff}}\nabla w_s\right)=0.
\]

Storage is in kg m^-3 of bulk porous medium. Area-integrated Box-face fluxes are in kg s^-1 per unit out-of-plane thickness. Reconstructed water quantities are exact algebraic differences within the discrete implementation:

\[
M_w=M_{total}-M_s,\qquad F_w=F_{total}-F_s.
\]

The retained Darcy closure is

\[
\mathbf q_l=-\frac{Kk_{rl}}{\mu_l}(\nabla p_l-\rho_l\mathbf g).
\]

## Diffusion and face rules

The retained constitutive diffusivity is

\[
D_{\mathrm{eff}}=\phi^{1.5}S_l^{3.5}D_{mol}.
\]

The physical diffusive salt mass flux uses \(\nabla w_s\), not a relabelled \(\nabla x_s\):

\[
\mathbf j_s^m=-\rho_lD_{\mathrm{eff}}\nabla w_s,\qquad
\frac{dw_s}{dx_s}=\frac{M_sM_w}{\bar M^2}.
\]

The discrete face rules are:

- arithmetic density for the gravity term;
- upwind relative permeability for Darcy mobility;
- upwind density and mass fraction for advective total/salt mass fluxes;
- arithmetic face density and arithmetic \(D_{\mathrm{eff}}\) for diffusion;
- the Box basis-function gradient of \(w_s\) for the diffusive gradient.

The residual and diagnostics call the same storage and transport kernels, while the independent face audit separately sums physical boundary and internal-face contributions.

## Mixture-velocity assumption

The implementation interprets the retained Darcy volumetric velocity as the mixture mass-average advective velocity. This is an explicit closure approximation, not an exact derivation of multicomponent mixture mechanics. It neglects effects such as partial-molar-volume corrections and a separate water diffusive velocity. The equal-and-opposite component diffusion implied by the mass-average frame is imposed through \(F_w=F_{total}-F_s\). Numerical component conservation is exact for this declared closure; the closure itself is not experimentally validated here.

## Evaporation modes

The frozen Case 36877 evaporation parameter is \(E=3.065544359023483\times10^{-9}\) m s^-1.

- `mass_demand` (Mode A, primary): outward H2O mass demand \(j_w=998.203E=3.0600355758103174\times10^{-6}\) kg m^-2 s^-1; net salt flux is zero.
- `equivalent_liquid_volume` (Mode B): outward total mass flux \(F_{total}=\rho_{l,face}E\); net salt flux is zero.
- `prescribed_components`: verification-only independent total and salt mass fluxes.
- `no_flux`: verification-only closed boundary.

Mode A and Mode B are never mixed in one run. With zero net salt flux, the outward water mass flux equals the outward total mass flux even though interior advective and diffusive salt terms need not separately vanish.

## Diagnostics and future interface contract

Each accepted step records the requested and accepted time step, retry state, total/NaCl/H2O inventories, independently reconstructed boundary fluxes, cumulative transfer and budget residuals. The bottom Dirichlet reaction is reconstructed from the bottom control-volume storage change and independently summed internal SCVF flux; it is not inferred as the whole-domain residual remainder.

For future Prat coupling, the solver can expose per top patch: H2O mass flux, NaCl mass flux, pressure, liquid saturation, NaCl mass fraction and evaporation demand. The future interface identities are

\[
F_{s,soil}=\frac{dM_{s,crust}}{dt},\qquad
F_{w,soil}=j_{atm}+\frac{dM_{w,crust}}{dt}.
\]

These identities do not determine a crust growth law. No arbitrary surface salt sink has been introduced.

## Build compatibility note

The isolated CMake target defines `SUPERLU_TYPEDEF_COMPLEX` because the available sal stack combines SuperLU 7 with a pinned Dune-ISTL interface that otherwise fails at compile time. This is a build compatibility definition and does not alter the governing equations or input physics.
