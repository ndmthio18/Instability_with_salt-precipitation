# P1 model equations and units

## Scope

P1 is an isothermal, homogeneous, two-phase, three-component porous-medium model. The phases are liquid brine and gas; the components are H2O, air and NaCl. There is no solid NaCl, precipitation reaction, surface inventory or moving geometry.

## Component balances

For component \(k\), DuMuX solves a molar balance

\[
\frac{\partial N_k}{\partial t}+\nabla\cdot\mathbf J_k=0,
\qquad
N_k=\phi\sum_{\alpha\in\{l,g\}}S_\alpha c_\alpha x_\alpha^k.
\]

Here \(N_k\) has units mol m\(^{-3}\) bulk, \(c_\alpha\) is phase molar density in mol m\(^{-3}\) phase, and \(\mathbf J_k\) is mol m\(^{-2}\) s\(^{-1}\). Reported inventories are integrated and divided by horizontal projected area, giving mol m\(^{-2}\).

The phase Darcy velocity is

\[
\mathbf q_\alpha=-\frac{Kk_{r\alpha}}{\mu_\alpha}
(\nabla p_\alpha-\rho_\alpha\mathbf g).
\]

The component flux contains phase advection and Fickian diffusion. The audited upstream `2pnc` implementation uses a mass-averaged Fick reference system internally. With `UseMoles=true`, its compositional local residual converts each diffusive mass flux to molar flux by division by the component molar mass before adding it to the molar equation. P1 does not add another conversion.

For liquid H2O--NaCl, the molecular coefficient is fixed at \(1.5\times10^{-9}\) m\(^2\) s\(^{-1}\). The porous-medium effective coefficient follows upstream Millington--Quirk:

\[
D_{\mathrm{eff},\alpha}=\phi S_\alpha^3(\phi S_\alpha)^{1/3}D_\alpha.
\]

At liquid saturation one, the apparent diffusion coefficient in the concentration equation is \(D\phi^{1/3}=1.1104972\times10^{-9}\) m\(^2\) s\(^{-1}\), which is the V4 analytical reference.

## Constitutive relations

- Homogeneous \(K=1.261386079996293\times10^{-11}\) m\(^2\), \(\phi=0.4057687912195745\).
- Van Genuchten capillarity and relative permeability: \(\alpha=3.121789225428109\times10^{-4}\) Pa\(^{-1}\), \(n=2.98602455003942\), Mualem \(L=-0.8395681012108449\), \(S_{wr}=0.1282679353072967\), \(S_{nr}=0\).
- Upstream `BrineAir` supplies salinity- and pressure-dependent liquid density/viscosity, gas properties, water activity/fugacity and phase equilibrium. Only the salt molecular diffusion coefficient is overridden.
- Temperature is 293.15 K and gravity is 9.8 m s\(^{-2}\) downward.

These functions are not identical to the frozen Richards model, which uses reference viscosity 0.0011 Pa s and its own quadratic density law. P1 does not tune `BrineAir` to force onset agreement.

## Boundaries

- Bottom: all primary variables are Dirichlet in benchmark `Reservoir` mode, using the unperturbed hydrostatic pressure, composition and phase state.
- Sides: zero component flux.
- Top: water-component loss is \(j_w=e\rho_w/M_w\); NaCl and air fluxes are zero. Thus salt is nonvolatile.
- Verification `Closed` and `Flux` modes provide, respectively, impermeable boundaries and an exact component molar-flux test.

For the two-phase hydrostatic initial state, liquid pressure follows the Case 36877 water-level convention. Gas pressure uses a second-order depth integration. Its two coefficients were determined from the same `BrineAir` state and then discretely balanced on the V3 grid; the resulting values and the residual velocity are reported, rather than treated as physical fit parameters.
