# Bottom kinetic closure audit

## Independent balances in the moving frame

Let L+ = m+ − rho_l eps_b vb and L− = m− − rho_l vb be the relative total liquid-mass fluxes, and S+, S− the corresponding salt fluxes. Let Db=rho_cr phi_b vb. The jump relations are

\[
S_+-S_-=D_b,\qquad L_+-L_-=D_b,
\]
\[
(L_+-S_+)-(L_--S_-)=0.
\]

Thus water is continuous in the moving frame; dissolved salt increases by the crystal mass removed at the boundary. The laboratory-frame liquid mass jump is m+−m−=(1−rho_l/rho_cr)Db. These are distinct statements; continuity of diffusion alone is not valid at this moving reactive interface.

For the liquid-only interval, independently reconstruct Q_l=sum rho_l C_i V_i and W_l=rho_l zb−Q_l. Accepted integration ledgers yield

\[
\Delta Q_l=I_{soil,salt}-I_{S-},\quad
\Delta W_l=I_{soil,mass}-I_{soil,salt}-I_{L-}+I_{S-}.
\]

For the crust, reconstruct Q_c and W_c from physical pore volumes. The salt budget is

\[
\Delta Q_c=I_{S+}-I_{top}-I_{R+}+I_{R-}.
\]

The solid budget is ΔMcr=Itop−IDb+IR+−IR−. Summing all three gives total NaCl change equal to soil input. The code stores these fluxes as BDF-integrated state ledgers; it does not approximate their integrals using coarse output trapezoids. Diagnostics reconstruct inventories independently. Shared face fluxes guarantee cancellation by construction; small budget residuals therefore test implementation consistency, not empirical validity of the kinetic closure.

## Additional physics, not an analytical identity

The constitutive assumption is exposed reactive area per projected area κ_b phi_b. It gives Db=rho_l kr κ_b phi_b max(Ceq−Cb,0) and vb=(rho_l/rho_cr)kr κ_b max(Ceq−Cb,0). The solid fraction is the instantaneous trace immediately above the bottom front, not the theoretical core value. No division by a tiny phi_b is performed. A finite trace-loss screen stops cases that lose the assumed solid-bearing front.

Internal dissolution is an integral of negative volumetric reaction over existing crust cells. Surface dissolution is only Db at the moving boundary. They occupy different support and are counted once. Liquid-only transport has no crystal reaction. A connected exactly depleted prefix can be conservatively transferred to the liquid mesh; isolated zero-solid cells remain in the crust. This transfer is not itself a source or sink.

Crystal mass loss does not in general determine geometric boundary motion. The kinetic exposed-area assumption is needed in this implementation; κ_b=1 has no independent experimental calibration here. No κ value may be selected from its agreement with the target growth law. Sensitivity tests are controlled perturbations of this assumption, not its physical validation.

The analytic core constant is absent from the ALE right-hand side, surface closure, geometry update and initialization. Fixed-bottom recovery, constant-bath dissolution and conservative prefix-transfer checks are independent of that analytical growth law. Final run-specific budget and sensitivity tables are linked from the final review.


The diffusion jump is consequently F+−F−=−(1−Cb)Db at this reactive moving boundary. Diffusion continuity from a stationary, nonreactive soil/crust matching problem cannot be carried unchanged into the present bottom-front closure. Bsoil is imposed below the liquid-only layer, while F− and F+ are solved at the moving front. Their separation is explicitly plotted and budgeted.


The model has no rigid-body settling/contact-mechanics law. Its absolute front coordinates use the assumed stationary-solid frame and a fixed soil boundary below the liquid interval. The active crust thickness hc=zt−zb excludes that liquid-only gap and is invariant under a common translation of both fronts. Experimental height above the soil (zt) and active crust thickness (hc) are therefore different observables; no settling assumption is silently added to convert between them.
