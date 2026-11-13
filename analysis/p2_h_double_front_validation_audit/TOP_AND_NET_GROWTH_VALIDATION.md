# Independent top and geometric net-growth comparison

H1 compares the prescribed-soil-flux finite-thickness top law directly with vt=dot(zt). H2 compares P/(rho_cr phi_core) with vnet=vt−vb=dot(zt−zb). No numerical height is reconstructed from the theoretical phi_core.

The finite-H top law uses the same birth parameters and prescribed Bsoil as the analytical derivation. The present two-region ALE also contains liquid-only storage and a moving reactive bottom, so Bsoil need not equal F+ at the crust boundary. Substituting F+ into the law is a boundary diagnostic, not the original predictive law. Any model-to-model discrepancy combines different boundary conditions, evolving coefficients and approximation errors; it is not automatically a discretization error.

With Mdot=Gtop−Db+Rplus−Dminus and phi_bar=Mcr/(rho_cr hc),

\[
\dot h_c=\frac{\dot M_{cr}}{\rho_{cr}\bar\phi}
-\frac{h_c}{\bar\phi}\dot{\bar\phi}.
\]

The first term is conversion at the instantaneous mean density; the second accounts for density evolution. The identity-derived phi derivative is explanatory and tautologically reconstructs vnet. An independent finite-difference derivative of saved phi_bar is also reported, with output-time truncation error; it is not advertised as an additional conservation proof. The difference Mdot−P includes actual soil-inflow and dissolved-storage terms, so replacing Mdot by P is another approximation.

Define Δvt=vt_ALE−vt_corr, vb_diag=vt_corr−vnet_corr, and Δvb=vb_ALE−vb_diag. Then Δvnet=Δvt−Δvb exactly. The inferred analytical bottom speed is solely a diagnostic: it is not independently specified bottom kinetics. Small net error accompanied by comparable large Δvt and Δvb is cancellation, not accurate prediction of both boundaries.

Velocity errors are normalized by the positive nonzero reference Sref/(rho_cr phi_core), not by a near-zero net speed. Actual positions, thickness, absolute errors and relative discrepancies are retained. Plotting conditional P0 results does not pass H1/H2: numerical convergence and sensitivity to additional bottom physics must first be established.


An exact diagnostic decomposition of the net-law error is also retained:

\[
v_{net}-P/(\rho_{cr}\phi_{core})=
\frac{J_{soil}-P}{\rho_{cr}\phi_{core}}
-\frac{J_{soil}-\dot M_{cr}}{\rho_{cr}\phi_{core}}
+\frac{\dot M_{cr}}{\rho_{cr}}(1/\bar\phi-1/\phi_{core})
-\frac{h_c}{\bar\phi}\dot{\bar\phi}.
\]

These are actual-supply, dissolved-storage, mean-density and structural-evolution terms. Supply and storage can be individually large and cancel; it is incorrect to call all error a core-density effect. Newly created brine storage is included in the ALE and neglected in the simple analytical closure. All four terms are computed from independently evolved geometry and physical inventories; the identity itself is not an independent validation result.
