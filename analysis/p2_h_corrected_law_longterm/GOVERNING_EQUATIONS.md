# Task 26 governing equations and diagnostic contract

Baseline: c6fb6277c4cdf1689d80f4bbdda7a421060d874a. The shared ALE/BDF
operators are used unchanged. All new code is a standalone experiment driver.

For a moving cell, V=h Delta-xi, E=epsilon V, Q=rho C E. With
H=-w epsilon_upwind, w=xi v_top, F the physical liquid mass flux and S the
relative salt flux, the conservative equations are

    dE/dt = H_left - H_right - R V/rho_cr
    S = F C_upwind + rho H C_upwind - rho epsilon D_eff grad(C)
    dQ/dt = S_left - S_right - R V.

Internal diffusion uses the two half-cell resistances in series. Physical
advection and grid motion use the SAME donor selected by F+rho H so that their
sum exactly reproduces the existing numerical flux. Separate donor selection
would change the operator and is not used. Boundary diffusion is the residual
needed to reproduce the actual boundary salt flux and kinetic top closure.

The saved `storage` diagnostic is (dQ/dt)/V at fixed ALE index; it includes
changes of pore volume and cell volume. It is NOT rho epsilon dC/dt at fixed
physical z. `advection`, `diffusion`, `mesh` and `reaction` are the signed
divergences/source in kg m-3 s-1, with storage equal to their sum by the RHS.
That algebraic equality is a decomposition, not an independent budget test.

Independent numerical checks compare reconstructed Q, solid and water
inventories against accepted-step integrated face fluxes and reaction ledgers.
Each cell has an integrated salt-face and reaction budget. Dissolved storage
plotted as a rate is a secant of reconstructed inventory over the output
interval; integrated conservation does not depend on this differencing.

Material laws: D_eff=sqrt(epsilon) D_s; a_v=36(1-epsilon)/d_grain.
The surface kinetic coefficient remains k_s=beta(1-epsilon_birth)k_r even
when the internal porosity changes. For Models A/B, mathcalD=epsilon D_eff/
(1-Ceq), ell=sqrt(mathcalD/(a_v k_r)), Da=(h/ell)^2 and
Da*=k_s^2/(mathcalD a_v k_r). The task-26 Pe convention is
Pe_ell=j ell/(rho epsilon D_eff); it differs by 1/(1-Ceq) from the
mathcalD-denominator convention in the earlier derivation.

Model A uses initial properties and measured evolving h. Model Bmean uses
volume-averaged epsilon; Bupper uses overlap-weighted epsilon in the upper
3 initial ell. The latter choices are additional closures, not exact results.
Model C freezes cell porosity and solves the linear, spatially heterogeneous
reaction-diffusion BVP with prescribed bottom B and top kinetic Robin balance.
Its linear system uses no ALE top concentration or velocity as input.

The runs begin at h=3.5 um, epsilon=.43 and C=Ceq. Initial non-equilibrium
velocity is reported, not discarded from raw files. Model comparisons exclude
t<1/jfactor s when reporting comparison-window maxima (jfactor=1 in P0).
This nominal clock-scaled cutoff is not a universal relaxation time in P1.
Initial transient exclusion is not
proof of local quasi-steadiness; the term diagnostics are retained.

P0 stopping screens are epsilon_max>=.999 (remaining local solid <=.001)
and epsilon_min<=.05 (severe pore filling). They are diagnostic thresholds,
not calibrated mechanical percolation limits. The crossing is bracketed by
the preceding accepted output and first threshold output (at most 60 s apart).
The threshold state is saved, then integration stops. CVODE may internally
step beyond an output time; no later state is used as an intact crust result.
No bottom recession, liquid-gap drainage or mechanical collapse is introduced.

All long-term predictions are actual geometric top velocities/heights.
No theoretical core solid fraction is used to convert crystal mass to height.

P1 multiplies a_v(epsilon) by the registered av multiplier, without changing
the surface kinetic coefficient. Bupper uses each case's own initial ell.
Evaporation stress tests use jfactor=1000 and 10000 and compress the output
clock and target duration by that factor (180 s and 18 s). These are numerical
Pe stress tests, not representative laboratory evaporation estimates.
An additional 1 s output sequence detects pore-filling screens in selected
P1 event reruns; it changes observation times, not reaction or transport laws.
