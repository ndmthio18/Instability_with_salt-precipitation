# P2-G: solubility-limited surface inventory and Prat comparison

## Purpose

The physical simulation is conservative Richards--NaCl transport followed by
an instantaneous surface-equilibrium projection. The independent Prat height
is a diagnostic hypothesis. It neither prescribes salt extraction nor changes
the gas-film boundary. Frozen P2-E/F are retained as historical evidence.

## Inventories and surface projection

At each accepted transport step, assemble BOX dual volumes belonging to top
vertices. For each of these volumes compute liquid salt and water masses using
the frozen density law and mole-to-mass conversion. With Ceq a mass fraction,

    p = max(0, Ms - Ceq/(1-Ceq)*Mw)
    Ms_new = Ms-p; Mw_new = Mw

The surface solid inventory increases by p. Recover x from the resulting mass
fraction, then Sw=(Ms_new+Mw_new)/(phi*rho(x)*V) and invert the retention law
for pressure. Merely clipping x at fixed pressure would not conserve water.
All candidate primary variables are assembled in a temporary solution before
assignment. Newton iteration and rejected transport attempts do not add solid.
There is no prescribed salt Dirichlet boundary: transport uses zero external
salt flux, followed by this explicitly split surface sink.

Units: soil integrals are kg per unit out-of-plane thickness; surface inventory
is kg/m2, obtained by dividing summed precipitate by the horizontal top area
(width times unit thickness). Thus h_mass=Ms_surface/[rho_cr*(1-eps)].
For a closed bottom, independent combined budgets are

    Rs = Ms_soil + A*Ms_surface - Ms_initial
    Rw = Mw_soil + integral(A*j*dt) - Mw_initial.

The old P1-R2 soil-only residual contains the precipitation sink and must not
be interpreted as a combined budget. `surface_precipitation.csv` records the
combined budgets and separate projection salt/water residuals.

## Evaporation and growth comparison

Gas flux uses the existing ideal H2O-air mixture, diffusion coefficient and
Pitzer water activity from frozen P2-D. Only this thermodynamic section is
extracted; no old reaction solver is compiled into this application. The
activity parameters are at 298.15 K and their use at 293.15 K remains an
explicit approximation. Film thickness .0083 m and RH .4 are benchmark inputs.
For the first-order split scheme, film flux is held at the previous accepted
surface composition during a step; endpoint projection never retroactively
changes its evaporated-water accounting.

Solid accumulates from actual precipitation starting at zero. At the declared
continuum threshold (default 1 um, sensitivity .5/2 um), h_Prat is initialized
to the actual mass height. No salt is debited at this switch. Thereafter

    dh_Prat = Ceq*j/[rho_cr*(1-eps)*(1-Ceq)] * f(Da,Da*) * dt
    dM_demand = rho_cr*(1-eps)*dh_Prat.

The finite-Da factor uses eps=.43, rho_cr=2165, kr=.0023, beta=2,
Ds=1.3e-9, D_eff=sqrt(eps)*Ds and av=36*(1-eps)/(25 um), following the
existing benchmark. Finite Da implicitly assumes internal reaction. Its
application here is exclusively a comparison with a model having no internal
reaction. The reaction-free P2-F factor is also recorded as a reference at
the same instantaneous diagnostic height.

Demand is compared to actual precipitation over exactly the steps on which
the diagnostic is active, excluding the initial mass-height assignment.
A mass mismatch alone cannot identify its cause: internal-reaction allocation,
soil supply, transient storage, spatial resolution, or the equilibrium
approximation can all contribute. Agreement likewise does not validate a
hydraulic crust model.

## Limits

No water or dissolved salt is assigned to crust pores. h_mass is a porous
skeleton equivalent; h_Prat is an independent prediction. Calling either a
fully coupled saturated-crust height would exceed the evidence. No dissolution,
hydraulic resistance, crust grid or soil-internal precipitation is included.
Top dual volumes have finite depth, so a surface-only interpretation requires
grid and time-step convergence. No 2D instability inference follows from this
one-dimensional benchmark. The switch scale is declared, not physically fitted.
