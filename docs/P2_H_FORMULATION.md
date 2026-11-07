# P2-H: conservative ALE inventories and conditional soil exchange

## Scope and identities

P2-H starts from P2-G commit `a3dae1374705f55a08d9d146c55567a3bb1becd0`.
Its external equation/data reference is salt_grow commit
`8c0df6eddeea95a8a76c61f245cef2f891255c86`, specifically `dynamic.py`,
`M3_METHODS.md`, `M4_MATERIAL_HISTORY.md` and the archived `grid_640` trajectory.
The current external working-tree HEAD is not substituted for that reference.
No upstream branch or frozen P2-G implementation is modified.

The executable supports one horizontal soil cell. The soil grid remains fixed;
the crust has its own normalized faces xi and physical faces z=xi*h. It has
no horizontal crust gradients and is not a 2D discretization.

## Conservative moving-grid equations

For each crust cell, V_i=h*Delta_xi_i, E_i=eps_i*V_i and Q_i=rho*C_i*E_i.
The B2 closures are `D_eff=sqrt(eps)*D_s` and
`a_v=36*(1-eps)/d_sub`. Reaction is signed:
`R_i=rho*a_v*k_r*(C_i-Ceq)`; dissolution is not silently clipped.
Let w=xi*dh/dt and H=-w*eps_upwind. The SAME H enters pore-volume transport
and relative liquid/salt transport:

    dE_i/dt = H_left-H_right-R_i*V_i/rho_cr
    F_left-F_right = (1-rho/rho_cr)*R_i*V_i
    L = F+rho*H
    dQ_i/dt = S_left-S_right-R_i*V_i.

Internal S uses upwind advection with L and the series resistance of the two
adjacent half cells for diffusion. H uses the material velocity -w, including
the reversed-velocity case. At the top:

    Vfront = beta*rho/rho_cr*k_r*(C_top-Ceq)
    Jfront = rho_cr*(1-eps_birth)*Vfront
    F_top = j + [rho*eps_birth + rho_cr*(1-eps_birth)]*Vfront
    S_top = Jfront.

The top diffusion condition is reduced to the same quadratic as frozen M3,
evaluated in excess concentration with its cancellation-safe root formula.
The admissible root nearest equilibrium is selected; absence of a physical
root is an error, not a fallback concentration. When gas-film feedback is
enabled, this root and j(C_top) are solved consistently within the RHS.

M3 regression retains the specified 1 mm initial crust and Ceq=.264,
rho=1200, rho_cr=2170, eps_birth=.43, d_sub=25 um, beta=2, kr=.0023,
Ds=1.3e-9, j=1200*.01/86400. This externally supplied regression is explicitly
separate from the zero-crust soil experiment.

## Implicit integration and independent checks

The C++ module uses pinned SUNDIALS 6.7.0 CVODE BDF/Newton with a direct dense
linear solver. The official API description is available in the
[CVODE usage documentation](https://sundials.readthedocs.io/en/v6.2.0/cvode/Usage/).
No explicit time integration of the stiff crust operator is used.

The first implementation integrated M3 primitive eps/C variables at the exact
reference tolerances. It reproduced the requested observables but failed the
stricter accumulated mass budget. That result is retained as a failed gate.
The accepted formulation advances E/Q directly, along with h and cumulative
top solid, internal solid, bottom salt, bottom solution and vapor inventories.
Its global balances are linear invariants. No mass is redistributed afterward.

The reference tolerance coefficients remain rtol=1e-8, field/inventory
atol=1e-12 and height atol=1e-14. Absolute tolerances are mapped to the current
volume: a_E=atol_field*V, a_Q=atol_field*rho*E. Relative weights use E and Q.
This is a different error norm from primitive SciPy BDF and is disclosed rather
than called an identical integrator. Field, full-history, endpoint and budget
comparisons independently bound the effect.

M4 distinguishes fixed material z from fixed ALE index xi. The displayed
porosity maps are Eulerian snapshots of ALE cell averages; tracking their
minimum does not follow a fixed crystal/material location. No new material
cohort validation is inferred from these maps.

## Zero-crust initialization

Start h=0 and solid=0. Bare-soil evaporation and the conservative P2-G surface
projection create an external solid inventory. Its equivalent height only
determines when activation is attempted. At activation, all existing solid is
transferred exactly once into a uniform eps_birth skeleton:

    h0 = Msolid / [rho_cr*(1-eps_birth)].

Pore solution requires `rho*eps_birth*h0*A` liquid mass. Only top BOX volumes
donate it, proportionally to their available liquid, at their actual salt
composition. The corresponding water and dissolved salt are removed from soil;
the soil pressure and composition are reconstructed from the remaining masses.
The initial ALE dissolved/water inventories equal these transfers individually.
Failed resource or saturation checks leave both the soil and pre-crust solid
inventory unchanged by the attempted activation. The preceding accepted
bare-soil transport/precipitation step remains valid.

This instantaneous local pore charging is a declared initialization closure,
not a resolved capillary filling process. Its dependence on top-cell volume
must be exposed by grid/activation sensitivity. A 1 or 2 um threshold is below
the approximately 3.5 um estimated crystal scale and is numerical only.

## Conditional two-way mass coupling

After activation, the old P2-G precipitation projection is disabled. Given
the old jointly accepted state, a candidate DuMuX surface mass fraction Cb
drives the ALE bottom boundary. Cb replaces the unconditional Ceq reservoir:

    S_bottom = Cb*F_bottom - rho*eps_0*D_0*(C_0-Cb)/distance_0.

The demonstration admits positive upward total-liquid supply only. ALE
computes the required liquid flux from its saturated pore-volume balance.
Its integrated liquid and salt exchanges over the candidate step become
DuMuX's top Neumann total-liquid and NaCl fluxes, respectively. A scalar
interface solve enforces Cb equal to the resulting DuMuX endpoint composition.
The shared mass fluxes have opposite signs in the two domains.

The crust's exposed-top concentration controls the gas-film water loss, using
the existing P2-D thermodynamics. Crust-state-dependent evaporation and the
coupled liquid/salt sink feed back to the soil. Cb is held across an ALE macro
step at its implicit endpoint value; this outer coupling is first-order in
time and is checked by halving the macro step.

There is no measured crust permeability, pressure-drop relation or unsaturated
crust equation. The executable is therefore a CONDITIONAL saturated-crust mass
coupling, with a declared retained soil-surface Sw>=.99 supply screen. That
screen and a converged soil solve establish numerical admissibility, not a
derived hydraulic continuity law. The result must not be called a complete
soil/crust hydraulic model or extended to drying conditions without new physics.

Coupled cases retain P2-G Ceq converted from mole fraction .1 and rho_cr=2165.
The crust's fixed liquid density is the soil density evaluated at that Ceq.
These differ from the frozen standalone regression and are recorded separately.
The density approximation is not silently presented as M3's 1200 kg/m3 value.

## Transactions and budgets

Each coupled trial starts from immutable accepted soil and crust states and
owns a fresh BDF solver. Newton calls only read held interface fluxes. A reduced
soil timestep invalidates the ALE trial; both domains are retried from the old
state. Coupling failure, supply-screen failure, or an injected rejection also
discards both trials. Only a converged, screened, budget-checked pair is committed.
Trial activation events are labeled as trials; accepted history is authoritative.

The output separates soil salt, crust dissolved salt, top solid (including
activation solid), internal solid, soil water, crust water and cumulative vapor.
With closed soil bottom, independent reconstructed balances are

    Rs = Msoil_s + Mcrust_dissolved + Mcrust_solid - Ms_initial
    Rw = Msoil_w + Mcrust_w + Mvapor - Mw_initial.

Crust solid is independently reconstructed from rho_cr*sum(V-E); a separate
ledger residual compares it to initial skeleton + integrated top/internal
precipitation. No prescribed-interface reservoir is included in these global
budgets. The standalone regression includes its explicitly logged reservoir.
