# Conservative two-boundary derivation

## Target and status

Target: independent laboratory-frame top/bottom positions and their difference,
with NaCl and water conserved in the liquid and solid-bearing regions.
Status: COHERENT AFTER EXTRA ASSUMPTION. Conservation determines interface jumps;
an explicitly assumed interfacial kinetic law supplies the missing recession
closure. Neither corrected growth formula is an evolution equation here.

## Invariant object, assumptions and notation

The invariant objects are total NaCl (dissolved plus crystalline) and water per
unit horizontal area. Isothermal, constant liquid/crystal densities, saturated
continuous liquid contact, stationary rigid skeleton, no air or mechanical
deformation. These conditions do not prove mechanical support of a real crust.
The lower liquid interval is [0,b], the crust is [b,t], h=t-b. Upward physical
liquid mass flux is m, downward-positive diffusive salt flux F=rho eps D C_z.
R>0 denotes internal precipitation; D_b>0 interfacial dissolution.
At the interface, minus means liquid below and plus means porous crust above.
Phi=1-eps. No singular surface storage is included. C is continuous at contact.

## Derivation map

1. Local phase mass balances -> moving-cell inventories.
2. Distributional salt/water/solid balances -> interface jumps.
3. Added kinetic closure -> bdot, without an inventory-to-height conversion.
4. Integration -> independent global/region budgets.
5. Exact crystal geometry identity -> assumptions needed for H2.

## 1. Moving cells (identities)

In the crust z=b+xi h and w=bdot+xi(tdot-bdot). In the liquid z=eta b,
w=eta bdot. For each crust cell V=h Delta xi, E=eps V, Q=rho C E:

    H=-w eps_upwind,
    Edot=H_left-H_right-R V/rho_cr,
    L=m+rho H,
    S=L C_upwind-F,
    Qdot=S_left-S_right-R V.

The crystal volume is V-E, evolved conservatively rather than clipping eps.
In a liquid-only cell E=V, R=0, Qdot=S_left-S_right. Its water inventory is
rho V-Q; it is never discarded when the front moves. Molecular D=Ds in liquid;
D_eff=sqrt(eps) Ds in crust. Physical liquid flux obeys

    m_left-m_right=(1-rho/rho_cr) R V.

This relation follows by eliminating pore-volume change between phase mass
balances; it does not constitute a pressure/permeability equation.

The delivered integrator advances crystal volume and full dissolved-salt mass
Q directly. A physical-scale finite-difference Jacobian is used with CVODE BDF;
error weights use cell volume, pore volume and physical Q. Accepted states are
not projected or clipped to repair conservation. The derivative of bottom
position is identically zero when kappa_b=0; dense linear-solver roundoff may
leave a position of order 1e-33 m, tested against the 1e-14 m position tolerance.

An affine excess-inventory formulation T=Q-rho C_ref E was investigated for
conditioning, but did not pass the combined two-interface tests consistently.
It is NOT the delivered production model. Its source snapshots and failures
are preserved in provenance. The original Q-variable pilot with a 0.1 nm liquid
seed also failed at 13.6836 s; new Jacobian runs have separate case identifiers.
The independent dissolution benchmark uses D=1e-6 m2/s and a 1e-4 relative gate
against its constant-bath recession asymptote. D=1e-4 stress-test failures in
the abandoned affine formulation are retained; production D remains 1.3e-9.

## 2. Bottom jumps (identities conditional on the stated phases)

Flux relative to interface speed v_b is L=m-rho eps v_b; salt S=L C-F and
water W=L(1-C)+F. Integration through the surface gives

    L_plus-L_minus=D_b,
    S_plus-S_minus=D_b,
    W_plus-W_minus=0,
    D_b=rho_cr phi_b v_b.

Therefore physical liquid flux is generally discontinuous:

    m_plus-m_minus=(1-rho/rho_cr)D_b.

The corresponding volume-flux jump is

    q_plus-q_minus=(1/rho-1/rho_cr)D_b,
    q=m/rho.

It reflects conversion of solid volume D_b/rho_cr into liquid volume D_b/rho.
For continuous C_b, the salt jump reduces exactly to

    F_minus-F_plus=(1-C_b)D_b.

The liquid and crust half-cell conductances g_l and g_c give the contact solve

    g_c(C_b-C_c)+g_l(C_b-C_l)-D_b(C_b)(1-C_b)=0.

This scalar root is solved independently of the analytical growth law. At zero
liquid thickness with prescribed bottom F=B, replace the liquid half-cell term
by B. For a moving skeleton, D_b=rho_cr phi_b(v_b-v_s), and the solid transport
equation must also contain v_s; that extension is not silently included.

## 3. Top (same constitutive boundary as the validated ALE)

    v_t=beta rho/rho_cr k_r(C_t-Ceq),
    G=rho_cr(1-eps_birth)v_t,
    m_top=j+[rho eps_birth+rho_cr(1-eps_birth)]v_t,
    S_top=G.

Use the original nonlinear boundary root and original top transport resistance.
The old fixed-bottom operator is recovered when v_b=0 and the liquid interval
has zero length. Production retains dissolved storage and the full nonlinear
fluxes, not the linearized analytical reaction-diffusion equation.

## 4. Global and separate-region identities

    dM_cr/dt = G-D_b+integral_crust R dz,
    dM_diss/dt = S_soil-G+D_b-integral_crust R dz,
    d(M_cr+M_diss)/dt = S_soil,
    dM_water/dt = m_soil-S_soil-j.

Separate dissolved budgets are

    dM_liquid_salt/dt = S_soil-S_minus,
    dM_crust_salt/dt = S_plus-G-integral R dz.

Independent cumulative flux/reaction ledgers are compared with reconstructed
inventories. Evaluating the same RHS in rearranged form is only an algebraic
interface check, not independent numerical validation. D_b enters once at the
interface; R excludes any surface delta source. No bottom cell reaction is
integrated a second time as D_b.

## 5. H2 geometry identity and additional approximations

Define actual mean solid fraction phi_bar=M_cr/(rho_cr h). Exactly,

    hdot=(G-D_b+integral R dz)/(rho_cr phi_bar)
         -h phi_bar_dot/phi_bar.

For hdot approximately P/(rho_cr phi_core), a sufficient closure is both
phi_bar approximately constant phi_core and crystal mass gain approximately P.
The latter is the condition S_soil-M_diss_dot approximately P; it does NOT
require the two terms separately to be P and zero. Indeed, eliminating the
outer liquid flux using water conservation gives the exact relation

    M_cr_dot=(j C_soil-B)/(1-C_soil)
             +C_soil/(1-C_soil) M_water_dot-M_diss_dot.

For nearly uniform, constant C approximately Ceq, the last two storage terms
can cancel despite substantial liquid storage and S_soil differing from P.
Both contributions are therefore diagnosed instead of attributing their sum
to either one alone. A locally stable core is not enough
if changing boundary zones contribute materially to the domain mean.
Rigid translation changes t and b equally, not h; deformation would require
new solid kinematics. These are hypotheses to test, not algebraic identities.

## Boundaries and open risks

There is no hydraulic prediction, mechanical stability model, homogeneous
nucleation in the liquid, or bottom precipitation branch. The stationary
skeleton and kinetic area closure are conditional assumptions. The liquid
seed and interfacial-area multiplier require sensitivity checks. Subgrain
initial thickness retains the continuum limitations of the preceding study.
