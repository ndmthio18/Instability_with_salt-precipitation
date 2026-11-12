# Choice of bottom closure, before implementation

## Option 1: exact depletion is not a finite-time front in this constitutive law

At a fixed stationary material position the existing law gives

    phi_dot = (rho/rho_cr)(36 k_r/d) phi (C-Ceq),
    phi(t)=phi(0) exp[integral (rho/rho_cr)(36 k_r/d)(C-Ceq) dt].

For bounded concentration and finite time, positive initial phi stays positive.
Thus the mathematical support boundary does not recede by exact depletion.
The earlier epsilon=.999 event is a diagnostic threshold, not the zero set.
A floating-point zero or a selected cutoff cannot repair this physical
underdetermination. Tracking disconnected low-solid cells would be still less
justified. Grain-size evolution or another area law could change this result,
but would introduce another constitutive model requiring independent inputs.

## Option 2 selected: a distinct exposed-surface reaction

Assume the exposed bottom crystal area per projected area is kappa_b phi_b.
This is an additional geometric closure, not measured data. Use the already
specified intrinsic k_r and local interfacial undersaturation:

    D_b=rho k_r kappa_b phi_b max(Ceq-C_b,0),
    v_b=(rho/rho_cr) k_r kappa_b max(Ceq-C_b,0).

The second expression is obtained analytically from the solid jump. The code
never divides a flux by a nearly zero numerical phi. The positive solid trace
immediately ahead of the boundary determines mass removal, not the speed.
Finite-volume trace accuracy and residual solid inventory are monitored.
The reference kappa_b=1 means projected exposed crystal area; .5 and 2 bracket
this assumption without fitting to an analytical velocity. kappa_b=0 is a
fixed-front control. This law is not determined by conservation alone.

Internal R and exposed-surface D_b describe disjoint volume and surface
measures. This option adds a surface process to the original volume model; it
does NOT reinterpret all old bottom-layer dissolution as that surface flux.
The changed assumption is why this is conditional validation, not a unique
deduction of the experimental bottom motion.

## Liquid contact and initialization

Use a small explicit liquid interval b0=1 nm, with C=Ceq and its full water/salt
inventory included initially. The initial crust thickness remains 3.5 um;
t0=b0+h0. Check b0=.1,1,10 nm before claiming insensitivity. This regularizes a
zero-volume newly created transport region and is a specified initial condition,
not an observed experimental gap. The zero-gap, zero-recession branch is used
for exact recovery of the prior fixed-bottom model.

## Exact-zero material conversion

A conservative overlap remap can transfer a connected prefix whose crystal
mass is exactly zero into the liquid region. It must keep dissolved mass and
pore volume, and leave an isolated interior zero-solid cell in place. Positive
crystal mass is never deleted under a porosity cutoff. The selected surface
kinetic model does not need such a threshold to move the interface; the remap
is independently tested for admissible zero-solid data rather than presented
as proof that the original exponential law reaches zero in finite time.

## Registered gates

P0/P1: derive all jumps and preserve original operators. P2: recovery, pure
dissolution, simultaneous motion, exact-zero remapping and non-circularity tests
must pass before 50 h production. Check global/region budgets below 2e-9 kg/m2,
position kinematics and mesh/time convergence. A failed gate is documented and
addressed before dependent production. No value of kappa_b or phi is calibrated
to obtain agreement with the corrected law.
