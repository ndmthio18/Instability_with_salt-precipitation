# Adaptive pore charging: numerical initialization, not a capillary model

This extension changes only activation pore-liquid allocation. `ale.hh`,
`bdf.hh`, the frozen soil application, Prat kinetics and gas-film source are
unchanged from 25a41470aac2fa52ae772d313f9ae2b35a49339f. The original saturated
crust and interface supply restrictions still apply after activation.

## Inventory and geometry

For global BOX DOF i, V_i is the sum of all incident element SCV volumes,
each visited once. The current model has constant phi=.41; the donor record
retains phi separately. A general heterogeneous phi field is not claimed.
One horizontal cell is enforced. Depth ordering uses DOF coordinates; bounding
depth is the deepest corner of any constituent rectangular SCV. Tensor-product
grading is supported through actual geometry. A donor is admitted only when its
entire BOX is above `P2H.MaxDonorDepth`; no partially cut BOX is invented. Thus
the depth-cut boundary is itself discretized and must be examined with refinement.

Let C_i be dissolved NaCl mass fraction, rho_i the existing nonlinear brine
density evaluated from the current salt mole fraction, and S_i saturation.
The capacity is rho_i phi_i V_i max(0,S_i-S_min), with S_min=.99. The code
retains an additional 1e-13 in saturation solely to avoid inverse-law roundoff
crossing the same .99 screen; this is a more restrictive numerical guard.

Actual surface precipitate gives h0=Msolid/[A rho_crystal(1-eps_birth)].
Required solution mass is rho_ALE eps_birth h0 A, where rho_ALE is the original
constant ALE **liquid** density, not crystal density. It is not replaced by a
donor density. Water and salt are withdrawn in each donor's existing ratio.
The ordered algorithm fills one physical depth level before the next; donors
at identical depth split that level's demand in proportion to their capacities.
This avoids arbitrary left/right DOF ordering breaking the 1D scope.

## Mixing and reconstruction

Transferred masses are S=sum(C_i dm_i), W=sum((1-C_i)dm_i). Initial ALE liquid
has uniform C0=S/(S+W). This is instantaneous, perfectly mixed initialization;
it is not a resolved capillary filling or interlayer diffusion process.
With S+W=rho_ALE eps_birth h0 A, the unmodified ALE initial state contains
exactly those transferred masses, up to floating-point arithmetic.

The soil retains its own density law. At each donating DOF, solve
C'=S'/(S'+W'), convert C' to mole fraction by the existing molar masses, then
S_w'=(S'+W')/(phi rho(C') V). Invert the existing retention law for pressure.
Proportional withdrawal leaves C and hence rho unchanged analytically; computed
values and the reconstructed water/salt inventories are checked independently.
Pressure and saturation do change, potentially substantially; small global
mass residuals do not establish small local physical perturbations.

## Transactional behavior and logs

Allocation operates on copied donor records and a trial solution. Insufficient
total capacity leaves the post-transport/precrust-projection state unchanged;
normal transport and actual precipitation for that accepted step are retained.
No partial extraction or ALE state is accepted. A failed coupled time step
restores the previous accepted soil and crust through the existing transaction.
Diagnostic attempts are intentionally retained even if the trial is rejected.

`donor_attempts.csv` contains every adaptive attempt, requirement, total
admissible capacity, outcome and transfer residual. `donor_profiles.csv`
contains the first refusal and every successful trial's full before/after soil
profile, individual DOF volume, capacity and actual extraction. An
`accepted_activation_attempt` event identifies which trial was committed; a
successful *trial* row alone is not proof of accepted activation. `top_only`
keeps its original allocation and diagnostics and uses attempt ID zero.
In the raw attempt CSV, `used_depth_m` is the proposed greedy traversal depth,
including a partially assembled but refused proposal. Refused rows have zero
actual transfer and zero donor count. The audited summary separates this as
`proposed_depth_m` and reports actual `used_depth_m=0` for refused allocation.
It must not be interpreted as partially withdrawn liquid on a failed attempt.

`PrecrustOnly` is an explicitly uncoupled diagnostic holding activation off.
`ControlledDonorTest` constructs a prescribed skeleton mass only inside an
allocation harness, copies the analytic initial soil state, records the trial,
and exits before transport. It never alters production initial solid inventory.

## Remaining physical interpretation

The depth bound is an assumption, not a capillary supply length. Uniform mixing,
instantaneous transfer, an arbitrary retained saturation and a constant-density
saturated ALE domain remain closures requiring assessment. Successful mass
transfer does not prove pressure continuity, available capillary power, absence
of air invasion, continuum applicability below crystal size, or convergence of
pre-activation precipitation. None of those conditions is repaired by loosening
a threshold in this development stage.
