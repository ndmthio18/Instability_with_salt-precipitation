# Candidate numerical correction and evidence gate

The historical method first solves transport with J_NaCl=0, then removes
P=max(0,Ms*-Ceq/(1-Ceq) Mw*) from only the top BOX. At fixed dt and bounded
trial concentration, P is proportional to the shrinking top BOX volume.
Thus Jprecip=P/(A dt) can vanish under surface refinement even while adjacent
BOXes become supersaturated. Taking dt→0 and top-volume→0 need not commute.
Conservation alone does not cure this splitting consistency problem.

The candidate preserves the existing **instantaneous surface equilibrium**
assumption, but enforces it during the transport solve:

    C_top <= C_eq; J_p >= 0; J_p (C_eq-C_top) = 0.
    J_total,out = j_evap + J_p
    J_NaCl,out = J_p
    J_water,out = J_total,out - J_NaCl,out = j_evap.
    M_surface,new = M_surface,old + A dt J_p.

The soil's two conservation equations are total liquid mass and dissolved
NaCl mass. Putting J_p only in the salt equation would inadvertently add
water; therefore BOTH Neumann components must include the solid salt mass.
No salt is removed by a second post-step projection in this mode. Interior
supersaturation is neither clipped nor prohibited. Top BOX storage remains
in the finite-volume equations and tends to zero normally with refinement.

Solve the scalar unknown nonnegative boundary flux outside the unchanged
DuMux Newton solve, using its converged endpoint C_top. For a wet active
surface, the previous accepted flux is a predictor. The no-precipitation
branch is retained when the endpoint is undersaturated. The gas film keeps
its original accepted-state evaluation; donor physics and ALE are unchanged.
The root solve is numerical, not a new kinetic coefficient or a fitted sink.

This is a numerical correction only conditional on the original instantaneous
surface-equilibrium assumption. It is NOT a validated nucleation model,
internal soil precipitation law or capillary supply model. Whether interior
supersaturation after activation is physically acceptable remains a separate
constitutive question. The diagnostic stage must establish evidence before
enabling this candidate in production tests.

Regression: retain `P2H.PrecipitationMode=projection`. Candidate:
`P2H.PrecipitationMode=implicit_surface`. No default switch before verification.

Physical cut diagnostics use piecewise BOX inventories. For a plane crossing
a BOX, alpha_i is that BOX's exact volume fraction above the plane. Internal
face flux contributes (alpha_out-alpha_in) F. To express this as cut influx
minus the complete physical top-boundary outflux, add (1-alpha_top) F_top to
the cut influx. This boundary reconstruction correction is reported separately
from the interior advective/diffusive face fluxes. It is not an independent
constitutive prediction inside an unresolved BOX.
