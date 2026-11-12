# Interpretation limits

1. Stationary skeleton and continuous liquid contact are imposed. No force
   balance proves that an unsupported experimental crust can remain suspended.
2. Conservation does not uniquely determine the dissolution velocity. The
   exposed-area rule kappa_b phi_b is additional constitutive information;
   sensitivity is part of the result, not a fit to a target growth law.
3. The old area law a_v proportional to phi has no finite-time exact depletion
   for bounded concentration. The new surface process is a different physical
   closure, not a threshold-based reinterpretation of old numerical zeros.
4. A finite liquid seed is an initialization regularization. Its water and
   salt are explicitly counted. A failed subnanometre numerical test is kept;
   successful larger seeds do not by themselves establish a zero-seed limit.
5. The liquid region has no homogeneous crystal nucleation and the bottom
   interface law only dissolves; supersaturation does not activate a separate
   downward-growing solid front. This is a conditional first model.
6. Flow is reconstructed from saturated volume and mass balance, without
   pressure, permeability, capillary continuity or a soil water-supply solve.
   x008 remains one-way imposed forcing, not a coupled DuMuX calculation.
7. Soil B and crust-side diffusion differ after an interfacial reaction and a
   liquid transport region are introduced. Comparisons must state which flux
   appears in the analytical boundary condition.
8. Initial crust and reaction-layer lengths are below the nominal grain
   parameter; continuum convergence is not pore-scale or experimental validity.
9. A stable local core alone does not fix the domain mean solid fraction or
   boundary-layer inventory. Net geometry must be tested using both positions.
10. State screens identify limits of this calculation, not measured fracture,
    collapse or percolation thresholds. A stopped case is never extended by
    plotting an analytical trajectory as if it were computed ALE data.

11. The u=0.6 strict temporal control (rtol=1e-9) failed at internal t=2.5201 s;
    the additional rtol=1e-10 control failed at 7.26963 s. Accepted output ends
    at 2 and 5 s, respectively. Both are numerical failures, not physical
    depletion events. They prevent a final long-time temporal-convergence
    claim for this case, even though the baseline and all spatial refinements
    reach 50 h. The u=0 and x008 stricter controls complete their full intervals.
    No unique floating-point/stiffness root cause has been established. The
    next numerical step should audit nondimensional scaling and interface
    Jacobians before further production; the present study reports NO-GO for
    final quantitative validation instead of treating failed controls as passes.

12. Early instantaneous contact-diffusion fluxes are not converged: x008 at1 s
    changes from+1.58e-5 (128 cells) to-1.06e-4 (256 cells) kg/m2/s, whereas the
    stricter128-cell integration gives+2.35e-6. Do not interpret these spikes
    as physical reversed diffusion. Full histories are preserved; Figure8
    separates startup and post-1 h fluxes. Windowed flux differences are in
    results/flux_resolution.csv. Small inventory/jump residuals do not certify
    instantaneous flux accuracy.
