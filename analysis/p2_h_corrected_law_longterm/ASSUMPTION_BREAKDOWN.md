# Assumption audit

The long-term study exposes distinct limits that must not be combined into a
single analytical-versus-ALE error.

1. **Spatial resolution.** The original 64-cell, power-1.5 mesh develops a
   substantial numerical top-speed bias as h grows. At 50 h it predicts
   4.070 nm/s, whereas the refined 512-cell power-2 result is 3.507 nm/s.
   A small mass residual does not remove this error. The 256/512 power-2
   common-time maximum top-speed difference is 0.506%.
2. **Constant properties.** On the refined u=0 trajectory, original Model A
   underpredicts the 50 h top speed by 5.05% (ALE denominator). The mean-porosity
   substitute has 12.44% error and the upper-3-ell substitute has 7.31% error.
   A single porosity replacement therefore does not reliably improve the law.
3. **Spatially variable reaction/diffusion.** Freezing the actual porosity
   profile in a linear independent BVP reduces the same endpoint discrepancy
   to 7.31e-6 relative. This identifies nonuniform material coefficients as the
   principal missing feature for this low-Pe reference trajectory. Model C
   uses the same finite-volume grid, so agreement alone does not independently
   validate that grid; the transient 256/512 comparison supplies that check.
4. **Quasi-steady concentration versus changing inventory.** At that endpoint
   the reconstructed concentration-storage/reaction magnitude ratio is
   8.23e-8. Nonzero dissolved inventory growth does not by itself imply a large
   concentration transient: pore-volume growth contributes to Q storage.
   The complete conservative term decomposition remains available. These
   statements do not generalize automatically to the high-evaporation P1 runs.
5. **Finite skeleton.** With back diffusion, bottom material approaches zero
   solid fraction while the top continues precipitating. The fixed-bottom
   saturated-domain closure supplies no mechanical contact, gap or moving-bottom
   law. The new runs therefore stop at epsilon_max=.999. They establish finite
   buffering before that screen, not indefinite support of an intact crust.
6. **Event resolution.** The P0 screen crossing is bracketed by at most 60 s.
   The u=.3 bracket varies with mesh; report that uncertainty instead of claiming
   a converged mechanical lifetime. Additional P1 event checks use 1 s outputs
   for cases that crossed a solver state bound between outputs.
7. **Continuum interpretation.** h0=3.5 um and ell0=0.779 um are smaller than the
   nominal 25 um grain parameter. Numerical convergence establishes the solution
   of the continuum equations, not the validity of this constitutive model at
   subgrain scales.

No term in Models A--C was calibrated to the computed top velocity. Model C
conditions on an ALE porosity profile and is therefore a mechanistic diagnostic,
not a standalone predictive reduced model with its own porosity evolution law.

P1 extends the audit: at initial Pe=.31485 the endpoint top concentration
reaches .29851 (Ceq=.26494); with u=.6 the bottom concentration falls to .13490.
Hence the high-evaporation operation changes advection AND the small-deviation
condition. The frozen linear RD BVP's remaining roughly 6% discrepancy cannot
be assigned uniquely to advection. No independent advective BVP was completed.
At Da*=.19844, even 256/512 comparisons retain 1.44--2.35% speed differences
near pore filling; those four parameter combinations remain unresolved.
