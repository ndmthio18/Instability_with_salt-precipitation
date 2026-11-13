# Numerical reliability audit — Task 28

## Frozen failure reproduction

Baseline is commit `04f9bb2d856512531ebcb6265e4b46ef35cc43b0`. The `source_legacy` reproduction changes diagnostics only, not parameters, equations or tolerances. The baseline u=0.6, 128/32-cell, 1 nm liquid-layer failures were reproduced at internal times 2.520100517189 s (rtol 1e-9) and 7.269625167084 s (rtol 1e-10). Last requested outputs were 2 and 5 s. They are not the exact failure times.

For rtol 1e-9 the proposed step fell to 2.39558e-6 s from a previous accepted 0.125792 s; weighted local errors were largest in liquid salt inventory indices 268 and 256. At rtol 1e-10, indices 282 and 275 were largest, with proposed step 3.11105e-6 s following 0.162959 s. Pore/solid/liquid volumes remained positive; no depletion or concentration-bound crossing explains these failures. Logs preserve nonlinear iterations, error failures, Jacobian counts and failure-state weights.

## Exact affine conservation-variable transformation

For pore volume per unit area E and dissolved NaCl inventory Q, define

\[
T=Q-\rho_l C_{ref}E,\qquad C=C_{ref}+T/(\rho_l E).
\]

With liquid total-mass flux L and salt flux S in the moving frame, define U=S−Cref L. Then

\[
\dot T_i=U_{i-1/2}-U_{i+1/2}-(1-C_{ref})R_iV_i.
\]

In the liquid-only region the reaction term is zero. Reconstruction Q=T+rho Cref E preserves the original salt and water equations exactly. Production Cref=Ceq; the constant-bath verification uses its bath concentration. This is not a new physical approximation or a fitted concentration.

Physical equivalence was checked against the original full-Q right-hand side at a nonuniform perturbed state (normalized maximum error 8.94e-11). A forward automatic-differentiation Jacobian passed an independent centered directional-difference comparison (normalized error 2.41e-10). The Jacobian includes nonlinear top and bottom root derivatives and geometry-dependent coefficients. It is not inferred from the analytical growth law.

## Scaling and error control

The dimensionless integration state is y/scale. Solid-volume scales are initial cell volumes, T scales are rho times initial cell volume times 1e-5, geometry scales are initial crust thickness, and ledger scales are rho_cr times initial thickness. Production material values are unchanged.

T error denominator is rtol[rho E (1e-5)+|T|]. Therefore the concentration perturbation is resolved on its own scale instead of relative to the O(0.265) equilibrium background. The value 1e-5 is an explicit numerical reference scale, not an absolute physical concentration floor. Both 1e-9 and 1e-10 tolerances are tested. Solid, geometry and ledger weights are documented in `two_bdf.hh`. Accepted states are not projected or clipped; the previous Newton-trial domain extension and physical stop screens are retained. T can have either sign, as required for dissolution.

The finite-difference ablation uses scaled centered differences, with one-sided differences when solid or bottom-position constraints would be crossed. The default uses automatic derivatives. Long-double intermediate arithmetic is an additional numerical experiment; CVODE state and matrix arithmetic remain double precision.

## Controlled evidence and limits

Six short actual-material runs completed at both strict tolerances. Reverting to full-Q relative error weights caused another failure. Centered finite-difference and unscaled-state ablations completed their short windows. Thus improved perturbation error control matters in this configuration, but these experiments do not prove that AD or scaling is individually necessary.

Six long-double P0 actual-material runs subsequently reached 50 h for u=0.6/0.9 and 44.70875 h for x008. Completion is separated from the matched observable differences in `results/temporal_convergence.csv`.

The high-D constant-bath dissolution stress test failed with both double and extended intermediates. It is not hidden, relabelled as a production case, or repaired by changing material values/tolerances. Updating the Jacobian every step together with tighter nonlinear convergence passed this same stress case. This motivates a separate nonlinear-convergence ablation; all versions and failed prefixes are retained. Success of a combined intervention does not identify a unique root cause.

Diagnostics report error weights/local estimates at the current internal accepted state. At a failed advance these are the solver's available diagnostics, not a saved rejected trial vector. Requested output times and internal times are stored separately. Large individual weighted components may coexist with a passed weighted RMS norm.


## Jacobian-age / nonlinear-solve ablation

The tighter nonlinear criterion alone (`FRESH_JAC=nonlinear`, coefficient .01 rather than .1, maximum 8 iterations rather than 3) completed the thin-layer early case but still failed the same pure-dissolution stress test at 0.0225947 s. Updating the Jacobian at every linear-solver setup (`FRESH_JAC=setup`, `CVodeSetJacEvalFrequency(1)`, ordinary setup frequency retained) together with those tighter nonlinear settings passed the original stress test and the formerly failing 0.1 nm case. Updating every step also passed, but is unnecessarily expensive for the verified cases. No error-test acceptance threshold was loosened and no accepted state was clipped.

This supports stale linearization / insufficient nonlinear convergence as additional contributors, beyond concentration-perturbation error scaling. The evidence identifies a working combined remedy within the tested envelope, not a universal unique cause. The `setup_v6` executable is frozen and used for all final matched cases, with `FRESH_JAC=setup` explicitly exported in every run script.

The stress case's constant-bath recession discrepancy was about 5.83e-5, below its pre-existing 1e-4 asymptotic tolerance. All tests, including pure dissolution, were run in the successful setup/fresh verification (no `SKIP_STRESS`). The earlier `SKIP_STRESS=1` run tested remaining independent checks only and was never a complete physical gate. Its log wording and environment are preserved.


## Arithmetic-only optimization and replay

Version sqrt_v7 replaces pow(eps,1.5) with eps*sqrt(eps), and its derivative with 1.5*sqrt(eps). This eliminates repeated general-power evaluations during forward automatic differentiation. The physical equations, state/error scales, nonlinear criteria, tolerances and bounds are unchanged. The entire original verification suite, all six P0 trajectories and all 30 P1 cases were repeated successfully. Final matched comparisons use this frozen executable throughout. `results/implementation_bridge.csv` compares its complete P0 histories with setup_v6 at the same settings; it is not a physical sensitivity experiment.


## Optional parallel Jacobian experiment

The independent Jacobian columns were evaluated with eight OpenMP threads in a separate v8 performance experiment. The full mathematical verification passed; 1-thread, 8-thread and v7 saved histories were byte-identical over the registered 128-cell x008 1000 s test. At 512 cells, 126 common saved rows through 7020 s were also byte-identical. This experiment did not show an observed throughput benefit on its allocated node, so its duplicate long run was intentionally censored, with Slurm cancellation and accepted prefix retained. The primary results remain the uniform serial v7 implementation. This is not a numerical failure and is not included among complete primary runs.


## Baseline finite-difference and active-set audit

The original Jacobian already used physically dimensioned perturbations: sqrt(machine_epsilon) times solid cell volume, equilibrium dissolved inventory, or geometric length; it used centered differences except near nonnegative solid/bottom-position constraints. Thus the evidence does not justify saying it had no physical scaling. Its dissolved-inventory perturbation corresponds to concentration of order sqrt(machine_epsilon)*Ceq, about 4e-9, and its relative error weight resolves the large equilibrium background rather than the small excess. The new FD ablation uses scaled differences; automatic differentiation removes that perturbation choice but is not shown to be individually necessary.

At exactly Cb=Ceq, max(Ceq−Cb,0) has different one-sided derivatives. The chosen inactive-set derivative is not a unique classical derivative. `verify_interface.cc` independently checks the active and inactive sides against the implicit jump-root derivative: dCb/de=(gc+gl)/(gc+gl+A(1−Ceq)) on the dissolving side and dvb/de=−(rho_l/rho_cr)kr*kappa*dCb/de; the non-dissolving derivative is zero. This supplements the directional Jacobian check away from the kink.

Diagnostic indices are: [0,n) solid volume; [n,2n) crust dissolved inventory (legacy Q, new T); [2n,2n+nl) liquid dissolved inventory; then top position, bottom position and the 11 integral ledgers. At n=128,nl=32 the troublesome 256–287 block is the liquid region, not an exhausted crystal cell.
