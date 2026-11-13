# Registered validation decisions

Reference numerical targets are <1% of nonzero vref=Sref/(rho_cr phi_core) for velocity differences, <1% geometric relative differences for positions/thickness, and <1% Sref for both interface diffusion fluxes. Absolute errors are always retained. Compare all matching saved times within each common requested interval; do not select a favorable endpoint. Conservation target is 2e-9 kg m−2 for each global/local inventory residual, separately from convergence. Tests are deterministic; no statistical confidence intervals apply.

P0: u=.6/.9 50h and x008 44.70875h at n=128, rtol1e-9/1e-10. P1: both .6 and x008, n128/256/512, liquid n/4, b0=.1/1/10nm, uniform and stationary-concentration initialization, tighter tolerance plus halved early internal-step caps. P2: .6/.9, kappa .5/1/2, n256, same physical duration. Long-time grid audit: n128/256/512 for .6/.9/x008. All improved matched cases use pinned setup_v6 with FRESH_JAC=setup. Earlier numerical ablations are not final comparisons.

P5 is permitted only after P0 and P1 pass and the P2 balance/sensitivity audit is addressed. It adds u0/.3 at matched n256/512 and any required further resolution. Kinetic sensitivity need not be zero to run a conditional comparison; a sensitive result prevents an unconditional physical-validation claim. The surface kinetic law lacks independent experimental calibration.

Analytical H1/H2 comparisons have no pre-existing empirical acceptance tolerance in Task28. Report full signed errors; do not retroactively choose a tolerance to call them validated. A 1% numerical target cannot be repurposed as an experimental accuracy guarantee. Overall scientific validation remains NO-GO if essential numerical convergence fails or bottom kinetics makes the claimed universal geometric law unsupported. Separate a resolved numerical gate from model-validity status.

The analytical phi is recomputed only during postprocessing. It never initializes or advances the ALE. Inferential bottom speed is diagnostic only. No tuning of kappa or phi to reduce growth-law discrepancies is allowed.
