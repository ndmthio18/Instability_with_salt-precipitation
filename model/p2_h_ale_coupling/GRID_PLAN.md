# Registered diagnostic sequence (2026-11-08)

Base: 301097e3d2de3541bf63741cd729bd5cda22329a. New branch only.
Physics, donor algorithm, tolerances and 3.5 um threshold are frozen.

1. Instrument projection control without altering accepted solutions. Re-run
   historical precrust 40/80/160/320 and 320 coupled failure; retain every accepted
   cut-region budget and the first bound-rejected ALE trial spatial state.
2. Fixed physical cuts: 0.1, 0.5, 1, 2, 5, 10 mm. Inventories use exact overlap
   of rectangular SCV volumes with the region. Fluxes use the same residual
   kernel with difference of global BOX overlap fractions: a conservative
   piecewise BOX reconstruction, not nearest-node plane snapping. Record
   transport endpoint flux BEFORE the equilibrium projection and accepted
   inventory AFTER projection separately. Projection removed salt has its own
   regional ledger. Flux positive into the upper region.
3. Family A: fixed geometric mapping (total grading ratio fixed at the old
   80-layer mapping), 40/80/160/320/640 cells. Per-cell grading changes with N;
   geometry recorded from actual SCVs. Family B: upper 2 mm at fixed 25 um,
   deeper region 40/80/160/320 cells; independent surface refinement 50/25/12.5 um.
   This intentionally does not label the total cell count as the deep count.
4. Stage A-C is a decision gate before implementing any precipitation closure.
   Candidate if supported: equilibrium surface complementarity Ctop<=Ceq,
   Jprecip>=0, Jprecip*(Ceq-Ctop)=0, solved within the transport step; outgoing
   total liquid mass = evaporation + solid salt flux, outgoing salt = solid
   salt flux. Water boundary remains evaporation exactly. No volumetric cap.
   Old projection remains an explicit regression mode. No commitment that
   this candidate will be justified or converge.
5. ALE ablations after diagnosis: actual versus equilibrium initialization with
   excess salt in a declared external inactive solid ledger; fixed geometry
   and pore volume require the removed salt's mass to become water, so that
   option is forbidden without an explicit compensating water reservoir.
   Instead reconstruct equilibrium with conserved water and reduced pore
   volume, recording changed porosity/solid geometry, or use paired boundary
   tests without claiming isolated H1 causality. Prefer an exact conservative
   equilibrium reconstruction solving salt/water at fixed height.
6. Final gates: combined budget 1e-8, activation transfer 1e-12, adjacent grid
   scalar 2%, dt-half 1%; every failed case retained. Full coupling only after
   preactivation audit. Numerical and constitutive judgments remain separate.

All computation and analysis on HPCC Slurm, DUMUX_NUM_THREADS=1. Source snapshots
and compiler logs retained. No expensive work on login nodes.
