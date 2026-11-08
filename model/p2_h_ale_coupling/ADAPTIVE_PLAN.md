# Adaptive donor experiment (registered before execution)

Base 25a41470aac2fa52ae772d313f9ae2b35a49339f, branch research/p2-h-adaptive-donor.
Frozen ALE operator, BDF, film, soil constitutive laws, solubility and saturation
screen remain unchanged. Reference runtime and M3 evidence are inherited with
hash verification; changed application receives a real HPCC build/regression.

## Declared choices

- Default adaptive donor depth 5 mm; factorial sensitivity 2, 5, 10 mm.
- Eligibility requires the complete aggregated BOX to lie within the depth.
  Depth uses actual rectangular tensor SCV extents, not uniform spacing.
- Aggregate all SCVs once per global DOF; sort by DOF depth. At equal depth
  share extraction proportionally to available mass to preserve 1D symmetry.
- Retain Sw >= .99; capacity uses .99+1e-13 for inversion roundoff headroom.
  Density is each donor's existing composition-dependent soil density. Salt
  and water are removed in that donor's existing ratio and pressure inverted.
- Mix transferred liquid uniformly in initial ALE pores by mass; ALE density
  remains its original constant saturated-brine value. No ALE equation changes.
- The old top_only algorithm is unchanged. Adaptive is explicitly selected.
- Coupled horizon 2500 s, dt .5 s, 64 crust cells, original graded soil input
  (Grading1=-.97), 0.2 m depth, 1 mm width and one horizontal cell.

## Stages

1. Compile donor unit tests and actual DuMuX; smoke adaptive 80/160 to400 s.
2. Full old-mode regression at 40/80/160/320 cells and 1 um; extra old80 at
   2/5/10 um. Compare old80/160 against frozen full accepted histories.
3. Adaptive factorial: 4 soil grids x 4 thresholds x 3 depth bounds =48 runs.
   Keep every outcome, including refused activation and stopped prefixes.
4. Controlled allocation at identical prescribed h=1/2/5/10 um, all four
   grids and three depth bounds, using the same analytic initial soil profile.
   This is a separate allocation-only harness: no artificial solid enters
   production evolution. Add uniform-mesh 80/160 controls for geometry.
5. PrecrustOnly 40/80/160/320: retain original precipitation to2500 s with
   activation disabled, separating threshold-crossing discretization.
6. dt=.25 s at160/320, h1um, depth5mm. Activation and active-step injected
   rollback at80, plus depth1nm refusal harness.
7. Before interpreting the matrix, supplement the controlled refusal with a
   real 400-s n80 transport pair: adaptive depth1nm versus PrecrustOnly. Their
   full accepted state histories must match exactly despite repeated attempts.

## Fixed gates and comparisons

Combined salt and water residual <=1e-8 of their initial total inventories.
Activation transfer residual <=1e-12 of corresponding initial inventory.
Compare h, solid-equivalent height (total solid), evaporation at exact common
physical times; successive40/80,80/160,160/320 differences <=2%, denominator
finer result. No activation, unequal endpoints or missing data fail the final
2500-s gate. Also retain common accepted prefix comparisons, clearly labeled.
Use 1% for the exploratory timestep gate. Inspect profile and trend evidence;
passing scalar gates is not asymptotic convergence or capillary validation.
No automatic adjustment of depth, Sw, threshold, physics, mesh or tolerances
after a failure. Event logs distinguish attempts, trials and accepted activation.
First-refusal and every successful-trial profiles retain each donor's volume,
depth, capacity, actual transferred mass, pressure, saturation and composition.

Supplementary descriptive analysis added during matrix execution (after some
40/80-cell results were visible): compare exact common *crust ages* at depth5mm
to examine the contribution of activation-time offsets. These are different
physical times, explicitly reported, and NEVER replace the registered
common-physical-time convergence gates. This is an error-attribution diagnostic,
not an additional preregistered acceptance test or a causal proof.

Operational amendment during the healthy matrix: the cluster denied extending
running job16371535 beyond its original60-minute allocation. A dependent
completion job16371664 (4 CPUs,90-minute allocation) checks all commands against
recorded originals. It is a no-op if the primary completed. Only a primary
Slurm TIMEOUT permits archiving and replaying interrupted/unstarted cases from
the same initial conditions; numerical failures are not automatically retried.
Each case retains the original1800-s wall timeout and all scientific settings.
Audit job16371564 waits for this completion check. This changes scheduling only.
