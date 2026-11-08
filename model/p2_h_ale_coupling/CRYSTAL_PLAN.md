# P2-H crystal-size activation and boundary correction

Registered before runs, 2026-11-08. Base c1fb2e56ccee892206486283e5262cd41354a814;
branch research/p2-h-crystal-boundary. Do not modify old branch or frozen soil.

- Default activation is 3.5e-6 m, a declared estimated crystal-size criterion,
  not a physically validated nucleation law and not an ALE cell size.
- Production starts with zero height and zero solid. Actual precipitated
  inventory sets birth height, including threshold overshoot; epsilon_birth=.43.
- Adaptive transfer, mass-weighted mixing, Sw>=.99 and 5 mm maximum donor depth
  are unchanged. Adaptive mode is explicit; top_only remains available.
- Correct the isolated P2-H application's top Neumann face predicate using
  upward normal plus position. Original frozen soil source remains unchanged.
  Log every actual boundary SCVF and independently check its normal, coordinate,
  area and classification, including the 320-cell legacy lateral false positives.
- Same 0.2 m x .001 m domain, one horizontal cell, grading -.97, original
  initial conditions, gas film and solver tolerances. 64 ALE cells. TEnd2500 s,
  max dt .5 s. No ALE/BDF/Prat/kinetic/solubility change. Unused legacy P2F seed
  parameters are omitted from new input files, not read by P2-H.
- Build/units then controlled allocation at default3.5um on all four grids,
  320-cell real smoke to120 s to cross the former boundary failure time.
- Production: default3.5um adaptive40/80/160/320; precrust-only40/80/160/320;
  dt/2 at160/320; explicit1um top_only80/160 regression against old histories.
  Main runs omit ActivationHeight from CLI to test the actual default.
- Each case bounded to3600 wall s, batch limit150min, eight independent cases
  with one thread each. No automatic numerical retries; retain failure prefixes.
- Gates: global water/salt <=1e-8 relative initial inventory, activation transfer
  <=1e-12; final2500 s adjacent-grid h/total-solid/evap differences <=2%; dt gate1%.
  Inspect internal precipitation, porosity and soil profiles separately. Missing
  activation or final accepted state fails the corresponding final gate.
- Preserve separate event times (precipitation, attempt, accepted activation).
  Do not infer asymptotic convergence from this nonuniform refinement family.

## Recorded operational observation during execution

The intended one-thread process setting only constrained OMP/BLAS. Inspection
of actual processes in job16374135 showed 11 OS threads per solver; the pinned
DuMux uses TBB and reads DUMUX_NUM_THREADS independently. No running solver was
restarted or retuned. A separate120-s n320 check (job16378274) set
DUMUX_NUM_THREADS=1, used the identical binary/input/CLI, and reproduced the
entire accepted history byte-for-byte (178.43 wall seconds). This validates
that prefix only, not a full-horizon thread convergence claim. The production
matrix retains its original environment so its reproduction script is unchanged.
