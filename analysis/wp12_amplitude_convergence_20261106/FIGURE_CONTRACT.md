# Figure contract

## Results-level question

How do the intermediate `1e-5` and enhanced `1e-4` Case 36877 eigenprofile
amplitudes change nonlinear timing, R24/R48/R96 agreement and the resolved
transport pathway relative to `1e-6`?

## Figure-level claims and evidence roles

1. **Timing and alignment:** physical-time trajectories establish acceleration;
   the separately shown preregistered `t-t_nl` coordinate tests whether later
   differences persist after a fixed nonlinear-state alignment. No fitted time
   shift is allowed.
2. **Transport convergence:** conservative SCVF peak, positive integral and
   signed integral at 1, 2, 5 and 10 mm provide the decisive grid-consistency
   comparison. Same-grid amplitude contrasts and same-amplitude grid contrasts
   must occupy distinct panels/tables.
3. **Mechanism and numerical symmetry:** modal energy, resolved local downflow,
   odd-mode energy and mirror error determine whether the enhanced perturbation
   changes the physical pathway or only its timing and symmetry-breaking
   exposure.
4. **Validity boundary:** salt/water closure, stopping reason, executable/input
   hashes and cost records bound every scientific conclusion.

## Visual contract

- Archetype: quantitative grid with a decisive transport panel and subordinate
   timing/mechanism panels. The `1e-5` and `1e-4` grid-sensitivity panels remain
   separate so that amplitude and grid effects are not pooled.
- Backend: Python/Matplotlib only.
- Export: 180 mm wide PDF and 300 dpi PNG; editable PDF text; minimum rendered
  glyph size 5 pt.
- Encoding: grid by line style/marker, amplitude by restrained colour family;
  each fixed-amplitude grid comparison uses one consistent three-grid mapping.
- Statistics: these are deterministic single trajectories, not independent
  replicates. No confidence interval or hypothesis-test annotation is valid.
- Source data: only generated CSV/JSON reductions from the nine registered
  trajectories.
- QA: render-time panel-alignment JSON, PDF text audit, collision audit and
  panel-by-panel final-size inspection are blocking gates.
