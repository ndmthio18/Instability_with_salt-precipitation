# P2-H dual-moving-front verification artifacts

This directory contains the lightweight, reproducible outputs of the
prescribed-rate diagnostic prototype. Full per-step histories are generated
locally but ignored by Git; the committed summary tables contain the A--G test
endpoints and grid/time-step comparisons.

Baseline commit: `fdd4d1f` (`research/p2-h-flux-mechanism-x009`).

Reproduce from the repository root:

```powershell
python -m unittest model.p2_h_dual_fronts.test_dual_front -v
python model\p2_h_dual_fronts\run_verification.py
```

The result status is `PASS_DIAGNOSTIC_NOT_PREDICTIVE`. The bottom retreat rate
and dissolved-product partition are prescribed inputs, not a soil-derived
physical closure.

Figure semantics and accessibility notes:

- All curves are deterministic model evaluations; there are no replicates or
  statistical uncertainty intervals.
- Normalization is stated on each axis. Log axes contain strictly positive
  values; the signed velocity-ratio figure remains linear and shows the zero
  crossing explicitly.
- Color is redundant with line style or marker. PNG files are opaque previews;
  PDF companions preserve vector lines and text.
- The CSV tables are the accessible data alternative. In brief: increasing
  `Da` exponentially weakens the direct `u` sensitivity, while the dual-front
  panels show pure growth, pure retreat, constant-thickness translation, and
  upward translation with thinning.
- No target journal or final print width was specified, so these are
  provisional general scientific figures rather than publisher-compliant
  submission files.
