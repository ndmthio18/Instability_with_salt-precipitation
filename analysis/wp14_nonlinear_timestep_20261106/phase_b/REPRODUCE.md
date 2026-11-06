# WP14 Phase B reproduction

## Frozen execution

The authorized input and execution plan are in `inputs/`. The sole production simulation is Case 36877 on R48 x 200 with `A=1e-4`, the original A28 profile, `dt=56.842943238478475 s`, `OutputInterval=100`, and a maximum horizon of 6T. The `plan.json` and `input_diff_audit.json` files define the only allowed changes relative to the frozen WP12 R48 input.

The remote production, acceptance and reduction jobs were executed on the approved hpcc compute partition. Their immutable records, scheduler accounting, logs, hashes and stopping event are copied to `provenance/`. Raw product archives in `raw/` are retained locally for audit but intentionally excluded from Git.

The extracted 69.6 MB per-step `transport_budget.csv` and redundant 2.8 MB `surface_events.csv` are likewise retained locally and ignored. Their SHA-verified source is `raw/refined_R48_dt2_metrics.tar.gz`; the committed scalar/field metrics and matched comparison tables are the lightweight scientific deliverable.

## Reduced analysis

From the `paper/2` directory:

```powershell
python scripts/wp14_phase_b_compare.py
python scripts/wp14_phase_b_figures.py
```

`wp14_phase_b_compare.py` pairs original step `k` with refined step `2k`. Field snapshots pair the original 50-step interval with the refined 100-step interval. It writes the matched scalar, modal, conservative-flux, inventory, profile and representative-field tables plus the finite-resolution error budget under `phase_b/comparison/`.

`wp14_phase_b_figures.py` regenerates Figures 5--7 as PDF, SVG, PNG and 600-dpi TIFF. TIFF is a transient source-validation artifact and is not committed. The figure audit uses the local `nature-figure` checks for source semantics, PDF font size, panel alignment and rendered collision detection.

## Independent acceptance

After generating the report and figure QA records:

```powershell
python scripts/verify_wp14_phase_b.py
```

The verifier rechecks input scope, frozen reference archive hashes, the independent archive gate, the refined archive hash, expected row counts, exact physical-time matching, conservation flags, required deliverables and figure QA. It does not reinterpret numerical acceptance as physical validation.
