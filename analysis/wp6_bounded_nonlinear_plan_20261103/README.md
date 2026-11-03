# WP6 planning package

Status: **plan submitted; no simulations launched**. This package is a proposal and preregistration, not a new numerical result.

- [运行计划](RUN_PLAN_ZH.md): canonical nx28 stages, seven conservative observation planes, stopping/censoring, symmetry controls, separate nx60 preparation.
- [预登记GO/NO-GO](GO_NO_GO_ZH.md): quality, signal-trigger, sensitivity, event timing and full-domain gates.
- [成本预算](COST_BUDGET_ZH.md): measured anchors, extrapolation assumptions, conditional caps and resource admission.
- `preregistration.json`, `stage_schedule.csv`, `cost_budget.csv`: machine-readable schedule, thresholds and budget.
- `observation_planes.csv`, `vertical_grid_*_vertices.csv`: seven fixed BOX faces and an odd-subdivision nested grid; solver integration remains pending.
- `source_hashes.json`, `plan_checks.json`, `artifact_manifest.json`: provenance and offline checks. Paper1 and WP5 sources/results unchanged.

Reproduce offline with `python scripts/prepare_wp6_plan.py` from the Paper2 root. Requires the already archived WP5 `tmp/single_wavelength_validation_20261103/reference/M1/post/profiles.npz`; it never invokes SSH, compilation, scheduling or a solver. Source profiles and input differences are audited, but no executable WP6 run inputs or new ELF are claimed ready.

After registration is accepted for execution, preserve this version. Record any protocol changes in a separate dated amendment before inspecting dependent new results.
