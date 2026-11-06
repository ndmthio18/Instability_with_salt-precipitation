# P1-R2 isolated conservative Richards-salt model

This directory is the isolated implementation requested by `task/3.md`.

- `source/`: buildable DuMuX source; no frozen source is modified.
- `source_diff_from_wp7_standard_2d.patch`: exact patch against the selected frozen standard-2D Newton source.
- `verification/inputs/`: V1-V3, V5, V7 and the single authorized V8 input.
- `verification/`: unit tests and deterministic post-processing scripts.
- `results/`: lightweight accepted-step, budget, comparison and figure products; raw VTU/PVD fields are excluded.
- `provenance/`: source/build/archive hashes and dependency revisions.

Local pure-Python checks:

```powershell
python analysis\p1_r1_unit_audit_20261106\test_conversion_audit.py -v
python model\p1_r2_conservative_richards\verification\test_physics.py -v
```

The build requires the pinned DuMuX/Dune stack recorded in `provenance/source_manifest.json`. All scientific inputs must explicitly provide the two molar masses and three density-law constants; rounded fallbacks are intentionally unavailable.

Interpretation is bounded: the model is discretely conservative under a declared mixture-mass-average closure. That is not equivalent to physical validation or a converged Prat-coupled prediction.
