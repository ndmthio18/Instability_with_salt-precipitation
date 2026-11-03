# WP6 first validated stages

Execution authority: `command/7.md`; frozen parent plan: commit `314d75d9`.
Completed: isolated seven-plane diagnostics, two threshold fixtures, T bridge,
and qualified 2T M0/M1 pair. **Stopped at the decision gate; no 4T or broader runs.**

- [Scientific findings](RESULTS_ZH.md)
- [Verification and actual costs](VERIFICATION_AND_COST_ZH.md)
- [GO/NO-GO decision](DECISION_GATE_ZH.md)
- [Multiscale validity design](MULTISCALE_VALIDITY_DESIGN_ZH.md)
- [Infrastructure, restoration and amendment](INFRASTRUCTURE_AND_AMENDMENT_ZH.md)

Figures 1–3 are available as PNG/PDF. CSVs retain signed fluxes, initial-inventory
correction, relative excess transport, and separate closure diagnostics. The
background endpoint screen is flagged, even though core numerical checks pass.
Unknown discretization errors are not zero; t_eff is undefined.

Read-only verification: `python scripts/prepare_wp6_plan.py --verify`.
Ignored-reference restoration: `python scripts/wp6_plan_integrity.py --restore-reference`.
Local analysis: `python scripts/analyze_wp6.py bridge` then `... extend`;
figures/results: `python scripts/report_wp6.py`;
cost report: `python scripts/report_wp6_cost.py`.
Restored ignored NPZ/field archives are needed for spatial figures. Frozen
reference products are never rewritten by restoration.

Run drivers are preserved for provenance, not as a request to rerun. Existing
start/launch claims prevent duplicate launches. Each new execution needs its own
authorization and registered scope. Raw archives are in ignored
`tmp/wp6_execution_20261103` and the dedicated frodo directory; manifests and
per-member hashes are committed.
