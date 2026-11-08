# P2-H interface concentration audit: registered execution plan

The audit starts from commit `1fe9d1f7f9fee8094b671367992d463a8e8b96c9`
on the isolated branch `research/p2-h-interface-concentration-audit`.
Historical result branches are read-only references.

## Evidence order

1. Reconstruct the accepted-state concentration and flux histories from the
   committed 0.1 and 3.5 micrometre failures and the successful
   `x0=0.09999/no_flux` reference.
2. Run a matched four-case cross comparison at `k_r=1e-3 m/s`,
   `alpha_av=1`, 80 soil cells, 64 ALE cells and `dt=0.5 s`.
3. Run 64/128/256-cell ALE and one 160-cell soil diagnostic refinement.
4. Run standalone, source-accounted concentration counterfactuals for 500 s.
5. Independently rebuild the ALE interface operator and the soil top flux from
   accepted storage changes and the separately reconstructed reservoir influx.

## Predeclared interpretation rules

- Completion is not evidence of physical correctness.
- A fixed-`Ceq` diagnostic carries an explicit bottom-salt ledger.
- Equilibrium initialization transfers salt from a declared solid ledger and
  returns the displaced pore water to a soil-side ledger; it is not a silent
  overwrite and its fixed-geometry inconsistency is reported.
- Failed trajectories are compared only over common accepted physical times.
- No 50 h scan is authorized in this audit.
