# P2-E — DuMuX–Prat 1D coupling validation

## Decision

**NO-GO for the first production Richards–Prat coupling.**

The C++ lifecycle defect is fixed and the frozen Python-to-C++ transfer is
verified, but the coupled surface layer and seed perturbation are not
grid-converged. The 160-layer case leaves the Prat validity envelope before
one accepted step, and the independent crust transport and inventory routes
miss the preregistered `1e-6` relative gate. No 2D instability run is
authorized by this result.

## What passed

- The frozen `pratcrust.hh` compiled with GCC 15.2.1 on `sal`; all 64 reference
  vectors passed. The worst relative difference was `9.712556e-14` against a
  `1e-11` gate.
- `pratcoupling.hh` separates pure trial evaluation from commit. Repeat
  evaluation, rollback, reduced-time-step retry, stale-trial rejection,
  double-commit rejection, and forced root failure all passed as compiled C++.
- Gas-film and crust-top subroots are rechecked after return. The lifecycle
  test residuals were `-2.03e-20` and `-4.24e-22 kg m-2 s-1`.
- The patched P1-R2 DuMuX executable compiled in the pinned DuMuX/DUNE stack.
  `Problem::neumann` only evaluates a trial from the current Newton iterate;
  the persistent crust is replayed and committed once after DuMuX accepts the
  step. A bounded 80-layer case completed 20 accepted one-second steps.
- Seed formation is an explicit soil-to-crust transaction. For seed heights
  5, 10, and 20 um, combined salt and water residuals were zero to floating
  precision (worst salt residual `3.47e-18 kg`, water `2.78e-17 kg`). The top
  boundary degree of freedom remains at the converged saturation composition;
  the debit is taken from the nearest underlying control volumes.
- Closed-bottom and reservoir-bottom configurations are separate. The
  reservoir run records nonzero bottom component influx; the closed-bottom
  run records zero bottom influx.

## Failed or incomplete gates

1. **Surface refinement:** after the conservative seed transfer, the 40- and
   80-layer cases give opposite-sign salt inventory fluxes at 20 s
   (`+2.04e-6` versus `-3.73e-5 kg m-2 s-1`). The 160-layer Newton path exits
   because the CVFE boundary state becomes more than 2% undersaturated. It has
   zero accepted steps. Top salinity and evaporation agreement between the two
   completed grids therefore cannot establish convergence.
2. **Solubility-crossing time:** the bounded tests start from a converged
   saturated interface to isolate activation and lifecycle behavior. A
   pre-activation activity-film refinement study was not completed, so the
   crossing-time gate remains failed rather than being assigned the trivial
   value zero.
3. **Three-route salt closure:** the initially active test already gives
   `4.84e-6` relative transport-versus-inventory mismatch. With conservative
   activation it is `5.91e-5` on 40 layers and `3.05e-5` on 80 layers. These
   exceed the preregistered `1e-6` gate and are consistent with the bounded
   nonzero route mismatch already reported by P2-D.
4. **Accepted-step budgets:** the 40-layer activated run closes total and salt
   inventories to `1.11e-9` and `9.34e-9` relative, but the 80-layer run gives
   `1.52e-8` and `1.09e-7`; the required `1e-8` gate is not grid-independent.
5. **Seed sensitivity:** the first accepted salt inventory flux changes from
   `+9.32e-7` to `-5.01e-6` to `-1.63e-5 kg m-2 s-1` as seed height changes
   from 5 to 10 to 20 um. This is a strong, unresolved seed perturbation, not a
   nucleation-rate model.
6. **Hydraulic supply:** the model still lacks equality of soil and crust
   hydraulic potentials. The only defensible numerical supply envelope is the
   bounded observation that 40/80-layer cases can sustain 20 s while the
   160-layer case sustains zero accepted time. This is neither mesh-converged
   nor a physical validation of continuously saturated crust supply. The code
   refuses states outside the declared concentration envelope instead of
   forcing the crust demand through them.

## Reproduction boundary

The integration is an overlay patch series in `integration/`; frozen P1-R2 and
P2-A/B/C/D files are unchanged. Raw VTU and detailed run logs remain on `sal`
under `/temp-hd/zhixin/p2_e_20261107/`. Git contains only source, criteria,
small result summaries, hashes, and the decision.
