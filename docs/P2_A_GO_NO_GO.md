# P2-A GO/NO-GO decision

> **Update after P2-B (see `P2_B_EVOLVING_POROSITY.md`).** Blocker 1 below has
> been addressed as a *modelling* gap: `surface/evolving_crust.py` evolves the
> crust porosity from Prat equation (2) on a depth-resolved material grid, is
> conservative to 1e-15, and reproduces equation (28) at second order in the
> frozen-uniform limit. The physical conclusion is sharper than before, not
> softer: with constant specific area the crust seals after growing by only
> `beta * eps_0 / a_v`, independent of the evaporation flux, which is less
> than one crystal diameter in the paper's own cubic construction. Blocker 1
> is therefore replaced by **gate 1'** in the gate list at the end of this
> document. Blocker 2 and every other item below are unchanged, and the
> overall decision is unchanged.

## Decision

**GO** for accepting the P2-A surface module as a verified, conservative,
standalone specification of physical film evaporation plus Prat (2026) crust
growth, with the one-way diagnostic chain of item 7 implemented and checked.

**NO-GO** for coupling it to the Richards solver, and therefore **NO-GO** for
any density-instability run with salt precipitation.

The split is necessary because item 7 explicitly scopes this stage as a one-way
diagnostic, and because two independent blockers were found during the work.

## What is established

1. The heterogeneity gas-film formulation is reused unchanged in structure and
   driven only by physical inputs: temperature, ambient relative humidity, film
   thickness, vapour diffusion coefficient and atmospheric pressure. No
   calibration against a prescribed flux was performed (item 2).
2. No local evaporation weight is normalised. The total evaporation rate falls
   41.6 % between fresh water and halite equilibrium, and 5.1 % along the
   frozen P1-R2 Case 36877 surface-salinity history (item 4).
3. Before surface salt exists, the water activity is computed from the local
   Richards soil-surface NaCl concentration, and the resulting vapour
   concentration gradient gives the local evaporation flux (item 3). The
   saturated-NaCl activity is 0.75121, inside the declared sanity band.
4. The Prat law is implemented in a standalone module with the finite-`Da`
   expression (28), the large-`Da` limit (29)/(32), crust porosity, crust
   height, explicit time integration and water/salt inventories (item 5).
   Three independent routes to `v_s` agree to 1.2e-14, `Da_top^2/Da = Da*`
   holds, and the reduced profile satisfies its ODE and boundary conditions.
5. Growth is gated behind an explicitly declared seed/precipitation event. At
   `h = 0` every Prat quantity raises a model-domain exit instead of returning
   a number, and the solubility crossing is recorded as eligibility, separate
   from activation (item 6).
6. The soil-side fluxes are reconstructed from equations (6)-(7) and the
   reduced crust transport solution. `F_s,soil` is strictly positive whenever
   solid accumulates, and `C_top` is never assumed equal to the soil-top
   concentration (item 8).
7. Salt and water budgets close to 1e-15 relative across all four
   interface/law combinations, under time-step refinement, and across the
   activation event. Explicit Euler converges at order 1.00 and Heun at order
   2.00.
8. Option A is acceptable exactly when `C*_top <= 0.01/S` with
   `S = -d ln j/dC = 4.4267` at the declared film; the measured boundary
   (2.211e-3 acceptable, 2.314e-3 not) brackets the predicted 2.259e-3
   (item 9).
9. Insufficient supply and loss of connectivity are reported and can stop the
   model; they are never used to clip `v_s`. The growth velocity is bit-identical
   across all supply fractions tested (item 10).
10. 61 registered tests pass. No frozen source was modified, no P1 result was
    altered, and no simulation was run.

## What is not established, and why coupling is held

### Blocker 1 - the constant-porosity assumption fails at the declared parameters

Prat equation (2) forces the crust porosity to change wherever `C != C_eq`, but
the reduced solution holds `eps_c` fixed. The module reports the
porosity-consumption timescale and flags the inconsistency. **At the declared
reference parameters the flag is false at 58 of 61 heights examined**, with
internal precipitation exceeding top growth by factors of 4 to 19.

This is not an implementation artefact. The derivation in
`P2_A_SURFACE_MODULE.md` section 4.2 gives, in the large-`Da` limit,

```
internal / top = 1 / sqrt(Da*) ,     f_inf = 1 / ( 1 + internal/top )
```

so a self-consistent constant-porosity crust requires `Da*` of order one, which
is exactly the regime where `f_inf -> 1` and the efficiency factor carries no
information. Prat (2026) states that `f_inf` stays below one over the relevant
parameter range, i.e. the relevant range is the one where internal
precipitation dominates and porosity cannot be constant. This tension is a
property of the published closure applied as a fixed-porosity subgrid, and it
must be resolved - by evolving `eps_c`, by a defensible `a_v`/`k_r`
measurement, or by a declared restriction of the usable window - before the
crust is allowed to feed a soil simulation.

### Blocker 2 - no C++ implementation

The module is Python. No C++ compiler was available in this environment, so the
DuMuX surface subgrid has not been written or compiled. The Python module is a
verified reference specification; it is not the coupled code.

### Further items that are not established

1. `a_v`, `k_r`, `eps_c`, the seed height and the film thickness are declared
   inputs. The audit supplies no value for `a_v` or `k_r`. Every absolute
   growth rate is parameter-conditional.
2. The declared film gives 5.82 times the legacy P1 prescribed Mode-A flux.
   Reproducing the legacy flux would need a 4.83e-2 m gas film. Which
   atmospheric control is right for Case 36877 is unresolved and is gate 1 of
   the P1-R2 GO/NO-GO.
3. No crust forms on any committed P1-R2 trajectory: the soil surface reaches
   at most 23.8 % of halite equilibrium in the frozen window. The growth part
   of the chain was exercised only on a declared synthetic concentration ramp.
4. The chain is one-way. The film flux was never fed back into the Richards
   solution, and the two sides use different atmospheric boundary conditions.
5. Lateral crust transport is zero by the independent-patch approximation,
   which is weakest exactly where horizontal variability is the object of
   study.
6. The activity correlation is a 298.15 K fit used at 293.15 K, inherited with
   the upstream file's own caveat. Nothing here validates it.
7. Numerical conservation of a subgrid does not validate the Prat closure, the
   film closure, the activation rule, or any onset prediction.

## Evidence-conclusion alignment

Supported:

> P2-A supplies a reproducible, conservative, parameter-transparent surface
> subgrid that reuses the heterogeneity film evaporation physics without
> calibration or normalisation, implements the Prat (2026) growth law with its
> finite- and large-`Da` forms, gates growth behind a declared seed event, and
> reconstructs soil-side water and NaCl fluxes whose budgets close to round-off.
> It also quantifies, with a closed first-order criterion, when the soil-top
> evaporation-surface approximation may be used instead of the Prat crust-top
> concentration.

Not supported, and must not be claimed:

> that the Prat closure is valid at any physical parameter set for this soil,
> that the declared film thickness represents the Case 36877 atmosphere, that
> crust growth occurs at a predicted time, or that the coupled model conserves,
> since no coupled model exists yet.

## Gates before coupling

1. ~~Resolve the constant-porosity tension.~~ **Done in P2-B.** Porosity now
   evolves from equation (2); see `P2_B_EVOLVING_POROSITY.md`.

1'. **(replaces gate 1)** Supply a defensible `a_v(eps)` closure, or a
   measurement of `a_v`, that permits a crust to reach a thickness worth
   modelling. With constant `a_v` the resolved model seals after
   `dh = beta eps_0 / a_v`, which is 8.4 % of the seed height at the declared
   inputs and under one crystal diameter in the paper's own cubic
   construction. P2-B shows the growth before sealing is
   `(beta eps_0^m / a_v0) * integral eps^-m d eps` for `a_v ~ eps^m` and is
   independent of the evaporation flux for every `m`, so no atmospheric or
   kinetic parameter can rescue it: only crystal coarsening can. That closure
   is the audit's own unresolved "changing crystal size" item and P2-B does
   not supply it. Until it exists, a crust that grows appreciably cannot be
   claimed, and neither can one that seals, since the sealing threshold is
   declared rather than derived from a percolation or permeability criterion.
2. Obtain or declare defensible `a_v` and `k_r` with provenance, and re-run the
   regime tables.
3. Decide the atmospheric control variable and film thickness for Case 36877,
   closing gate 1 of the P1-R2 GO/NO-GO.
4. Port the module to C++ against the P1-R2 isolated source, preserving the
   interface contract of `P1_R2_IMPLEMENTATION.md` section "Diagnostics and
   future interface contract", and reproduce every identity test in C++.
5. Specify the top-patch sign convention and the time discretisation of the
   soil/crust exchange, then demonstrate that the coupled budgets close on at
   least one 1D column before any 2D run.
6. Only then revisit the 2D density-instability question. Item 10's instruction
   stands: no large-scale density-instability run before conservative
   soil/crust coupling is verified.

Until gates 1', 2, 4 and 5 are met, coupling remains held.
