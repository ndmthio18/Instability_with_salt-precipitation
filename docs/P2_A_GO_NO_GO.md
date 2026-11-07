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

> **Update after P2-C (see `P2_B_SALT_GROW_PARAMETER_COMPARISON.md`).** The
> sealing conclusion above is a property of the **declared placeholder
> parameters**, not of the formulation. Transferring `eps_0 = 0.43`,
> `k_r = 2.3e-3 m/s` and `a_v = 820800 1/m` from the sibling `salt_grow`
> project moves `sqrt(Da*)` from 0.053 to 2.56, across the sealing threshold
> `(1-eps_0)/(eps_0-eps_stop) = 1.390`, and the crust then grows steadily
> instead of sealing. With every physical input and both closures matched, the
> P2-B quasi-steady evolving-porosity model reproduces the salt_grow transient
> ALE growth velocity to 0.05 %. Gate 1' is therefore **narrowed, not closed**:
> the parameters now have provenance but still no measurement, and sealing
> would return if `a_v` were 3.4 times larger or `k_r` 3.4 times smaller. Gate
> 2 is partly addressed for the same reason. Blocker 2 and the overall decision
> are unchanged. P2-C also found and fixed a defect in the P2-B **soil-side
> flux diagnostics** (not in any growth result): see the new item under "What
> is not established".

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

### Corrected during P2-C - the P2-B base-flux diagnostics

The committed P2-B `salt_flux_base_route` used a quasi-steady solution-mass
balance `d(rho U)/dz = -R` and evaluated the moving-interface pore volume at
the current top-cell porosity. Both are wrong: precipitation frees pore volume
that refills with solution, so `d(rho U)/dz = -R(1 - rho S/rho_cr)`, and the
new pore volume is created at the **deposition** porosity. This is the bulk of
the 7.5-10.4 % transport-versus-inventory route difference reported in the P2-B
verification, and the committed values of that specific diagnostic should be
read as a closure defect rather than as a physical discrepancy. It is fixed in
`graded_crust.py`, which drops the gap by three to four orders of magnitude; the
remainder is a closed-form residual of Prat's own linearised sink. **No growth
velocity, porosity, height or budget in P2-A or P2-B is affected**, because the
diagnostic is never fed back into the solution, and no committed file was
edited.

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

1'. **(replaces gate 1; narrowed by P2-C)** Supply a defensible `a_v` and
   `k_r`, measured rather than declared, together with a crystal-coarsening
   closure `a_v(eps)` if one is needed. P2-B showed that at the P2-A declared
   placeholders the crust seals after `dh = beta eps_0 / a_v`, 8.4 % of the
   seed height. P2-C showed that this is a statement about those placeholders:
   `k_r = 1e-7 m/s` is 23000 times smaller than the value salt_grow uses, and
   with the salt_grow values the crust does not seal at all. The gate is
   therefore no longer "find a coarsening closure or admit the crust cannot
   grow"; it is:

   - whether a crust seals is decided by the single criterion
     `sqrt(Da*) <= (1-eps_0)/(eps_0-eps_stop)`, which P2-C derived and
     confirmed on four of four cases, so the needed evidence is a measured
     `Da*`, not a trajectory;
   - the salt_grow values are a sibling project's declared baseline with
     provenance, not a measurement, and that project's own RESULTS.md lists
     independent kinetic measurements as still required;
   - the margin is 3.4 in `a_v` or `k_r`, so a measurement an order of
     magnitude away from the salt_grow baseline can put the model back on the
     sealing side;
   - `eps_stop` is still a declared porosity threshold rather than a
     percolation or permeability criterion, so a sealing *time* still cannot
     be claimed even when the criterion says sealing occurs;
   - the two projects disagree by a factor of 6 on the cubic-construction
     specific area (`36(1-eps)/d_sub` versus `6(1-eps_0)/d`), and the audit
     transcribes neither, so `a_v` is not traceable to the paper through this
     repository at all. The discrepancy errs away from sealing and does not
     change the conclusion, but it must be resolved before any `a_v` is
     quoted as Prat's.
2. Obtain or declare defensible `a_v` and `k_r` with provenance, and re-run the
   regime tables. **Partly addressed by P2-C**: the salt_grow baseline supplies
   values with traceable provenance and the regime tables have been re-run at
   them (`P2_B_SALT_GROW_PARAMETER_COMPARISON.md` sections 3, 5 and 6). What
   remains is a measurement rather than a declaration.
3. Decide the atmospheric control variable and film thickness for Case 36877,
   closing gate 1 of the P1-R2 GO/NO-GO.
4. Port the module to C++ against the P1-R2 isolated source, preserving the
   interface contract of `P1_R2_IMPLEMENTATION.md` section "Diagnostics and
   future interface contract", and reproduce every identity test in C++. The
   port must carry the corrected base-flux closure, the connectivity stop and
   a top-refined grid: P2-C shows that at realistic parameters the reaction
   layer is about 0.5 micrometres, so a surface subgrid using the committed
   P2-B uniform 32-cell discretisation would be under-resolved by a factor of
   about 23.
5. Specify the top-patch sign convention and the time discretisation of the
   soil/crust exchange, then demonstrate that the coupled budgets close on at
   least one 1D column before any 2D run.
6. Only then revisit the 2D density-instability question. Item 10's instruction
   stands: no large-scale density-instability run before conservative
   soil/crust coupling is verified.

Until gates 1', 2, 4 and 5 are met, coupling remains held. P2-C strengthens
the case that the *formulation* is sound -- it reproduces an independent
transient ALE implementation to 0.05 % on five diagnostics -- without changing
the decision, because agreement between two models sharing one closure is
verification and not validation, and because both projects record that the
reaction layer is thinner than the estimated crystal size.
