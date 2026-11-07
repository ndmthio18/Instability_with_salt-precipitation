# P2-A isolated surface subgrid

Standalone external porous-salt surface module: heterogeneity gas-film
evaporation plus the Prat (2026) growth law. It does **not** modify P1-R2, any
frozen source, or any existing result, and it runs no DuMuX simulation.

```
surface/            the module
  thermo.py           source-traced thermophysical correlations
  water_activity.py   NaCl water activity (Pitzer), reused from heterogeneity
  film_evaporation.py physical gas-film evaporation flux
  prat_crust.py       P2-A: Prat closed growth law, frozen uniform porosity
  coupling.py         P2-A: one-way chain, base fluxes, interface options A/B
  evolving_crust.py   P2-B: depth-resolved crust with porosity from eq (2)
  graded_crust.py     P2-C: the same physics on a top-refined, optionally
                      coarsened grid, plus the salt_grow B1/B2 closures and
                      corrected soil-side flux diagnostics
verification/
  config/acceptance_criteria.json       P2-A tolerances and declared inputs
  config/acceptance_criteria_p2b.json   P2-B tolerances, inputs and the
                                        predictions registered before
                                        measurement
  config/acceptance_criteria_p2c.json   P2-C transferred parameters, declared
                                        differences from salt_grow and the
                                        predictions registered before
                                        measurement
  reference_config.py        the declared P2-A reference configuration
  reference_config_p2b.py    the declared P2-B reference configuration
  reference_config_p2c.py    the declared P2-C configuration and case ladder
  test_p2a_surface.py        61 tests, maps onto task-4 item 10
  test_p2b_evolving.py       32 tests, evolving porosity
  test_p2c_graded.py         44 tests, graded grid and parameter transfer
  run_verification.py        regenerates the P2-A evidence tables
  run_p2b.py                 regenerates the P2-B evidence tables
  run_p2c.py                 regenerates the P2-C evidence tables
  drive_from_p1r2.py         one-way chain driven by frozen P1-R2 data
results/            generated evidence (CSV + JSON)
```

`prat_crust.py` and `coupling.py` are the audited P2-A closed-form path and are
left unmodified by P2-B; `evolving_crust.py` must reproduce them in the
frozen-uniform-porosity limit, and a test enforces that.

Reproduce:

```bash
cd verification && python -m unittest test_p2a_surface test_p2b_evolving test_p2c_graded
```

```bash
cd verification && python run_verification.py
```

```bash
cd verification && python run_p2b.py
```

```bash
cd verification && python run_p2c.py
```

```bash
cd verification && python drive_from_p1r2.py
```

Python 3.12, standard library only; numpy is not required by the module.

Documentation:

- `docs/P2_A_SURFACE_MODULE.md` - equations, provenance, derivations,
  assumptions and open tensions
- `docs/P2_A_VERIFICATION.md` - numerical evidence for task item 10
- `docs/P2_A_INTERFACE_COMPARISON.md` - evaporation-interface options A and B
- `docs/P2_B_EVOLVING_POROSITY.md` - removing the frozen-porosity assumption
- `docs/P2_B_SALT_GROW_PARAMETER_COMPARISON.md` - P2-C: transferring the
  salt_grow crust parameters, and the benchmark against its transient ALE model
- `docs/P2_A_GO_NO_GO.md` - decision and gates

Read `docs/PRAT_PHYSICS_AUDIT.md` and `docs/INTERFACE_CONSERVATION.md` first;
every Prat equation used here is taken from the audit, not from memory.

**Status: one-way diagnostic only.** Growth is gated behind an explicitly
declared seed event and the evaporation flux is never fed back into a Richards
solution. P2-B shows that with porosity evolving as equation (2) requires, the
crust seals after growing by only about `beta * eps_0 / a_v` **at the P2-A
declared parameters**; P2-C shows that this sealing is a property of those
placeholder values, not of the formulation, and that it disappears once the
salt_grow parameters are used. See the GO/NO-GO before using any number here.
