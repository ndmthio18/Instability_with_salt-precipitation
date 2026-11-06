# P2-A isolated surface subgrid

Standalone external porous-salt surface module: heterogeneity gas-film
evaporation plus the Prat (2026) growth law. It does **not** modify P1-R2, any
frozen source, or any existing result, and it runs no DuMuX simulation.

```
surface/            the module
  thermo.py           source-traced thermophysical correlations
  water_activity.py   NaCl water activity (Pitzer), reused from heterogeneity
  film_evaporation.py physical gas-film evaporation flux
  prat_crust.py       Prat growth law, reduced crust transport, inventories
  coupling.py         one-way chain, base-flux reconstruction, options A/B
verification/
  config/acceptance_criteria.json   tolerances and declared parameters,
                                    registered before the runs
  reference_config.py               the single declared reference configuration
  test_p2a_surface.py               61 tests, maps onto task item 10
  run_verification.py               regenerates the evidence tables
  drive_from_p1r2.py                one-way chain driven by frozen P1-R2 data
results/            generated evidence (CSV + JSON)
```

Reproduce:

```bash
cd verification && python -m unittest test_p2a_surface
```

```bash
cd verification && python run_verification.py
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
- `docs/P2_A_GO_NO_GO.md` - decision and gates

Read `docs/PRAT_PHYSICS_AUDIT.md` and `docs/INTERFACE_CONSERVATION.md` first;
every Prat equation used here is taken from the audit, not from memory.

**Status: one-way diagnostic only.** Growth is gated behind an explicitly
declared seed event, the evaporation flux is never fed back into a Richards
solution, and the constant-porosity assumption of the reduced crust solution
does not hold at the declared parameters. See the GO/NO-GO before using any
number here.
