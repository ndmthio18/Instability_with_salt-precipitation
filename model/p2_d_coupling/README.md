# P2-D first conservative Richards-Prat coupling

Bounded 1D proof of concept for coupling the Prat (2026) porous-salt crust to
the P1-R2 conservative Richards-NaCl model, plus the isolated C++ module the
production solver would use.

```
cpp/                        item 3: isolated C++ crust module
  pratcrust.hh                header-only, dependency-free C++17 port
  verify_vectors.cc           standalone checker against the reference vectors
  reference_vectors.txt       64 vectors from the verified Python model
  CMakeLists.txt, README.md
column/                     items 4-9: the 1D coupling
  soil.py                     1D Richards-NaCl column, a mirror of P1-R2's closures
  interface.py                seed, coupled flux conditions, compatible-flux solve
verification/
  config/acceptance_criteria_p2d.json   tolerances and the predictions
                                        registered before measurement
  reference_config_p2d.py     declared soil, crust and film configuration
  test_p2d_coupling.py        47 tests
  run_p2d.py                  regenerates the evidence tables
  generate_cpp_vectors.py     regenerates cpp/reference_vectors.txt
results/p2d_coupling/       generated evidence (CSV + JSON)
```

## Status

**The C++ module has never been compiled.** This environment has no C++
toolchain (`g++`, `clang++`, `cl`, `cmake`, `make`, `ninja` are all absent)
and the pinned DuMuX/Dune stack is not installed, so `model/p1_r2_...` cannot
be built or run either. Items 4 to 9 are therefore demonstrated against
`column/soil.py`, an isolated 1D mirror of P1-R2's closures, and item 3 ships
with the reference-vector gate that would transfer verification to the C++ on
a machine that can build it. See `cpp/README.md` and
`docs/P2_D_COUPLING.md` sections 2 and 10.

## Reproduce

```bash
cd verification && python -m unittest test_p2d_coupling
```

```bash
cd verification && python run_p2d.py
```

```bash
cd verification && python generate_cpp_vectors.py
```

On a machine with a compiler:

```bash
cd cpp && g++ -std=c++17 -O2 -o verify_vectors verify_vectors.cc && ./verify_vectors reference_vectors.txt
```

Python 3.12, standard library only. No scipy, no numpy, no Dune.

## What is coupled, and how

The soil surface evaporates through a **water-activity-controlled gas film**:
nothing prescribes a rate, and the flux falls as the surface concentrates.
When the surface reaches halite saturation a crust seed is formed **out of the
soil's own salt**, so the combined inventory is unchanged. From then on the
interface NaCl mass fraction `w_i` is the single unknown of

```
g(w_i) = F_s,soil(w_i) - F_s,crust(w_i) = 0
```

solved each step by a safeguarded Illinois iteration. Neither the interface
concentration nor the interface salt flux is prescribed. The soil is handed
the crust's **inventory-route** fluxes, so the pair conserves NaCl exactly;
the crust's independent **transport-route** reconstruction is reported
alongside as an audit, and differs only by the audited linearised-sink
residual.

The model stops rather than extrapolating on loss of hydraulic connectivity,
on an interface outside the declared halite-equilibrium envelope, on an
arrested top closure, on leaving the Pitzer activity range, on `Pe` not being
small against `Da`, and on a failed Newton solve. The growth velocity is never
capped.

## Relationship to the rest of the repository

* `model/p2_a_prat_surface/surface/graded_crust.py` is reused unchanged except
  for an optional base-concentration argument whose default reproduces the
  committed model exactly; the 44 P2-C tests enforce that.
* No frozen P1-R2 source, input or result is modified.
* `docs/PRAT_PHYSICS_AUDIT.md` section 3.3 now carries Prat equations (16) and
  (17) and the parameter provenance, which resolves the factor-six
  specific-area discrepancy between this repository and `salt_grow`.
