# P2-D isolated C++ crust module

`pratcrust.hh` is a header-only, dependency-free C++17 port of the minimum
crust functionality needed to couple the Prat (2026) growth law to the P1-R2
conservative Richards-NaCl solver.

## Status: written, not compiled

**This code has never been compiled or executed.** The environment it was
written in has no C++ toolchain at all - `g++`, `clang++`, `cl`, `cmake`,
`make` and `ninja` are all absent - and the pinned DuMuX/Dune stack recorded
in `model/p1_r2_conservative_richards/provenance/source_manifest.json` is not
installed. Nothing in this directory may be treated as verified until the
check below passes.

## The check that transfers verification

```bash
g++ -std=c++17 -O2 -o verify_vectors verify_vectors.cc
./verify_vectors reference_vectors.txt
```

or

```bash
cmake -S . -B build && cmake --build build && ctest --test-dir build
```

`reference_vectors.txt` holds 64 quantities produced by the Python model that
P2-A, P2-B and P2-C verified: the thermophysical correlations, the Pitzer
water activity, the gas film, Prat's equation (16)/(17) specific-area closure
and the `sqrt(eps)` tortuosity closure, the seeded graded grid, three
quasi-steady profile solves at different base concentrations, and a ten-step
integration with its inventory and transport-route fluxes. `verify_vectors.cc`
recomputes every one from the header and reports the worst relative
difference against a 1e-11 gate.

Regenerate the vectors with:

```bash
cd ../verification && python generate_cpp_vectors.py
```

The declared vector configuration is written into both files and is kept
adjacent in each, because a silent drift between them would make the
comparison meaningless.

## What is ported, and what is not

Ported, because the coupling needs it:

| | |
|---|---|
| thermo | IAPWS-IF97 saturation pressure, H2O-air gas diffusion, ideal binary gas density, mass/mole fraction conversions, the frozen brine density law |
| activity | Pitzer osmotic coefficient for NaCl and `a_w`, with the declared molality limit of 8 |
| film | stagnant gas-film evaporation, signed, no calibration and no normalisation of local weights |
| crust | Prat equation (17) `a_v = 36(1-eps)/d_sub`, `D_s,eff = sqrt(eps) D_s`, the graded top-refined grid, the quasi-steady tridiagonal transport solve with a **generalised base condition**, the exact top closure, equation (2) porosity evolution, the moving top, inventories, and the corrected base fluxes |
| coupling | the gas-film fixed point, the interface root solve, and `InterfaceFluxes` carrying exactly the two Neumann values P1-R2 needs |

Deliberately not ported: wake coarsening (`graded_crust.coarsen`), the P2-B
power-law closures, the legacy base-flux forms kept in Python only to measure
the correction, and the P2-C evidence runners. None is needed to couple, and
each would be dead code in the solver.

## Interface contract with P1-R2

P1-R2's two equations are total liquid mass (`conti0EqIdx`) and dissolved NaCl
mass (`conti1EqIdx`), and its `Problem::neumann` values are positive when mass
**leaves** the domain. `prat::InterfaceFluxes::totalMass` and `::saltMass` are
exactly those two numbers for the top boundary, in kg/(m^2 s):

```cpp
prat::CrustState crust = /* seeded once, when the surface reaches C_eq */;
const auto f = prat::solveInterface(
    crust, params, film,
    [&](double totalMassFlux, double interfaceMassFraction) {
        return /* the soil's own top-face salt flux reconstruction */;
    },
    dt);
values[conti0EqIdx] = f.totalMass;
values[conti1EqIdx] = f.saltMass;
```

The soil is given the **inventory** route, `d(M_solid + M_dissolved)/dt`, so
the coupled pair conserves NaCl exactly rather than to the linearisation
error of the reduced sink. `f.saltFluxTransportRoute` is the independent
transport reconstruction of the same flux and is there to be audited, not to
be used as the boundary value.

Neither the interface concentration nor the interface salt flux is
prescribed: `solveInterface` solves for the concentration at which the two
sides' salt fluxes agree, and the flux follows from it.

## Declared stops

`solveInterface` and `step` throw `prat::DomainExit` rather than extrapolating
when the crust loses hydraulic connectivity (the primary stop) or clogs, when
no compatible interface concentration exists inside the declared
halite-equilibrium envelope, when the top closure has no positive
supersaturation root, or when the crust-top molality leaves the activity
correlation's range. The growth velocity is never capped to match available
supply.
