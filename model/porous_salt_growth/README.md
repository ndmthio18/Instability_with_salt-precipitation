# Standalone porous-salt growth

## Status

**COHERENT AFTER REFRAMING / EXTRA ASSUMPTION**

This module isolates the porous-salt layer from the soil model. It does not
compute soil saturation, soil concentration, soil diffusion, evaporation, or
CT coverage. The caller supplies `j(t)` and either the net salt-supply ratio
`s(t)` or downward salt flux `B(t)`.

The conserved object is

```text
P_net = G_top + R_positive - D_internal.
```

Top position, bottom position, thickness, and net solid mass remain separate.

## Finite-layer response

With `u=B/(Ceq*j)`, `s=1-u`, `H=h/ell`, and `X=sqrt(Da*)`, the module solves

```text
y'' - y = 0,
y'(0) = u,
y'(H) + X y(H) = 1.
```

It returns top growth, positive internal precipitation, internal dissolution,
and net solid growth, and checks

```text
G_top + R_positive - D_internal = s * Ceq*j/(1-Ceq).
```

This fixed-geometry response does not specify bottom motion. Prat's finite
thickness result is recovered at `u_c=1/[cosh(H)+X sinh(H)]`.

## Thick moving-front closure

The extra closure requires `H >> 1`, a wet connected layer, reversible
dissolution, and a thin bottom front. It gives

```text
epsilon_w = epsilon_top - (1-epsilon_top)/X,
v_bottom = u * v_Prat,
dh/dt = v_top - v_bottom = s * v_Prat.
```

The module rejects this closure when `epsilon_w <= 0`. This moving-front law
has been checked inside the reduced model; it is not a validated general law
for thin, dry, fractured, or detached crusts.

## API

```python
from porous_salt_growth import PorousSaltModel

salt = PorousSaltModel(
    equilibrium_mass_fraction=0.2634,
    crystal_density_kg_m3=2165.0,
    top_porosity=0.43,
    reaction_length_m=2.4729683670962094e-6,
    sqrt_da_star=156.0**0.5,
)

finite = salt.finite_layer_response(
    height_m=20e-6,
    evaporation_flux_kg_m2_s=8e-5,
    supply_ratio=0.70,
)

thick = salt.thick_layer_response(
    evaporation_flux_kg_m2_s=8e-5,
    supply_ratio=0.70,
)

states = salt.integrate_thick_layer(
    times_s=[0, 3600, 7200],
    evaporation_flux_kg_m2_s=[8e-5, 8e-5, 8e-5],
    supply_ratio=[0.70, 0.70, 0.70],
)
```

Use `from_transport_parameters(...)` to derive `ell` and `X` from one
consistent parameter set. Supplying reduced `ell` and `X` directly supports
reproduction of a reported `Da*`, such as Prat's 156; the two routes are not
silently mixed.

## Validation

```bash
cd model/porous_salt_growth
python -m unittest -v test_porous_salt_growth.py
```

The six tests cover Prat recovery, parameter reduction, signed internal
reaction, salt and kinematic closure, constant-forcing integration, and
rejection of a disconnected thick layer.

## Non-claims

- No soil, evaporation, coverage, or CT observation model is embedded.
- The finite-layer response does not provide a bottom velocity.
- The thick moving-front closure is not used outside its declared regime.
- Initial nucleation and the `H=O(1)` to thick moving-front transition remain
  external modelling choices.
