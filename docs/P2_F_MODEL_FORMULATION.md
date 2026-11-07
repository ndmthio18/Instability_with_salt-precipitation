# P2-F reaction-free analytical surface-growth formulation

## Evidence status

The new law below is **derived for P2-F** from Prat's moving-interface balances
and surface kinetic postulate. It is not stated explicitly as such in Prat
(2026). The paper's equations (28), (29) and (36) include the consequences of
distributed internal precipitation and cannot simply be reused after declaring
that process absent.

## Reaction-free derivation

Let the already activated crust have constant porosity `eps`, constant liquid
density `rho`, effective diffusivity `D`, height `h`, and no volumetric source.
Under quasi-steady, negligible-advection transport,

```
d2C/dz2 = 0,
C(0) = C_eq,
C(h) - C_eq = DeltaC_top.
```

The gradient is `DeltaC_top/h`. Prat's surface kinetic postulate, equation
(20), gives

```
v = beta k_r rho DeltaC_top / rho_cr.
```

Prat's moving total and salt balances, equations (6) and (7), combined with
`C approximately C_eq` and negligible vapor swept volume, give equation (8):

```
rho eps D dC/dz|h
  = C_eq j - rho_cr (1-eps) (1-C_eq) v.
```

Eliminating `DeltaC_top` yields

```
v = [C_eq j / (rho_cr (1-eps) (1-C_eq))]
    [Da_top / (1 + Da_top)],

Da_top = beta (1-eps) (1-C_eq) h k_r / (eps D).
```

Algebraically this is also the joint limit of equation (28) when internal area
`a_v -> 0`: `Da -> 0`, `Da* -> infinity`, while
`sqrt(Da Da*) = Da_top` and `tanh(sqrt(Da))/sqrt(Da) -> 1`. That identity alone
does not validate the transport approximation, because the paper's original
condition `Pe << Da` becomes impossible as `Da -> 0`.

## Correct reaction-free validity conditions

Retaining liquid advection gives

```
U dC/dz = eps D d2C/dz2,
Pe_U = |U| h / (eps D).
```

The exact top-gradient multiplier is
`B(Pe_U)=Pe_U exp(Pe_U)/(exp(Pe_U)-1)`. The corresponding implicit velocity is

```
v = G Da_top / [Da_top + B(Pe_U(v))],
G = C_eq j / [rho_cr (1-eps) (1-C_eq)],
rho U = j + [rho eps + rho_cr(1-eps) - rho_v] v.
```

The P2-F closed law requires `Pe_U << 1`. Quasi-steady transport in the growing
domain independently requires

```
Pe_move = v h / D << 1,
```

which is the ratio of diffusion time `h^2/D` to growth time `h/v`. P2-F uses
predeclared numerical limits `Pe_U <= 0.1` and `Pe_move <= 0.1`; these are model
acceptance choices, not thresholds stated by Prat.

## State and inventories

The production state is only activation status, height and revision. Porosity
is immutable. There is no crust grid, concentration profile, dynamic porosity,
specific internal area, or internal precipitation source.

For a saturated equilibrium pore liquid,

```
M_solid = rho_cr (1-eps) h,
M_s,d   = rho eps C_eq h,
M_water = rho eps (1-C_eq) h.
```

Therefore the upward soil-to-crust component demands are

```
F_s = [rho_cr(1-eps) + rho eps C_eq] v,
F_w = j - rho_v v + rho eps (1-C_eq) v.
```

The water convention treats `j` as the laboratory-frame vapor mass flux; set
`rho_v=0` for the prescribed-flux benchmark. These fluxes are inventory
identities. They do not prove that an unsaturated soil can hydraulically supply
the assumed saturated crust.

## Coupling and limitations

Before activation, DuMuX uses prescribed water evaporation and zero external
salt flux. Reaching the solubility mass fraction is an eligibility event, not a
nucleation-rate prediction. A declared seed is transferred from soil inventory
to crust inventory at an accepted step. After activation, `neumann()` evaluates
only a trial; height is committed once after the global soil step is accepted.

The model monitors interface composition, saturation, `Pe_U`, `Pe_move`, and
component budgets. It flags or stops outside the declared domain and never
clips `v`. It does not contain a crust hydraulic-pressure equation, film flow,
dissolution/recession law, pore filling, or soil precipitation. Consequently a
successful numerical run is only conditional evidence inside this envelope.
