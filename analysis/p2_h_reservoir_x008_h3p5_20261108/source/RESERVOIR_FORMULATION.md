# P2-H reservoir boundary and accepted-step budgets

This revision changes only the soil bottom boundary and the corresponding
diagnostic ledger. The soil, gas-film, surface precipitation, donor transfer,
and ALE equations are unchanged.

## Geometry and Dirichlet state

The column coordinate is `y=0...0.20 m`, increasing upward. With
`Problem.WaterLevel=0.01 m`, the hydrostatic reference water table is at
`y=-0.01 m`. The bottom Dirichlet pressure is therefore

```text
p_b = p_n - rho(x_b) g (0.01 m + y_b),  y_b = 0.
```

This parameter does **not** place the water table 0.01 m above the column
bottom. The accepted initial field and bottom Dirichlet state both use
`Problem.SaltMoleFraction=0.05`; `reservoir_config.csv` records both values and
the accepted initial extrema.

## Bottom component flux reconstruction

Flux is positive into the soil. The bottom BOX control volumes are isolated by
the conservative internal SCVF cut halfway between the first two vertex
levels. At each accepted backward-Euler endpoint,

```text
B_liquid = (M_bottom,new - M_bottom,old)/dt - F_cut,down
B_salt   = (S_bottom,new - S_bottom,old)/dt - F_salt,cut,down.
```

`F_cut,down` is summed from the same DuMuX `transportFlux` operator used by the
residual. This is a local boundary-face reconstruction, not the whole-domain
residual rearranged to force closure. The water influx is
`B_water = B_liquid - B_salt`.

Only accepted steps update

```text
I_liquid += dt B_liquid
I_salt   += dt B_salt.
```

These cumulative values are members of the transactional `SurfaceState`.
Rejected soil, ALE, activation, or injected rollback trials therefore leave all
accepted inventories unchanged.

## Coupled balances

For area-integrated inventories, the checked balances are

```text
M_water,soil + M_water,ALE + M_vapour - (I_liquid - I_salt)
    = M_water,initial,

M_salt,soil + M_salt,ALE,pore + M_salt,solid - I_salt
    = M_salt,initial.
```

The activation-transfer and internal solid quadrature checks remain separate.
`open_boundary_budget.csv` contains the instantaneous and cumulative component
fluxes, bottom-control-volume terms, full coupled inventories, and residuals.

## Verification-only switch

`P2H.VerificationZeroEvaporation=true` sets the gas-film flux to zero only for
the registered hydrostatic reservoir test. It defaults to false and is not
used in any pilot or production case.
