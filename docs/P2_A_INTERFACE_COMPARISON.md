# P2-A: evaporation-interface options A and B

This answers P2-A item 9. Evidence:
`model/p2_a_prat_surface/results/p2a_verification/interface_comparison.csv`
(200 rows) and `summary.json -> interface_comparison`.

## 1. The two options

**Option A - soil-top concentration.** The water activity driving the gas film
is evaluated at the soil-top (crust base) NaCl mass fraction. Under Prat's base
equilibrium that is `C_eq`. This is the simplified approximation: it ignores
the supersaturation that the crust transport solution builds between base and
top.

**Option B - Prat crust-top concentration.** The water activity is evaluated at
`C_top = C_eq + C*_top` from the reduced crust transport solution. Because
`C*_top` depends on the evaporation flux through `Pe`, option B is implicit:

```
j = film( a_w( C_eq + C*_top(j) ) )
```

It is solved by bracketed bisection. `g(j) = film(C_top(j)) - j` is strictly
decreasing (`C*_top` grows linearly with `j`, the film flux falls with
`C_top`), so the bracket `[0, j(C_eq)]` is guaranteed and no relaxation
constant is tuned. The residual converges below 1e-16 and the recovered flux
reproduces `film(C_top(j))` to better than 1e-9.

Before activation there is no crust, `C*_top` is undefined, and the two options
return bitwise identical fluxes.

## 2. Direction and size of the difference

Option B always gives the **smaller** flux: the crust top is supersaturated
relative to its base, so its water activity is lower and the vapour driving
force is weaker. Over the 194 evaluated sweep rows
(`h` from 1e-6 to 3.16e-2 m, `a_v` from 1e3 to 1e7 m^-1, `k_r` from 1e-9 to
1e-6 m/s):

| quantity | value |
|---|---|
| rows evaluated / model-domain exits | 194 / 6 |
| acceptable at the declared 1 % level | 143 (73.7 %) |
| maximum relative evaporation-flux difference | 13.13 % |
| maximum relative growth-velocity difference | 13.13 % |
| worst case | `h = 3.16e-3 m`, `a_v = 1e3 m^-1`, `k_r = 1e-9 m/s`, `Da = 0.0194` |

The six exits are rows where `C*_top` drives the crust-top molality past 8,
the validity limit of the inherited NaCl activity correlation. They are
recorded as model-domain exits, not extrapolated.

## 3. When option A is acceptable

### 3.1 The controlling quantity is `C*_top`, not `Pe` and not `Da`

Option A omits exactly `C*_top`. To first order its relative evaporation-flux
error is therefore

```
|dj|/j  ~=  S * C*_top ,        S = -d ln j / dC  evaluated at C_eq
```

At the declared film (`T = 293.15 K`, `RH = 0.4`, `p = 1e5 Pa`),
**`S = 4.4267`**, so the 1 % acceptance level corresponds to

```
C*_top  <=  0.01 / S  =  2.259e-3          (about 0.85 % of C_eq)
```

The sweep confirms this sharply:

| boundary measured in the sweep | value |
|---|---|
| largest `C*_top` that is acceptable | 2.2110e-3 |
| smallest `C*_top` that is not acceptable | 2.3137e-3 |
| predicted threshold `0.01/S` | 2.2590e-3 |

The predicted threshold falls between the two, and the first-order estimate
reproduces the measured difference with a median relative error of 0.07 % and a
worst case of 5.5 %.

`Pe` and `Da` do **not** separate the cases on their own. The sweep contains
acceptable rows at `Pe = 1.235` and unacceptable rows at `Pe = 0.0123`, and
both acceptable and unacceptable rows occur in every `Da` band:

| `Da` band | acceptable | not acceptable |
|---|---|---|
| `Da < 1` | 91 | 18 |
| `1 <= Da < 10` | 10 | 7 |
| `10 <= Da < 1e3` | 16 | 18 |
| `Da >= 1e3` | 26 | 8 |

The reason is that for `Da >> 1`, `C*_top ~ C_eq Pe / sqrt(Da)`: it is the
*combination* that matters, not either group alone.

### 3.2 Practical reading for this problem

Acceptance by crust height, at the declared parameter sweep:

| `h` [m] | acceptable rows |
|---|---|
| 1e-6 to 1e-4 | 20 / 20 at every height |
| 3.16e-4 | 11 / 20 |
| 1e-3 to 3.16e-2 | 8 / 17-20 |

So option A is safe for a **thin** crust and becomes parameter-dependent once
the crust exceeds roughly 1e-4 m. Because `C*_top` grows with `h` while
`Da` also grows, no single height threshold is universal; the defensible rule
is the `C*_top` criterion of section 3.1, which the module evaluates directly.

### 3.3 Recommendation

1. Use the `C*_top <= 0.01/S` test at run time, not a height or `Da` rule.
   `S` must be recomputed if `T`, `RH` or `p_atm` change, since it is a
   property of the film, not of the crust.
2. Where the test fails, use option B. Its cost is a scalar bisection per patch
   per step, which is negligible next to the Richards Newton solve.
3. Do not mix: a run must record which option it used, as P1-R2 records its
   evaporation boundary mode.
4. Note that a 13 % evaporation-flux difference propagates directly into
   `v_s`, which is linear in `j`. The worst-case growth-velocity difference is
   the same 13 %. For a study whose target observable is a growth rate, that is
   not a negligible modelling choice.

## 4. Caveats

* The acceptance boundary is quantified **at the declared reference film**.
  `S = 4.4267` is specific to `T = 293.15 K`, `RH = 0.4`, `p = 1e5 Pa`. A drier
  atmosphere raises the flux but also raises `Pe` and hence `C*_top`, so the
  comparison must be redone if the atmospheric control changes.
* `a_v` and `k_r` are declared, not measured. The sweep spans four and three
  decades respectively to make the conclusion structural rather than
  point-valued, but the location of a *physical* crust within that sweep is
  unknown.
* Option B is still built on the reduced quasi-steady, uniform-property,
  `Pe << Da` crust solution. It is a better interface concentration within that
  closure, not an independent check on it.
* Both options inherit the 293.15 K caveat of the activity correlation and the
  constant-porosity tension documented in `P2_A_SURFACE_MODULE.md` section 7.1.
