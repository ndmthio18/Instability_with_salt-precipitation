# P2-B: removing the frozen-porosity assumption

> **Correction (P2-D, 2026-11-07).** Section 6.5 calls `a_v = 6(1-eps_0)/d`
> "the paper's own cubic construction". It is not: Prat equation (16) is
> `a_v = 5(1-eps_s)/d_c` and equation (17) is `a_v = 36(1-eps_s)/d_sub`. The
> corollary should read `dh/d_c = beta eps_0/(5(1-eps_0))`, or
> `dh/d_sub = beta eps_0/(36(1-eps_0))` per substrate grain. The qualitative
> statement -- that a constant-area self-replicating crust seals before adding
> one crystal diameter of its own microstructure -- is unchanged and in fact
> strengthened. See `PRAT_PHYSICS_AUDIT.md` section 3.3.
>
> **Read with P2-C (`P2_B_SALT_GROW_PARAMETER_COMPARISON.md`).** Every sealing
> result below is evaluated at the P2-A **declared placeholder** `k_r = 1e-7`
> and `a_v = 1e5`, for which `sqrt(Da*) = 0.053`. P2-C transfers the sibling
> `salt_grow` project's values (`k_r = 2.3e-3`, `a_v = 820800`,
> `eps_0 = 0.43`), which give `sqrt(Da*) = 2.56`, and the crust then does not
> seal at all: it settles to a travelling reaction layer and grows steadily.
> The two regimes are separated by the single criterion
> `sqrt(Da*) <= (1-eps_0)/(eps_0-eps_stop)`, derived in P2-C section 4. The
> law `dh = integral beta/a_v d eps` below is correct, but it only governs the
> thin-crust regime `delta_r > h` that the declared parameters happen to sit
> in; for `Da >> 1` the controlling result is instead the wake plateau
> `eps_0 - (1-eps_0)/sqrt(Da*)`. Nothing in this document is withdrawn; its
> scope is narrower than its wording suggests.

## 1. Why

`docs/P2_A_GO_NO_GO.md` holds coupling behind two blockers. This document
addresses the first: Prat equation (2) forces the crust porosity to change
wherever `C != C_eq`, while the reduced solution that produces the closed
growth law holds `eps_c` fixed (`PRAT_PHYSICS_AUDIT.md` section 3,
reduction 3). P2-A quantified the tension - at the declared parameters the
depth-integrated internal reaction exceeds the top skeleton rate by
`1/sqrt(Da*)` = 18.9 - but could not resolve it, because a single uniform
porosity has nowhere to put the information.

P2-B evolves porosity and keeps everything else. `prat_crust.py` and
`coupling.py` are untouched; `evolving_crust.py` must reproduce them in the
frozen-uniform limit, and a test enforces that.

The paper itself lists porosity reduction and changing crystal size among its
unresolved items (audit section 4), so what follows quantifies a gap the
author already identified; it is not a claim that the paper is wrong within
its stated scope.

## 2. The material coordinate is already `z`

In Prat's geometry the base is fixed at `z = 0` and new material appears only
at the moving top `z = h(t)`. Existing material neither moves nor deforms.
Therefore `z` is itself a material coordinate, and equation (2) applies
pointwise with no Lagrangian mapping:

```
d eps(z,t) / dt  =  - a_v rho k_r C*(z,t) / rho_cr ,       0 < z < h(t)
```

New material enters at `z = h(t)` with the declared deposition porosity
`eps_0`. A material point's porosity is thus the time integral of the local
supersaturation it has seen since it was laid down.

This also supplies, for free, the deposition-age bookkeeping that
`INTERFACE_CONSERVATION.md` section 4 requires before a top/internal solid
split may be reported. P2-A had to track total solid only; P2-B can report the
split, and a test confirms it closes to round-off.

## 3. Generalised reduced transport problem

With depth-dependent `eps`, the quasi-steady, `Pe << Da`, linearised
reaction-diffusion problem becomes variable-coefficient:

```
d/dz [ eps(z) D_s,eff(z) dC*/dz ]  =  (1 - C_eq) a_v k_r C* ,
C*(0) = 0 ,   C*(h) = C*_top .
```

It is linear and homogeneous in `C*`, so `C*(z) = C*_top phi(z)` with `phi`
solving the same equation under `phi(0) = 0`, `phi(h) = 1`. One tridiagonal
solve per step gives `phi`, and hence the top gradient coefficient

```
A  =  rho (eps D)_top  dphi/dz|_h .
```

When `eps` is uniform this reduces to the analytic
`phi = sinh(sqrt(Da) z/h)/sinh(sqrt(Da))` of P2-A.

## 4. Closing relations, kept exact

P2-A inherited equation (8), derived with `C ~ C_eq` at the top. P2-B combines
equations (6) and (7) without that step. Eliminating `rho U(h)` and writing
`C_top = C_eq + s`:

```
rho eps_top D_top dC/dz|_h
    =  C_top j  -  rho_cr (1-eps_top)(1-C_top) v_s  -  C_top rho_v v_s
```

and the outer-surface kinetic relation of equation (20), in the form P2-A
verified to be the one compatible with equation (28),

```
v_s = beta k_r rho C*_top / rho_cr  ==  K s .
```

Substituting gives one scalar equation for `s`,

```
g(s) = A s - (C_eq+s) j + rho_cr (1-eps_top)(1-C_eq-s) K s + (C_eq+s) rho_v K s = 0,
```

quadratic and concave in `s`, with `g(0) = -C_eq j < 0` for `j > 0`. Bracketed
bisection solves it; no Newton tuning and no relaxation constant.

### 4.1 Algebraic reduction to equation (28)

Drop the `s^2` terms (i.e. put `C_top -> C_eq` exactly where equation (8)
does), set `rho_v = 0`, and take `eps` uniform. Then
`A = rho eps D sqrt(Da) / (h tanh(sqrt(Da)))` and

```
v_s = K s = C_eq j K / ( A + rho_cr (1-eps)(1-C_eq) K ) .
```

Dividing through by `rho_cr (1-eps)(1-C_eq) K` and using

```
A / ( rho_cr (1-eps)(1-C_eq) K )  =  1 / ( sqrt(Da*) tanh(sqrt(Da)) )
```

- which follows from `sqrt(Da) = h sqrt((1-C_eq) a_v k_r/(eps D))` and the
definition of `Da*` - gives

```
v_s = G / ( 1 + 1/(sqrt(Da*) tanh(sqrt(Da))) )
    = G sqrt(Da*) tanh(sqrt(Da)) / ( 1 + sqrt(Da*) tanh(sqrt(Da)) )
```

with `G = C_eq j / (rho_cr (1-eps)(1-C_eq))`. That is equation (28) exactly.
The numerical counterpart of this proof is in section 6.1.

## 5. Declared closures, default off

Reducing porosity physically changes both the effective diffusivity and the
specific crystal area. Prat (2026) supplies no closure for either. The module
exposes

```
D_s,eff(eps) = D_ref (eps/eps_0) ** diffusivity_exponent
a_v(eps)     = a_v0  (eps/eps_0) ** specific_area_exponent
```

with **both exponents zero by default**, so the bare effect of equation (2) is
isolated. Non-zero exponents are declared sensitivities, never fits.

## 6. Numerical evidence

Tolerances and three predictions were registered in
`model/p2_a_prat_surface/verification/config/acceptance_criteria_p2b.json`
before the runs. 28 tests in `verification/test_p2b_evolving.py` pass.
Evidence tables are in `model/p2_a_prat_surface/results/p2b_evolving_porosity/`.

Reproduce:

```bash
cd model/p2_a_prat_surface/verification && python -m unittest test_p2b_evolving
```

```bash
cd model/p2_a_prat_surface/verification && python run_p2b.py
```

### 6.1 Reduction to equation (28)

Table: `eq28_reduction.csv`.

With the same `C ~ C_eq` top linearisation and uniform frozen porosity, the
variable-coefficient solver converges to Prat equation (28) at clean second
order:

| cells per seed height | relative difference | observed order |
|---|---|---|
| 8 | 2.614e-4 | - |
| 16 | 6.535e-5 | 2.0000 |
| 32 | 1.634e-5 | 1.9999 |
| 64 | 4.085e-6 | 2.0000 |
| 128 | 1.021e-6 | 2.0000 |
| 256 | 2.553e-7 | 1.9999 |
| 512 | 6.383e-8 | 2.0000 |

Without that linearisation the converged answer sits 3.6337e-3 away from
equation (28), against a predicted offset `C*_top/C_eq` = 3.5927e-3, i.e.
agreement to 1.1 %. So the residual is not a discretisation error: it is the
top-linearisation error that equation (28) itself carries, and the new solver
removes it.

This is the strongest single check in P2-B. It exercises the transport solve,
the top closure, the kinetic relation and the parameter plumbing at once, and
it ties the new model to the audited closed form.

### 6.2 Reference trajectory: the crust seals

Tables: `reference_trajectory.csv`, `porosity_snapshots.csv`.

At the declared reference inputs, with constant `a_v` and `D_s,eff`:

| quantity | value |
|---|---|
| sealing time | 5.60e4 s (15.6 h) |
| height at sealing | 1.0838e-4 m, from a 1.0e-4 m seed |
| growth before sealing | 8.384e-6 m, i.e. **+8.4 %** |
| growth velocity, start -> sealing | 1.061e-10 -> 3.061e-10 m/s (**x2.88**) |
| porosity at the base / minimum / top | 0.3956 / 0.0098 / 0.3474 |
| mean porosity at sealing | 0.2288 |
| cumulative internal / top precipitation | 3.69 |
| salt budget residual | 5.71e-15 |
| water budget residual | 8.01e-15 |
| solid top/internal split residual | 3.06e-16 |

The depth profile at sealing shows a sharp clogged layer and fresh material
above it:

| z / h_seed | porosity | C* |
|---|---|---|
| 0.016 | 0.3956 | 1.28e-5 |
| 0.328 | 0.3025 | 3.04e-4 |
| 0.641 | 0.1913 | 7.29e-4 |
| 0.891 | 0.0736 | 1.42e-3 |
| 0.953 | 0.0341 | 1.85e-3 |
| 1.016 | 0.0694 | 2.67e-3 |
| 1.073 | 0.3474 | 2.79e-3 |

Two results here are worth stating plainly.

**The seal forms at the original seed top, not at the base and not at the
current top.** That material is the oldest that was ever exposed to the full
top supersaturation, so it has integrated the largest reaction. The registered
prediction was `z_clog = h_seed`, and the measured location converges to it
exactly as the grid refines: 0.9687, 0.9844, 0.9922, 0.9961 for 16, 32, 64 and
128 cells, each being the centre of the last seed cell.

**Growth accelerates as the crust tightens.** `v_s` nearly triples before
sealing. The mechanism is that lower porosity reduces `eps D`, which impedes
the diffusive back-flux of salt from the supersaturated top towards the base,
so `C*_top` rises and `v_s = beta k_r rho C*_top / rho_cr` rises with it. This
is the opposite of the intuition that a clogging crust must grow more slowly,
and it means a frozen-porosity model under-predicts growth right up to the
moment the real crust stops.

### 6.3 Grid and time-step convergence

Table: `grid_convergence.csv`. Sealing time at `dt = 200 s`:

| cells per seed height | sealing time [s] | growth [m] | z_clog / h_seed |
|---|---|---|---|
| 16 | 57200 | 8.574e-6 | 0.9687 |
| 32 | 56200 | 8.431e-6 | 0.9844 |
| 64 | 55600 | 8.314e-6 | 0.9922 |
| 128 | 55400 | 8.277e-6 | 0.9961 |

The 64-to-128 change is 0.36 %, inside the registered 5 % gate. Time-step
sensitivity is far smaller: over `dt` = 400, 200 and 100 s the sealing time
moves by at most 0.5 %.

### 6.4 The growth-before-sealing law

Table: `sealing_law.csv`, 108 combinations of `a_v` (1e4 to 1e7 m^-1), `k_r`
(1e-8 to 1e-6 m/s), relative humidity (0.2, 0.4, 0.6) and deposition porosity
(0.3, 0.4, 0.6). All 108 seal.

Registered before measurement:

```
dh  =  beta * eps_0 / a_v
```

Derivation: a material point at the top seals after
`t = eps_0 rho_cr / (a_v rho k_r C*_top)`, while the top advances at
`v_s = beta k_r rho C*_top / rho_cr`. The product removes `k_r`, `rho`,
`rho_cr` and `C*_top` - and therefore the evaporation flux as well.

Nine rows were excluded because the grid did not resolve the reaction
boundary layer `h/sqrt(Da)` (fewer than 4 cells across it). That is a
resolution criterion on the *boundary layer*, not on the growth increment, and
a test confirms those cases move towards the law on refinement. Over the 72
remaining thin-crust rows:

| quantity | value |
|---|---|
| ratio of measured to predicted growth | 0.991 to 1.147, mean 1.052 |
| spread of growth across RH = 0.2, 0.4, 0.6 | **0.54 %** |
| spread of growth across `k_r` over two decades | 16.5 % |

The evaporation-independence is essentially exact, which is the sharpest part
of the claim: tripling the vapour driving force changes the growth achieved
before sealing by half a per cent. The residual `k_r` sensitivity and the
mean ratio slightly above one both come from the same source - once `dh`
stops being small compared with `h`, the clogging point no longer sits at the
full `C*_top`, so it seals later and the estimate becomes a lower bound.

### 6.5 Corollary in the paper's own cubic construction

The audit records `beta = 2` for the ideal cubic construction of the paper's
Appendix. For that same construction, a packing of cubes of side `d` at solid
fraction `1 - eps_0` has `a_v = 6(1-eps_0)/d`, so

```
dh / d  =  beta eps_0 / ( 6 (1 - eps_0) )  =  0.222   at eps_0 = 0.4 .
```

Across the whole sweep the measured `dh/d` lies between 0.14 and 1.03. **With
constant `a_v`, the crust seals before it has added one crystal diameter of
its own microstructure.** A self-replicating constant-porosity crust is
therefore not consistent with equation (2) at any evaporation rate. This is a
statement about the closure as audited, and it sits squarely on the item the
paper itself lists as unresolved.

Note that `a_v = 6(1-eps_0)/d` is a geometric corollary of the same idealised
construction, derived here; the audit does not record the paper's own `a_v`
value or pore geometry.

### 7. Declared closures: what would rescue growth

Table: `closure_sensitivity.csv`. Six combinations of the declared exponents,
all at the reference inputs. The analytic prediction generalises to any power
law, because the same cancellation survives:

```
a_v(eps) = a_v0 (eps/eps_0)**m
  =>  dh = (beta eps_0**m / a_v0) * integral_{eps_seal}^{eps_0} eps**-m d eps
      m = 0 :  beta eps_0 / a_v0
      m = 1 :  (beta eps_0 / a_v0) ln(eps_0/eps_seal)
      m = 2 :  (beta eps_0 / a_v0) (eps_0/eps_seal - 1)
```

still independent of `j`, `k_r` and `C*_top`. Measured:

| `m` for `a_v(eps)` | `m` for `D(eps)` | sealed | growth [m] | predicted [m] | ratio | growth vs constant |
|---|---|---|---|---|---|---|
| 0 | 0 | yes, 5.62e4 s | 8.431e-6 | 7.800e-6 | 1.08 | 1.0 |
| 0 | 1.5 | yes, 3.94e4 s | 8.884e-6 | 7.800e-6 | 1.14 | 1.05 |
| 1 | 0 | yes, 1.334e5 s | 3.767e-5 | 2.951e-5 | 1.28 | 4.5 |
| 1 | 1.5 | yes, 7.08e4 s | 4.636e-5 | 2.951e-5 | 1.57 | 5.5 |
| 2 | 0 | yes, 5.606e5 s | 6.919e-4 | 3.120e-4 | 2.22 | 82 |
| 2 | 1.5 | **no**, 4e6 s run | 1.348e-2 | 3.120e-4 | 43 | 1599 |

Every ratio is at or above one, as the lower-bound prediction requires, and
the ratio grows with `m` exactly because the crust then grows far beyond the
seed height and the thin-crust assumption lapses.

The decisive closure is `a_v(eps)`, not `D(eps)`. A tortuosity exponent alone
changes the growth before sealing by 5 %; it mostly accelerates everything, so
the crust seals sooner (3.94e4 s instead of 5.62e4 s) while reaching a similar
thickness. A shrinking specific area is what buys thickness: at `m = 2` the
crust reaches 0.69 mm, and with a tortuosity exponent as well it reaches
13.5 mm without sealing at all.

Physically this says that a crust can only keep growing if its crystals
coarsen as its porosity falls, so that the internal reaction shuts down faster
than the pore space does. That is precisely the "changing crystal size" item
the audit lists as unresolved in the paper. **P2-B does not supply that
closure.** Both exponents are declared sensitivities with no measurement
behind them, and none of the numbers in the `m > 0` rows should be quoted as
a prediction.

### 8. Frozen versus evolving porosity

Table: `frozen_versus_evolving.csv`. Both models from the same seed, the same
constant evaporation flux, over the window ending at the resolved model's
sealing time:

| quantity | value |
|---|---|
| window | 5.60e4 s |
| growth, evolving | 8.384e-6 m |
| growth, frozen (P2-A) | 6.072e-6 m |
| ratio | 1.381 |
| growth-velocity ratio at sealing | 2.753 |
| maximum relative height difference | 27.6 % |

So within its own window the frozen-porosity model **under**-predicts growth
by 38 %, and its instantaneous velocity is 2.75 times too low by the end. More
importantly, it has no sealing event at all: past 5.6e4 s it keeps growing
indefinitely, and the comparison stops being between two admissible models.
The P2-A error is therefore not a bounded bias that could be absorbed into a
parameter; the two models diverge qualitatively.

## 9. What P2-B does not establish

1. `a_v` and `k_r` are still declared inputs with no measured value. The
   sealing law is structural, but where a real crust sits on it is unknown.
   **P2-C update:** the salt_grow baseline supplies values with provenance,
   and they land on the non-sealing side of the criterion, with a margin of
   only 3.4 in either `a_v` or `k_r`. They are still not a measurement.
2. Sealing is declared at a porosity threshold, not derived from a percolation
   or permeability criterion. The true loss of hydraulic connectivity will
   occur at a porosity above zero and is material-specific.
3. The quasi-steady transport reduction, the `Pe << Da` neglect of advection,
   the linearised reaction term and the independent-patch (no lateral
   transport) approximation are all inherited from the audit unchanged.
   **P2-C update:** the first two are now benchmarked. Against salt_grow's
   transient ALE model with every input matched, the quasi-steady reduction
   costs less than 0.1 % in growth velocity, and the linearised reaction term
   is shown to be exactly the residual `- s L_top` between the two base-flux
   routes.
4. Nothing here is coupled to a Richards solution, and the second P2-A
   blocker - no C++ implementation - is untouched.
   **P2-C correction:** the `salt_flux_base_route` diagnostic in this module
   omits the pore-volume released by precipitation and uses the top-cell
   rather than the deposition porosity at the moving interface. That is most
   of the 7.5-10.4 % route difference reported in section 6. It is a
   diagnostic defect only - no growth velocity, porosity, height or budget
   here is affected - and it is corrected in `surface/graded_crust.py`.
5. A sealed crust is reported as a model-domain exit. What a real crust does
   after its pore network closes (fracture, detachment, efflorescence on a new
   site, vapour-mediated growth) is outside every closure used here.
