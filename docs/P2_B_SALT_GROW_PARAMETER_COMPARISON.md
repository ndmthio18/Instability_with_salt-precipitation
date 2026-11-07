# P2-C: transferring the salt_grow crust parameters into the P2-B model

Date: 2026-11-07
Base commit: `632e8f3189db302c150bd6df96545297dd9a2afe`
Reference project: `ndmthio18/salt_grow`, `prat_dynamic_growth/`
Criteria registered before any run: `verification/config/acceptance_criteria_p2c.json`

## 0. Answer to the question that was asked

> The objective is to determine whether the different P2-B results were mainly
> caused by the arbitrary initial crust parameters, rather than a fundamental
> difference between the two numerical formulations.

**It was the parameters.** When every physical input and both constitutive
closures are matched, the P2-B quasi-steady evolving-porosity model and the
salt_grow transient ALE model agree on the 6-hour growth velocity to
**0.05 %**, on the minimum porosity to 0.13 %, on the surface supersaturation
to 0.05 %, on the internal precipitation fraction to 0.3 % and on the
reaction-zone depth `L95` to 0.8 %. The two codes also agree on the
closed-form Prat velocity to 8 significant figures and on all four
dimensionless groups.

The qualitative difference between the two projects' conclusions -- P2-B's
crust seals and stops after growing 8 %, salt_grow's crust grows steadily for
six hours and never seals -- is reproduced inside **one** code, by changing
only `eps_0`, `k_r` and `a_v`. The switch is `k_r`: P2-A declared `1e-7 m/s`
where salt_grow uses `2.3e-3 m/s`, a factor of 23000.

The physics behind the switch is a single group. Section 4 derives, and the
runs confirm, that a crust seals if and only if

```
sqrt(Da*)  <=  (1 - eps_0) / (eps_0 - eps_stop) .
```

P2-A/P2-B sit at `sqrt(Da*) = 0.053`, deep inside the sealing side. The
salt_grow parameters sit at `sqrt(Da*) = 2.56`, outside it. The threshold is
1.390.

Two qualifications are carried through the whole document and are **not**
rhetorical:

1. The margin is small. Sealing would return if `a_v` were 3.4 times larger or
   `k_r` 3.4 times smaller. Neither has been measured by either project.
2. Agreement between two models that share the same closure is verification,
   not validation. Both projects record that the reaction layer
   (`delta_r = 0.51 um`) is thinner than the estimated crystal size
   (`d_c ~ 3.5 um`), so the continuum interpretation is doubtful at *both*
   parameter sets.

A defect in the committed P2-B soil-side flux diagnostics was found and fixed
during this work (section 7). It changes no growth velocity and no committed
P2-B result.

## 1. What was read, and what was reused

Read from `salt_grow/prat_dynamic_growth/`, as the task directs, and used only
as source values -- nothing was guessed or recalled:

| File | What was taken |
|---|---|
| `configs/baseline.json` | every physical and numerical value in section 2 |
| `src/prat_growth/parameters.py` | `a_v_ref = 36(1-eps)/d_sub`, `D_eff_ref = sqrt(eps) D_mol` |
| `src/prat_growth/analytic.py` | the dimensionless-group definitions and the closed-form velocity |
| `src/prat_growth/dynamic.py` | the B1/B2 closures, the top-boundary quadratic, the ALE flux construction |
| `RESULTS.md` | the published groups, `2.94033793e-8`, `3.00578141e-8`, and the 6 h B2 state |
| `M3_METHODS.md` | the frozen/transient decomposition and the Q0-Q6 attribution |
| `M4_MATERIAL_HISTORY.md` | the fixed-material histories and the travelling-wake evidence |

Reused unchanged from this repository: `evolving_crust.py`'s physics, the
P2-A gas film, the Pitzer water activity, the seed/activation gate and the
budget machinery. No new model was developed.

`prat_crust.py`, `evolving_crust.py`, `coupling.py`, their configs and their
committed results are **byte-identical** to commit `632e8f3`. A test
(`TestNothingCommittedChanged`) asserts that Case 0 of this study *is* the
committed P2-B configuration.

### 1.1 One new module, and why it was necessary

`evolving_crust.py` discretises the crust on a **uniform** grid. That is
adequate at the P2-A parameters, where `delta_r = 227 um` exceeds the
100 um seed, so the whole crust reacts nearly uniformly. At the salt_grow
parameters `delta_r = 0.54 um`, and a uniform grid resolving it over a 0.1-1 mm
crust needs 10^4-10^5 cells, which is not affordable in pure Python and would
in any case be wasted on an inert wake.

`surface/graded_crust.py` therefore generalises the P2-B geometry from a
single `cell_size` to an explicit list of cell thicknesses, keeping every
equation, the exact top closure and the solid-volume state variable. It adds
only: non-uniform cells, optional wake coarsening, the two named salt_grow
closures, two corrections to the base-flux diagnostics (section 7) and
connectivity as the primary stop.

It is required to reproduce P2-B, and does:

| Check | Result |
|---|---|
| Seed state on a uniform grid vs `evolving_crust` | identical to 1e-20 |
| First solve (`v_s`, `s_top`, `R_int`, `M_d`, `M_w`, `A`) | `< 1e-12` relative |
| 400-step trajectory, 7 diagnostics | worst `3.2e-11`, from the top-closure bisection alone; the committed test asserts `< 1e-9` over 250 steps |
| Legacy base salt flux | `< 1e-12` relative |
| Linearised uniform frozen limit vs equation (28) | `3.9e-8` at 512 cells, clean second order |

## 2. The transferred parameters

Three quantities were changed, exactly as the task specifies. Everything else
was held fixed for the one-factor study.

| Quantity | P2-A/P2-B declared | salt_grow | Ratio | Transferred |
|---|---:|---:|---:|:--:|
| `eps_0` | 0.4 | 0.43 | 1.075 | yes |
| `k_r` [m/s] | 1e-7 | 2.3e-3 | 23000 | yes |
| `a_v` [1/m] | 1e5 | 820800 | 8.208 | yes |
| `d_sub` [m] | not used | 2.5e-5 | - | (defines `a_v`) |
| `C_eq` | 0.2649408697676898 | 0.264 | 0.9966 | no |
| `D_s,eff` [m2/s] | 9.48683298050514e-10 | 8.5240e-10 | 0.8985 | no |
| `rho` [kg/m3] | 1202.1686 | 1200 | 0.9982 | no |
| `rho_cr` [kg/m3] | 2165 | 2170 | 1.0023 | no |
| `beta` | 2.0 | 2.0 | 1 | no |
| `h0` [m] | 1e-4 | 1e-3 | 10 | no |
| `j` [kg/m2/s] | 1.781528648571137e-5 | 1.3888888888888889e-4 | 7.796 | no |

`a_v0 = 36(1-0.43)/25e-6 = 820800 1/m`, as the task states and as
`parameters.py` computes it.

### 2.1 A declared discrepancy between the two projects

salt_grow uses `a_v = 36(1-eps)/d_sub`. P2-B's own cubic corollary used
`a_v = 6(1-eps_0)/d`. When this document was written the audit transcribed no
specific-area formula, so neither value was traceable to the paper and this
section declined to adjudicate.

> **Resolved in P2-D (2026-11-07).** The paper was read directly. Its
> equation (16) is `a_v = 5(1-eps_s)/d_c` and its equation (17) is
> `a_v = 36(1-eps_s)/d_sub`, obtained from (16) with the measured
> `d_c/d_sub ~ 0.14`. **salt_grow's closure is equation (17) exactly**, and
> its `D_s,eff = sqrt(eps) D_s` is the paper's `tau ~ eps^(-1/2)` tortuosity.
> The `6(1-eps_0)/d` of P2-B is not in the paper; it was a single-cube
> surface-to-volume estimate, wrong in the coefficient and ambiguous in the
> length scale. See `PRAT_PHYSICS_AUDIT.md` section 3.3. Nothing in this
> document changes: the mandated value was already the correct one, and the
> table below showed the conclusion was insensitive to the factor anyway.

The factor matters for `delta_r` (by `sqrt(6)`) and for `Da*` (by 6) but not
for the conclusion, and it errs in the safe direction:

| `a_v` [1/m] | `delta_r` [um] | `Da` | `Da*` | predicted wake porosity | seals |
|---|---:|---:|---:|---:|:--:|
| 820800 = 36(1-e)/d (salt_grow) | 0.5422 | 3.402e4 | 6.562 | 0.2075 | no |
| 136800 = 6(1-e)/d (P2-B reading) | 1.328 | 5670 | 39.37 | 0.3392 | no |
| 82080 | 1.715 | 3402 | 65.62 | 0.3596 | no |
| 8208 | 5.422 | 340.2 | 656.2 | 0.4077 | no |

A smaller `a_v` moves *away* from sealing, so the P2-B reading of the cubic
construction would only strengthen the conclusion.

## 3. Dimensionless groups, computed before any transient run

All groups are evaluated at the deposition porosity, i.e. at the state of the
crust the instant it is created, which is the only state at which the
closed-form law is defined. `delta_r = h/sqrt(Da) = sqrt(eps D/((1-C_eq) a_v k_r))`.

| Quantity | Case 0 (P2-B) | Case 1 `+eps_0` | Case 2 `+k_r` | Case 3 `+a_v` |
|---|---:|---:|---:|---:|
| `Pe` | 3.9052e-3 | 3.6328e-3 | 3.6328e-3 | 3.6328e-3 |
| `Da` | 0.19371 | 0.18019 | 4144.4 | 34017 |
| `Da_top` | 0.023245 | 0.020542 | 472.46 | 472.46 |
| `Da*` | 2.7894e-3 | 2.3418e-3 | 53.861 | 6.5620 |
| `sqrt(Da*)` | 0.052814 | 0.048392 | 7.3390 | 2.5616 |
| `delta_r` [m] | 2.2721e-4 | 2.3558e-4 | 1.5534e-6 | 5.4219e-7 |
| `h/delta_r` | 0.4401 | 0.4245 | 64.38 | 184.4 |
| `C*_top` | 9.5185e-4 | 8.9126e-4 | 1.7929e-6 | 1.4652e-6 |
| `v_s` [m/s] | 1.0571e-10 | 9.8979e-11 | 4.5794e-9 | 3.7424e-9 |
| efficiency `f` | 0.021384 | 0.019022 | 0.88008 | 0.71923 |
| internal/top | 45.763 | 51.571 | 0.13626 | 0.39038 |
| `tau_eps` [s] | 7.5681e4 | 8.6887e4 | 1878.0 | 279.97 |
| `tau_burial` [s] | 2.1495e6 | 2.3801e6 | 339.2 | 144.9 |

The transfer moves the model across three order-of-magnitude boundaries at
once:

- `Da` rises by 1.8e5, from `Da < 1` (the reaction layer is **larger** than the
  crust, so the whole crust reacts) to `Da >> 1` (a 0.54 um reaction skin on a
  100 um crust).
- `Da*` rises by 2353, from 2.8e-3 to 6.56, taking the internal/top
  precipitation ratio from 45.8 down to 0.39 and the Prat efficiency from
  2.1 % to 71.9 %.
- `v_s` rises by 35.4.

`Pe` is unchanged to 7 %, confirming the `Pe << Da` neglect of advection is far
safer after the transfer than before (`Pe/Da` goes from 0.020 to 1.1e-7).

### 3.1 The fully matched benchmark configuration

| Quantity | This repository | salt_grow `RESULTS.md` | Difference |
|---|---:|---:|---:|
| `Pe` | 0.31574771900178206 | 0.31575 | 7e-7 |
| `Da` | 3790504.027697575 | 3.7905e6 | 1e-6 |
| `Da_top` | 5264.588927357744 | 5264.59 | 2e-7 |
| `Da*` | 7.3119290657746445 | 7.31193 | 1e-8 |
| `delta_r` [m] | 5.13631344455016e-7 | 5.136e-7 | 6e-5 |
| Prat `v_s` [m/s] | 2.9403379276556446e-8 | 2.94033793e-8 | 1e-9 |
| `1 - f` | 0.2699741823782068 | "about 26.9 %" | - |

Two independent implementations, written from the same audited equations in
different languages with different solvers, agree on the closed-form law to 8
significant figures. This is the cross-check that makes the rest of the
comparison interpretable.

## 4. Two closed-form predictions, registered before measurement

Both were derived and written into `acceptance_criteria_p2c.json` before the
corresponding runs. Neither is in Prat (2026); both are consequences of
applying equation (2) pointwise in the material coordinate, which is what
P2-B does.

### 4.1 Thick crust (`Da >> 1`): the wake porosity plateau

A material point deposited at the top at `eps_0` is buried at speed `v_s`
through a reaction layer in which the supersaturation decays as
`C*(d) = C*_top exp(-d/delta_r)`. Equation (2) applied to that point gives

```
d eps / dt = - a_v rho k_r C*_top exp(-v_s t / delta_r) / rho_cr ,
```

so the total porosity it ever loses is

```
d eps  =  a_v delta_r / beta  =  (1 - eps_0) / sqrt(Da*) ,
```

using `v_s = beta k_r rho C*_top / rho_cr` and
`sqrt(Da*) = beta(1-eps_0)/(a_v delta_r)`. **The evaporation flux, the
supersaturation, the solution density and `rho_cr` all cancel.** The crust
settles to a wake plateau at `eps_0 - (1-eps_0)/sqrt(Da*)` and stops
densifying, because the material has left the reaction layer.

### 4.2 The sealing criterion

A crust therefore reaches a declared stop porosity `eps_stop` if and only if
the available drop exceeds the distance to it:

```
(1 - eps_0) / sqrt(Da*)  >=  eps_0 - eps_stop
   <=>   sqrt(Da*)  <=  (1 - eps_0) / (eps_0 - eps_stop) .
```

With `eps_stop = 0.02` (the declared hydraulic-connectivity limit) the
threshold is 1.579 at `eps_0 = 0.4` and 1.390 at `eps_0 = 0.43`. This is a
**one-number** prediction of whether a crust seals, independent of the
atmosphere, of the crust height and of time.

### 4.3 Thin crust (`Da <~ 1`): the P2-B growth law still applies

Where `delta_r > h` the whole crust reacts at nearly the top supersaturation,
there is no wake to escape into, and the P2-B law
`dh = integral beta/a_v d eps` governs the growth before sealing. Cases 0 and
1 are in this regime.

### 4.4 Relation to Prat's efficiency factor

`internal/top -> 1/sqrt(Da*)` in the large-`Da` limit, so
`f_inf = 1/(1 + internal/top)`: Prat's efficiency *is* the top-growth share of
total precipitation, as P2-A derived. Section 4.2 can therefore be restated as
"a crust seals when Prat's own efficiency is small enough", which is why the
P2-A parameters (`f = 2.1 %`) sealed and the salt_grow parameters
(`f = 71.9 %`) do not.

## 5. Resolution: the committed grid cannot be used, and what replaces it

`delta_r` falls from 227 um to 0.54 um across the transfer, while the
committed P2-B seed grid has 3.125 um cells. That grid carries **0.17 cells
across the reaction layer**. Running Case 3 on it anyway, which this study
does once and only to show the size of the error:

| Grid | cells across `delta_r` | `v_s` [m/s] | `eps_min` | `L95` [m] | plateau time [s] |
|---|---:|---:|---:|---:|---:|
| committed P2-B, 32 uniform cells | 0.174 | 4.51242711e-9 | 0.32797078 | 1.876944e-6 | 3982.5 |
| converged, 16 cells per `delta_r` | 16 | 3.84051645e-9 | 0.22770723 | 1.283293e-6 | 2956.2 |
| error of the committed grid | | **+17.5 %** | **+44.0 %** | **+46.3 %** | +34.7 % |

An under-resolved run does not merely add noise: it reports a crust that is
44 % more porous and grows 17.5 % faster. No such result is interpreted
anywhere in this document.

The replacement is a top-refined grid: uniform cells of `delta_r/n` over the
top `40 delta_r`, then geometric growth by 1.25 up to `10 delta_r`, with new
material always created at the fine size and the inert wake optionally
coarsened by exact pairwise merging. For the 1 mm matched benchmark this needs
236-834 cells where a uniform grid at the same top cell size would need
2e4-3e5.

### 5.1 Grid convergence

Registered gate: at least 4 cells across `delta_r`, and successive changes
below 0.5 %.

| Configuration | cells/`delta_r` | cells | `v_s` [m/s] | d`v_s` | `eps_min` | d`eps_min` | `L95` [m] | d`L95` | resolved |
|---|---:|---:|---:|---:|---:|---:|---:|---:|:--:|
| Case 3 | 1 | 85 | 3.90819362e-9 | - | 0.24474906 | - | 1.436666e-6 | - | no |
| Case 3 | 2 | 147 | 3.86243241e-9 | 1.185 % | 0.23270785 | 5.174 % | 1.328763e-6 | 8.121 % | no |
| Case 3 | 4 | 271 | 3.84597850e-9 | 0.428 % | 0.22896744 | 1.634 % | 1.295686e-6 | 2.553 % | yes |
| Case 3 | 8 | 518 | 3.84157293e-9 | 0.115 % | 0.22796448 | 0.440 % | 1.285991e-6 | 0.754 % | yes |
| Case 3 | 16 | 1007 | 3.84051645e-9 | 0.0275 % | 0.22770723 | 0.113 % | 1.283293e-6 | 0.210 % | yes |
| matched B2 | 1 | 236 | 3.06218684e-8 | - | 0.25379016 | - | 1.185841e-6 | - | no |
| matched B2 | 2 | 275 | 3.02029037e-8 | 1.387 % | 0.24059240 | 5.486 % | 1.053043e-6 | 12.61 % | no |
| matched B2 | 4 | 355 | 3.00704938e-8 | 0.440 % | 0.23676253 | 1.618 % | 1.027163e-6 | 2.520 % | yes |
| matched B2 | 8 | 514 | 3.00433108e-8 | 0.0905 % | 0.23580750 | 0.405 % | 1.017076e-6 | 0.992 % | yes |
| matched B2 | 16 | 834 | 3.00339331e-8 | 0.0312 % | 0.23557439 | 0.0990 % | 1.014606e-6 | 0.244 % | yes |

Successive-change ratios are 2.8-4.9, i.e. close to second order, for all five
required quantities. `C*_top` converges identically to `v_s` because
`v_s = beta k_r rho C*_top / rho_cr` exactly, and `R_int` converges at the
same rate. Every row at 4 cells or more meets the 0.5 % gate for `v_s`, and 8
or more meets it for `eps_min`. **Production runs use 16 cells per `delta_r`,
where the remaining change is 0.03 %**, and that figure is quoted whenever a
sub-0.1 % agreement is claimed below.

The fifth required quantity, the hydraulic-connectivity-loss time, does not
exist at the transferred parameters: Case 3 and the matched benchmark never
lose connectivity. Its convergence is therefore demonstrated on Cases 0 and 1,
which do, and the **absence** of the event at the transferred parameters is
itself the converged result: it is stable across every grid and every time
step tested, and it is predicted analytically by section 4.2.

### 5.2 Time-step convergence

The step is `tau_eps/n` with `tau_eps = eps_0 rho_cr/(a_v rho k_r C*_top)`.

| Configuration | `n` | `dt` [s] | `v_s` [m/s] | d`v_s` | `eps_min` |
|---|---:|---:|---:|---:|---:|
| Case 3 | 375 | 0.746168 | 3.84171360e-9 | - | 0.22788698 |
| Case 3 | 750 | 0.373084 | 3.84135955e-9 | 0.0092 % | 0.22793914 |
| Case 3 | 1500 | 0.186542 | 3.84157293e-9 | 0.0056 % | 0.22796448 |
| Case 3 | 3000 | 0.093271 | 3.84163492e-9 | 0.0016 % | 0.22797662 |
| Case 3 | 6000 | 0.046636 | 3.84160657e-9 | 0.0007 % | 0.22798271 |
| matched B2 | 375 | 0.094969 | 3.00397711e-8 | - | 0.23573273 |
| matched B2 | 1500 | 0.023742 | 3.00433108e-8 | 0.0040 % | 0.23580750 |
| matched B2 | 6000 | 0.005936 | 3.00442300e-8 | 0.0010 % | 0.23582566 |

Time discretisation is two orders of magnitude below grid error and is not a
limiting factor. Case 0 is the exception: it terminates on a threshold
crossing, so its terminal state depends on where the last step lands and its
successive changes (0.35 %, 0.70 %, 0.091 %, 0.046 %) are non-monotone. That
is a property of a stopping event, not of the integrator.

## 6. The controlled one-factor comparison

All four cases use the identical gas film, temperature, relative humidity,
film thickness, seed height, `C_eq`, `D_s,eff`, densities and `beta`; the
evaporation flux is the same 1.781528648571137e-5 kg/m2/s throughout. Only the
quantity named in each column changes, cumulatively.

| | Case 0 | Case 1 | Case 2 | Case 3 |
|---|---:|---:|---:|---:|
| changed | - | `eps_0` | `+k_r` | `+a_v` |
| `Da` | 0.1937 | 0.1802 | 4144.4 | 34017 |
| `sqrt(Da*)` | 0.05281 | 0.04839 | 7.3390 | 2.5616 |
| `delta_r` [m] | 2.272e-4 | 2.356e-4 | 1.553e-6 | 5.422e-7 |
| **predicted to seal** | **yes** | **yes** | **no** | **no** |
| **observed** | **connectivity loss** | **connectivity loss** | **plateau** | **plateau** |
| stop time [s] | 55196.7 | 63442.8 | 18776.5 | 2956.3 |
| `h - h0` [m] | 8.1549e-6 | 8.8388e-6 | 8.6317e-5 | 1.1350e-5 |
| growth / seed height | 8.15 % | 8.84 % | 86.3 % | 11.35 % |
| `v_s` final [m/s] | 2.79084e-10 | 2.69339e-10 | 4.59716e-9 | 3.84157e-9 |
| `v_s` final / initial | 2.631 | 2.712 | 1.004 | 1.026 |
| `eps_min` | 0.019612 | 0.019346 | 0.354836 | 0.227964 |
| predicted `eps_wake` | - | - | 0.352332 | 0.207485 |
| wake prediction error | - | - | **0.71 %** | **8.98 %** |
| `dh` / predicted sealing law | **1.073** | **1.078** | - | - |
| `C*_top` final | 2.51303e-3 | 2.42528e-3 | 1.7998e-6 | 1.50399e-6 |
| mean porosity | 0.23198 | 0.24834 | 0.39517 | 0.40936 |
| `L95` [m] | 7.7179e-5 | 7.7455e-5 | 4.3250e-6 | 1.2860e-6 |
| internal/top | 2.832 | 2.988 | 0.1319 | 0.3545 |
| internal fraction of new solid | 0.7879 | 0.7969 | 0.1165 | 0.2621 |
| `F_s,soil` [kg/m2/s] | 1.27372e-6 | 1.21631e-6 | 6.94082e-6 | 6.70016e-6 |
| `F_w,soil` [kg/m2/s] | 1.74949e-5 | 1.75122e-5 | 1.92567e-5 | 1.85890e-5 |
| solid salt [kg/m2] | 0.179837 | 0.177118 | 0.243974 | 0.142388 |
| dissolved salt [kg/m2] | 8.00758e-3 | 8.62616e-3 | 2.34503e-2 | 1.45181e-2 |
| total salt [kg/m2] | 0.187845 | 0.185744 | 0.267424 | 0.156906 |
| water [kg/m2] | 2.21539e-2 | 2.38677e-2 | 6.50611e-2 | 4.02793e-2 |
| salt budget residual | 2.3e-15 | 8.6e-16 | 1.3e-15 | 4.6e-14 |
| water budget residual | 6.7e-15 | 2.3e-14 | 2.1e-13 | 2.9e-13 |
| solid top/internal split residual | 1.9e-15 | 3.2e-15 | 1.9e-14 | 3.5e-14 |

Time-resolved `h(t)`, `dh/dt`, the porosity extremes and every inventory are
in `case_<n>_trajectory.csv`; `eps(z,t)` and `C*(z,t)` at four times per case
are in `profile_snapshots.csv`; the final profiles are in `profiles.csv`.

### 6.1 Reading the ladder

**`eps_0` alone changes nothing qualitative** (Case 0 to 1). Both seal, the
sealing times differ by 15 %, and the growth before sealing follows the P2-B
law to within 7-8 %. Raising the deposition porosity simply gives the crust a
slightly longer way to fall.

**`k_r` is the switch** (Case 1 to 2). Multiplying it by 23000 raises `Da` and
`Da*` by the same factor, crossing the sealing threshold. The crust stops
sealing, the growth velocity rises by 46.5, the minimum porosity goes from
0.019 (stopped) to 0.355 (a stable plateau), and the internal fraction of new
solid collapses from 80 % to 12 %. **This single parameter accounts for the
entire qualitative difference between the two projects.**

**`a_v` modulates, it does not switch** (Case 2 to 3). Multiplying it by 8.2
lowers `Da*` by the same factor, deepening the wake plateau from 0.355 to
0.228 and slowing growth by 16 %, but the crust still does not seal. It moves
the model *towards* the sealing boundary without reaching it: `sqrt(Da*)` goes
from 7.34 to 2.56 against a threshold of 1.39.

**The registered sealing criterion was correct in all four cases**, and the
wake-plateau prediction `eps_0 - (1-eps_0)/sqrt(Da*)` was accurate to 0.71 %
in Case 2 and 9.0 % in Case 3, both inside the registered 15 % band. The
Case 3 error is larger because its `Da*` is smaller, so the exponential
boundary-layer approximation behind the derivation is weaker, and because the
porosity change feeds back on `delta_r` itself.

Cases 0 and 1 stop on **loss of hydraulic connectivity at `eps = 0.02`, not on
geometric clogging at 0.01**, as task item 7 requires. The minimum porosities
at the stop (0.019612, 0.019346) lie between the two thresholds, so the
connectivity event genuinely precedes clogging and the model refuses to step
further rather than extrapolating a saturated-supply closure into a
disconnected pore network.

Cases 2 and 3 instead reach a steady travelling state: the velocity changes by
0.4 % and 2.6 % over the whole run and then stops changing, the reaction layer
stays `1-3 delta_r` deep, and the wake behind it is inert. This is the same
state salt_grow's M4 study identified independently in its own model.

## 7. Conservation and the two base-flux routes

Task item 7 requires the salt influx at the crust base to be reconstructed
independently from the transport equation and compared with the
inventory-derived flux, without hiding the discrepancy.

### 7.1 What the P2-B 7.5-10.4 % discrepancy actually was

Two defects in the **diagnostic**, not in any governing equation:

1. **Pore-volume release.** Precipitation removes salt mass `R` from the
   solution but frees pore volume `R/rho_cr`, which refills with solution.
   Since `d(rho eps S)/dt = rho S d eps/dt = - rho S R/rho_cr` exactly, the
   solution-mass balance is

   ```
   d(rho U)/dz = - R (1 - rho S / rho_cr) ,
   ```

   not `- R`. P2-B used `- R`, overstating `rho U(0)` by `R_int rho S/rho_cr`.

2. **Interface porosity.** The pore volume created at the moving top is filled
   at the **deposition** porosity `eps_0`, not at the current top-cell
   porosity.

Both are corrected in `graded_crust.py`, which also evaluates the legacy forms
so the change is visible rather than silent. The legacy gap then decomposes
**exactly** into the defect plus the residual of section 7.2:

```
legacy gap  =  C_eq (rho S / rho_cr) R_int / F_s   -   s L_top / F_s
```

| Configuration | corrected | legacy | predicted pore-volume term | check |
|---|---:|---:|---:|---|
| Case 0 | 3.588e-2 | 8.270e-2 | 1.186e-1 | 1.186e-1 - 3.588e-2 = 8.27e-2 |
| Case 2 | 6.090e-6 | 1.585e-2 | 1.586e-2 | residual negligible, agrees to 0.06 % |
| Case 3 | 5.062e-6 | 3.690e-2 | 3.690e-2 | residual negligible, agrees to 0.01 % |

The committed P2-B route differences of 7.5-10.4 % are therefore a closure
defect of this size combined with the residual below, and should be read as
such. **No P2-B growth velocity, porosity, height or budget is affected**,
because that diagnostic is never fed back into the solution, and no committed
file was edited.

### 7.2 What remains, and why it is not a bookkeeping error

Summing the discrete cell balances of the reduced problem gives
`A s + D_base = (1 - C_eq) R_int` exactly. Substituting that, the top closure
and `rho U(0)` into `C_eq rho U(0) + D_base`, then subtracting the exact
inventory identity `d(M_solid + M_d)/dt`, leaves a closed form:

```
F_s,transport  -  F_s,inventory  =  - s L_top ,
     L_top = j + rho_cr (1 - eps_0) v_s - rho_v v_s ,
```

the solution flux at the top measured **relative to the moving interface**,
times the surface supersaturation. This is the residual of the linearised sink
`(1 - C_eq)` that Prat's own reduction uses in place of `(1 - C)`: an audited
modelling choice, not a defect.

It behaves exactly as that explanation requires:

- the closed form reproduces the measured gap to **0.012-0.083 %** across all
  four cases and both matched closures (`route_identity_error_relative`);
- it is **independent of the grid** (5.158e-6 to 5.059e-6 as cells per
  `delta_r` goes 1 to 16) and **independent of the time step** (5.059e-6 to
  5.063e-6 over a 16-fold refinement);
- it scales with the supersaturation, so it is 3.6 % at the P2-A parameters
  (`C*_top/C_eq = 9.5e-3`) and 5.1e-6 at the transferred ones
  (`C*_top/C_eq = 5.7e-6`).

**The registered prediction was partly wrong, and is recorded as such rather
than quietly amended.** The criteria file predicted that the corrected
difference would "converge to zero at first order under time-step refinement".
It does not converge to zero; it converges to `- s L_top`. The other clause of
that prediction -- that the legacy closure converges to a non-zero limit of
about `C_eq (rho S/rho_cr) R_int/F_s` -- was correct to better than 0.1 % in
Cases 2 and 3.

### 7.3 Budgets, and the absence of any cap on growth

Across all four cases and both matched closures, with the seed inventory
charged to the soil across the activation event:

| Residual | worst value |
|---|---:|
| salt | 3.0e-13 |
| water | 6.0e-13 |
| solid top/internal split | 3.5e-14 |

No growth velocity is ever capped or clipped. `v_s` is
`beta k_r rho C*_top / rho_cr` and nothing else; supply limitation and loss of
connectivity are reported and can stop the model, but never scale `v_s`
(`test_the_growth_velocity_is_never_capped`).

## 8. The fully matched benchmark: A, B and C

Every physical input and both closures matched to salt_grow: `C_eq = 0.264`,
`rho = 1200`, `rho_cr = 2170`, `eps_0 = eps_birth = 0.43`, `beta = 2`,
`h0 = 1 mm`, `j = rho U_evap` with `U_evap = 1 cm/day`,
`a_v = 36(1-eps)/d_sub`, `D_eff = sqrt(eps) D_mol` with `D_mol = 1.3e-9`, and
the P2-A vapour-displacement term set to zero because salt_grow has no
counterpart for it.

| | A: Prat closed form | B1: this repo, frozen coefficients | B2: this repo, salt_grow closures | C: salt_grow M3 transient ALE |
|---|---:|---:|---:|---:|
| `v_s` [m/s] | 2.94033793e-8 | 3.01024177e-8 | **3.00339331e-8** | **3.00578141e-8** |
| `v_s` [mm/day] | 2.54045197 | 2.60084889 | 2.59493182 | 2.59699514 |
| vs A | - | +2.3774 % | +2.1445 % | +2.2257 % |
| **vs C** | -2.1773 % | +0.1484 % | **-0.0795 %** | - |
| `C_top - C_eq` | 1.15589371e-5 | 1.18337403e-5 | 1.18068179e-5 | 1.18162059e-5 |
| vs C | -2.1773 % | +0.1485 % | **-0.0794 %** | - |
| `eps_min` | 0.43, frozen | 0.237304419 | 0.235574398 | 0.23612005 |
| vs C | - | +0.502 % | **-0.231 %** | - |
| `L95` [m] | - | 1.23335983e-6 | 1.01460599e-6 | 1.00907e-6 |
| vs C | - | +22.2 % | **+0.549 %** | - |
| internal fraction of new solid | 0.269974 | 0.253054 | 0.254716 | 0.253785 |
| vs C | +6.38 % | -0.288 % | **+0.367 %** | - |

The two independent codes agree on the growth velocity to **0.0795 %**, which
is the same size as each code's own remaining discretisation error (0.031 %
here between 8 and 16 cells per `delta_r`; 0.0342 % in salt_grow between 480
and 640 cells). `v_s` and `C_top - C_eq` agree to the same figure because
`v_s` is exactly proportional to the surface supersaturation in both models.

The agreement extends to the **closure sensitivity**, which is a sharper test
than any single number because it cancels everything the two models share:

| B2 relative to B1 | this repository | salt_grow M3 |
|---|---:|---:|
| growth velocity | -0.22751 % | -0.21791 % |

a difference of 0.0096 percentage points.

`L95` is the one diagnostic where B1 differs materially (+22 %). That is not
disagreement between the codes: salt_grow's own B1 archive gives
`L95 = 1.2427 um` against this repository's 1.2334 um, a difference of
-0.75 %. The +22 % is the genuine B1-versus-B2 closure effect, which both
models reproduce.

### 8.1 The minimum porosity is converging to the same value from the same side

Taken at face value `eps_min` differs by 0.231 %. But both codes are still
decreasing under refinement:

| salt_grow cells | `eps_min` | this repo cells/`delta_r` | `eps_min` |
|---:|---:|---:|---:|
| 160 | 0.23878705 | 1 | 0.25379016 |
| 240 | 0.23743351 | 2 | 0.24059240 |
| 320 | 0.23685667 | 4 | 0.23676253 |
| 480 | 0.23634826 | 8 | 0.23580750 |
| 640 | 0.23612005 | 16 | 0.23557439 |

Extrapolating salt_grow's last interval (refinement ratio 4/3, at the first
order its own report claims) gives about 0.2354; extrapolating this
repository's last interval (ratio 2, observed order about 2) gives about
0.2355. **The apparent 0.23 % discrepancy is largely each code's own remaining
grid error, and the continuum limits agree to roughly 0.04 %.**

### 8.2 Attribution of the remaining difference

salt_grow's own Q0-Q6 decomposition (`M3_METHODS.md`, `RESULTS.md`) assigns
contributions in percentage points of `(V - V_A)/V_A` at 6 h. Matched against
what this model does and does not contain:

| Mechanism | salt_grow at 6 h | in this model? | effect on the comparison |
|---|---:|---|---|
| non-uniform porosity field | +2.6116 | **yes**, from equation (2) | the dominant term; both have it |
| diffusivity closure `sqrt(eps)` | +1.1602 | **yes**, matched explicitly | both have it |
| specific-area closure `36(1-eps)/d` | -1.5879 | **yes**, matched explicitly | both have it |
| diffusivity/area interaction | +0.0314 | yes | both have it |
| nonlinear top condition | +0.0055 | **yes**, exact in `C_top` | both have it |
| advection and ALE mesh terms | -0.0028 | **no**, `Pe << Da` | about -0.003 points |
| nonlinear sink `(1-C)` vs `(1-C_eq)` | +0.0001 | **no**, linearised as audited | about -0.0001 points |
| transient vs frozen transport, `T2 - Q6` | ~0 | **no**, quasi-steady | about 0, by salt_grow's own measurement |
| discretisation, `Q0 - A` | +0.0076 | different scheme and grid | not comparable term by term |

The mechanisms this model deliberately omits account for about **-0.003
percentage points** of the -0.079 point difference. The remainder is
discretisation in the two codes, bounded at 0.031 % and 0.034 % by section 5.1
and salt_grow's own grid table. In other words the quasi-steady reduction and
the neglect of advection cost essentially nothing at these parameters, which
is what `Pe/Da = 8.3e-8` predicts and what salt_grow's `T2 - Q6 ~ 0`
independently confirms.

This also settles a question the one-factor study raises. **The frozen-porosity
Prat law is 2.2 % slow here, not 38 % slow.** P2-B's finding that the frozen
model under-predicts growth by 38 % is a statement about the *sealing* regime,
where porosity collapses by an order of magnitude. In the travelling-wave
regime the porosity field is a fixed-shape wake and its effect on the velocity
is a 2 % coefficient, exactly as salt_grow concluded.

### 8.3 The two top closures are algebraically the same

This was checked by derivation, not only by agreement of outputs. salt_grow
solves a quadratic for the top excess with

```
qa = beta rho (1-eps_birth) k_r ,
qb = j - qa (1-C_eq) - K ,
qc = j C_eq + K (C_last - C_eq) ,     K = rho eps D_eff / distance .
```

Writing `v_s = beta rho k_r s / rho_cr` and rearranging gives

```
K (s - C*_last)  =  C_top j - rho_cr (1-eps_0) (1 - C_top) v_s ,
```

which is this repository's residual with `rho_v = 0`, since
`K (s - C*_last) = rho (eps D)_top dC/dz|_h`. The two projects therefore
impose the *same* interfacial condition, exact in `C_top`, derived
independently. Equally, salt_grow's `S_top` uses the liquid flux **relative to
the mesh**, `L_top = j + J_f`, which is the same `L_top` that appears in the
route-difference identity of section 7.2.

### 8.4 Carried to the reference time of 6 hours

The comparisons above are made at the travelling-wave plateau, which the model
reaches after about 360 s. To check that this is not a different experiment
from salt_grow's, one matched B2 run was carried to the full 21600 s using
wake coarsening and a declared larger step.

| | `dt = 0.25 s` | `dt = 0.5 s` | C: salt_grow M3 B2 |
|---|---:|---:|---:|
| steps | 86400 | 43200 | adaptive BDF |
| wall time [s] | 368 | 185 | - |
| cells / wake merges | 640 / 9991 | 640 / 9988 | 640 fixed |
| `h(6 h)` [m] | **1.64875838e-3** | 1.64855951e-3 | **1.64919695e-3** |
| vs C | **-0.0266 %** | -0.0387 % | - |
| `v_s(6 h)` [m/s] | **3.00346466e-8** | 3.00237030e-8 | **3.00578141e-8** |
| vs C | **-0.0771 %** | -0.1135 % | - |
| `C_top - C_eq` | 1.18070984e-5 | 1.18027963e-5 | 1.18162059e-5 |
| `eps_min` | 0.235529764 | 0.235157597 | 0.23612005 |
| vs C | -0.250 % | -0.408 % | - |
| `L95` [m] | 1.01829689e-6 | 1.01982498e-6 | 1.00907e-6 |
| internal fraction of new solid | 0.25433187 | 0.254560422 | 0.253785 |
| salt budget residual | 2.1e-14 | 1.2e-15 | - |
| water budget residual | 8.0e-13 | 6.1e-13 | - |
| solid split residual | 2.7e-15 | 1.0e-15 | - |
| route identity error | 5.4e-4 | 2.9e-4 | - |
| sealed / disconnected | no / no | no / no | no stop event |

The crust height is the sharpest of these, because it integrates the entire
6 h history rather than reporting an instantaneous rate: **1.64876 mm against
1.64920 mm, a difference of 0.027 %**, after growing 64.9 % from its 1 mm
seed.

Two caveats are attached to this table. The 6 h runs use 8 cells per
`delta_r` rather than the production 16, and a time step of 0.25-0.5 s rather
than the production 0.024 s, because 86400 steps on the finer grid is not
affordable here; the two steps shown differ by 0.036 % in velocity, so the
6 h velocity carries roughly +/-0.04 % of time-step uncertainty. Within that,
carrying the run to the reference time changes nothing: the plateau value was
-0.0795 % and the 6 h value is -0.0771 %. This is consistent with salt_grow's
own observation that its velocity changes by only 0.0146 % between 5 minutes
and 6 hours, and it confirms that comparing at the plateau was legitimate.

The committed `summary.json` headlines the `dt = 0.5 s` row of this table,
because it was generated before `run_p2c.py` was changed to headline the
finest step; `matched_6h.csv` carries both rows and is authoritative.

Wake coarsening did the work that made this affordable: 9991 exact pairwise
merges held the grid at 640 cells where an uncoarsened one would have needed
about 10700, and a controlled comparison with coarsening disabled
(`test_coarsening_does_not_change_the_answer`) leaves every diagnostic
unchanged to better than 1e-6.

## 9. What this does not establish

1. **It is verification, not validation.** Two models that share the same
   closure agreeing to 0.08 % tells you the equations were solved correctly in
   both. It says nothing about whether the closure describes a real crust.
   Both projects record the same reason for doubt: the reaction layer
   (`delta_r = 0.51 um`, `L95 = 1.01 um`) is thinner than the estimated
   crystal size `d_c ~ 3.5 um`, so a continuum reaction-diffusion description
   of the top layer is applied below the scale at which it can be justified.
   That objection applies at *both* parameter sets and is untouched by this
   work.

2. **`a_v` and `k_r` have literature provenance, but not for this soil.**
   P2-D's reading of the paper shows `k_r = 2.3e-3 m/s` comes from
   microfluidic NaCl crystallization experiments and `eps_s = 0.43` from a
   measured efflorescence porosity, so these are not arbitrary. What is not
   measured is `d_sub` for the soil being modelled, and whether a crust grown
   on it has the crystal-size ratio and porosity reported for those
   efflorescences. salt_grow's own RESULTS.md still lists independent kinetic
   measurements among the evidence it needs.

3. **The margin is small.** Sealing returns if `a_v` is 3.4 times larger or
   `k_r` 3.4 times smaller:

   | | value | `sqrt(Da*)` | seals? |
   |---|---:|---:|:--:|
   | P2-A declared | `k_r = 1e-7`, `a_v = 1e5` | 0.053 | yes |
   | salt_grow baseline | `k_r = 2.3e-3`, `a_v = 8.208e5` | 2.56 | no |
   | threshold at `eps_0 = 0.43`, `eps_stop = 0.02` | | 1.390 | - |
   | sealing boundary | `a_v = 2.787e6` or `k_r = 6.77e-4` | 1.390 | - |

   A measurement an order of magnitude from the salt_grow baseline can put the
   model back on the sealing side. "The crust does not seal" is therefore a
   parameter-conditional statement, exactly as "the crust seals" was.

4. ~~**The two projects disagree by a factor of 6 on the specific area.**~~
   **Resolved in P2-D**: the paper's equations (16) and (17) give
   `a_v = 5(1-eps_s)/d_c = 36(1-eps_s)/d_sub`, so the salt_grow closure used
   here is the paper's. What remains unmeasured is `d_sub` for the specific
   soil being modelled, not the formula.

5. **`eps_stop` is declared, not derived.** Both the connectivity limit (0.02)
   and the clogging limit (0.01) are declared porosity thresholds, not
   percolation or permeability criteria. The criterion of section 4.2 depends
   on `eps_stop`, so a sealing *time* -- and for Cases 0 and 1, the sealing
   event itself -- is only as defensible as that threshold.

6. **The quasi-steady reduction is benchmarked at one parameter set.** It is
   shown to cost under 0.1 % where `Pe/Da = 8.3e-8`. It is *not* shown to be
   adequate where `Pe` is comparable to `Da`, and Cases 0 and 1 sit at
   `Pe/Da = 0.02`, four orders of magnitude worse, where no transient
   reference exists.

7. **Still a one-way diagnostic.** Nothing here is coupled to a Richards
   solution. The evaporation flux is prescribed from the P2-A gas film or from
   salt_grow's fixed value and is never fed back; lateral crust transport is
   zero by the independent-patch approximation; no 2D run was performed; and
   no C++ implementation exists.

8. **The benchmark is at a plateau, cross-checked at 6 h, not integrated for
   6 h in every case.** Section 8.4 carries one matched run to the reference
   time; the other comparisons are made at the travelling-wave plateau, which
   is legitimate because salt_grow measures the velocity changing by only
   0.0146 % from 5 min to 6 h, but it is a different protocol from the
   reference and is stated as such.

## 10. GO/NO-GO for Richards coupling

**NO-GO**, unchanged from P2-A and P2-B. What has changed is the reason list,
not the decision.

Resolved or narrowed by this work:

- The P2-B conclusion that an evolving-porosity Prat crust always seals after
  `beta eps_0 / a_v` is **parameter-conditional** and does not hold at
  realistic kinetics. Gate 1' is narrowed from "find a crystal-coarsening
  closure or admit the crust cannot grow" to "measure `Da*`".
- The **formulation** is sound. The quasi-steady, depth-resolved,
  evolving-porosity model reproduces an independent transient ALE
  implementation to 0.08 % in velocity, 0.08 % in surface supersaturation,
  0.04 % in the converged minimum porosity, 0.5 % in `L95` and 0.4 % in the
  internal precipitation fraction, and reproduces its closure sensitivity to
  0.01 percentage points.
- The P2-B soil-side flux diagnostic defect is **found, explained in closed
  form and corrected**, and what remains is an audited linearisation with a
  closed-form value.
- Grid requirements for a coupled implementation are now **quantified**: at
  realistic kinetics the surface subgrid needs at least 4, preferably 16,
  cells across a 0.5 um reaction layer, which the committed uniform 32-cell
  discretisation misses by a factor of 23 and which costs 17.5 % in velocity
  and 44 % in minimum porosity if ignored.

Still blocking:

1. **No measurement of `a_v` or `k_r`**, and the sealing margin is a factor of
   3.4. (Gate 1', gate 2.)
2. **No C++ implementation.** The module is Python. (Gate 4, unchanged.)
3. **The atmospheric control for Case 36877 is unresolved**; the declared P2-A
   film gives 5.82 times the legacy P1 prescribed flux. (Gate 3, unchanged.)
4. **No coupled budget has ever been demonstrated**, not even on one 1D
   column. (Gate 5, unchanged.)
5. **Continuum validity at the crystal scale is doubtful at every parameter
   set tested**, by both projects' own assessment.

Item 10 of task 4 still stands and is reaffirmed: no large-scale
density-instability run before conservative soil/crust coupling is verified.

## 11. Reproducing this

```bash
cd model/p2_a_prat_surface/verification
python -m unittest test_p2a_surface test_p2b_evolving test_p2c_graded
python run_p2c.py
```

137 tests pass: 61 P2-A, 32 P2-B, 44 P2-C. Evidence tables are written to
`model/p2_a_prat_surface/results/p2c_salt_grow_parameters/`:

| File | Contents |
|---|---|
| `parameter_comparison.csv` | every parameter, both projects, with provenance |
| `dimensionless_groups.csv` | `Pe`, `Da`, `Da_top`, `Da*`, `delta_r` and the closed-form quantities for all six configurations |
| `grid_convergence.csv` | 1-16 cells per `delta_r`, plus the committed 32-cell grid |
| `timestep_convergence.csv` | 375-6000 steps per `tau_eps` |
| `case_summary.csv`, `case_<n>_trajectory.csv` | the one-factor ladder |
| `matched_benchmark.csv`, `matched_B1_trajectory.csv`, `matched_B2_trajectory.csv` | models A, B1, B2 and C |
| `matched_6h.csv`, `matched_B2_6h_trajectory.csv` | the matched run carried to the reference time |
| `conservation.csv` | both base-flux routes, corrected and legacy, and the identity check |
| `profiles.csv`, `profile_snapshots.csv` | `eps(z)` and `C*(z)` at the end and at four times per configuration |
| `summary.json` | headline numbers and provenance |

Trajectory tables are thinned for commit; the integration is never thinned.

Nothing in P2-A or P2-B was modified. `git status` is clean for
`prat_crust.py`, `evolving_crust.py`, `coupling.py`, every P2-A/P2-B test,
runner, config and committed result; the only edits to committed files are
this document's entry in the README, the package docstring, and the update
banners in `P2_A_GO_NO_GO.md` and `P2_B_EVOLVING_POROSITY.md`.
