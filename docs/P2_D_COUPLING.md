# P2-D: preparing the first conservative Richards-Prat coupling

Date: 2026-11-07
Base commit: `f04c361169bd8271813fd384447b43d02dde9174`
Criteria registered before any run: `verification/config/acceptance_criteria_p2d.json`

## 0. What this delivers, and what it does not

Delivered:

| Task item | Status |
|---|---|
| 1. Resolve the factor-six specific-area discrepancy from the paper | **Done.** Prat equation (16) is `a_v = 5(1-eps_s)/d_c`, equation (17) is `a_v = 36(1-eps_s)/d_sub`. salt_grow was right; the `6(1-eps_0)/d` of P2-B was my own error. |
| 2. Preserve `graded_crust` and the salt_grow B2 closures | **Done.** The base condition was generalised behind a default that reproduces the committed model exactly; all 44 P2-C tests still pass unchanged. |
| 3. Port the minimum crust functionality to C++ | **Written, not compiled.** No C++ toolchain exists in this environment. The header is complete and comes with 64 reference vectors and a standalone checker, which is the gate that would transfer verification. |
| 4. Couple activity-controlled film evaporation to a 1D unsaturated Richards column | **Done**, in an isolated Python mirror of the P1-R2 closures, because P1-R2 itself cannot be built here. |
| 5. Independently verify water and NaCl budgets under the dynamic boundary | **Done.** Residuals at round-off. |
| 6. Minimal soil/crust interface with a physically consistent seed and a connected saturated crust | **Done.** The seed is formed out of the soil's own salt when the surface reaches halite saturation, and is exactly budget-neutral. |
| 7. Derive the coupled component flux conditions without prescribing both interface concentration and salt flux | **Done.** One unknown, one equation. |
| 8. Solve for compatible fluxes; crust transport flux, crust inventory change and soil component flux must close | **Done**, with one tolerance missed and explained. |
| 9. Preserve the option to stop on connectivity loss or invalid Prat assumptions | **Done.** Five declared stops, all exercised. |
| 10. No 2D runs, no change to frozen P1-R2 results | **Respected.** |

Not delivered, and these are the reasons the decision in section 12 is what
it is:

* **Nothing here ran inside P1-R2.** The coupling is demonstrated against a
  1D mirror of its closures, not against the solver itself.
* **The C++ module has never been compiled.**
* The crust remains a single column under a single patch; there is no lateral
  transport and no 2D geometry.

## 1. Item 1: the factor of six, resolved

The two sibling projects used specific areas differing by a factor of six.
The audit transcribed no formula for `a_v`, so neither was traceable. Reading
the paper directly (same PDF, SHA-256 `82E4F450...AB23AB3A`, page 5) settles
it:

```
(16)   a_v = 5 (1 - eps_s) / d_c
(17)   a_v = 36 (1 - eps_s) / d_sub        since d_c/d_sub ~ 0.14 and 5/0.14 = 35.7
```

and the diffusivity closure is `D_s,eff = D_s/tau` with `tau ~ eps_s^(-1/2)`,
i.e. `D_s,eff = sqrt(eps_s) D_s`.

So:

* **`salt_grow`'s `36(1-eps)/d_sub` is equation (17) exactly**, and its
  `sqrt(eps) D_s` is the paper's tortuosity closure. Both were right.
* **P2-B's `6(1-eps_0)/d` is not in the paper.** It was a single-cube
  surface-to-volume estimate of mine. The paper's coefficient on the crystal
  size is **5**, from numerical simulations of cube packings, not the naive 6;
  and the length scale that makes 36 appear is the *substrate grain* size, not
  the crystal size.
* The P2-B cubic corollary is corrected to `dh/d_c = beta eps_0/(5(1-eps_0))`,
  or `dh/d_sub = beta eps_0/(36(1-eps_0))` per substrate grain. Its
  qualitative conclusion is unchanged and in fact strengthened.

The paper also supplies values this repository had been treating as
undetermined: `k_r ~ 2.3e-3 m/s` from NaCl crystallization in a microfluidic
device, `eps_s = 0.43` from a measured efflorescence porosity, and
`d_c/d_sub ~ 0.14` from two independent image sources. Earlier P2-A/P2-C
statements that "the audit supplies no value for `a_v` or `k_r`" were true of
the audit and are now superseded for the paper. What is still unmeasured for
*this* project is `d_sub` for the soil being modelled, and whether a crust
grown on it resembles the efflorescences those references measured.

All of this is transcribed into `PRAT_PHYSICS_AUDIT.md` section 3.3, and the
documents that relied on the wrong reading carry corrections.

## 2. What this environment could and could not do

`model/p1_r2_conservative_richards` is DuMuX/Dune C++ and needs the pinned
stack in its `provenance/source_manifest.json`. This environment has **no C++
toolchain at all**: `g++`, `clang++`, `cl`, `cmake`, `make` and `ninja` are
all absent, which a test asserts rather than assumes
(`TestCppPort::test_no_cpp_toolchain_is_available`).

That forces a split:

* **Item 3** produces a complete, self-contained C++ module plus the means to
  verify it elsewhere. It is not verified here and is labelled as such in the
  header's first paragraph.
* **Items 4 to 9** are demonstrated in `model/p2_d_coupling/column/`, an
  isolated 1D Python mirror of P1-R2's closures. The alternative was to
  deliver no numbers at all.

The mirror is faithful where it matters and different where it has to be:

| | P1-R2 | this mirror |
|---|---|---|
| equations | total liquid mass, dissolved NaCl mass | identical |
| storage | `phi S_w rho`, `phi S_w rho w_s` | identical |
| Darcy flux | upwinded `krw`, face-averaged `rho` in gravity | identical |
| salt flux | upwinded advection + `-rho_face D_face grad(w_s)` | identical |
| `sw(pw)`, `krw`, `D_eff`, `rho(x)`, `w_s(x)` | van Genuchten / Mualem / Millington-Quirk / frozen brine law | identical expressions, identical Case 36877 parameters |
| discretisation | 2D vertex-centred box (CVFE) | 1D cell-centred two-point flux |
| time integration | DuMuX Newton with its own control | backward Euler, damped Newton, block-tridiagonal |

Nine transcription tests check the closures term by term against `main.cc`.
Discrete values must not be compared digit for digit with P1-R2 output; the
closures may.

## 3. The isolated 1D column

`column/soil.py` reproduces P1-R2's two balance equations and all five
constitutive closures exactly, on a cell-centred one-dimensional
two-point-flux grid. Conservation is not inherited from a discretisation
argument: `face_fluxes()` reconstructs every face flux from the converged
primary variables without calling the residual assembly, and `audit()`
compares the inventory change against the boundary transfer the scheme
actually applied.

Two checks establish the column before anything is coupled to it.

**The P1-R2 hydrostatic profile is an exact steady state.** With no-flux top
and bottom, the residual of the initial state is 7.06e-16 in the infinity
norm, Newton converges in zero iterations, and over 20 steps of 1000 s the
inventories drift by **0.000e+00**, against a registered gate of 1e-12.

**Prescribed component fluxes conserve.** With `F_tot = 1e-5` and
`F_s = 1e-7 kg/m2/s` for 2000 s, the liquid, salt and water residuals are
1.16e-15, 9.15e-16 and 1.05e-15, against a registered gate of 1e-10.

One correction was needed to get there. The first audit reported a 1.4e-4
water residual because the transfer was accumulated with the **start**-of-step
flux while backward Euler applies the **end**-of-step flux. `step()` now
returns the fluxes it actually applied and `run()` accumulates those; the
residual fell to 9.4e-16. A test pins this.

## 4. Item 4: water-activity-controlled evaporation

The top boundary is the P2-A gas film evaluated at the **soil surface** NaCl
mass fraction through the Pitzer water activity. No rate is prescribed, no
local weight is normalised, and the salt flux is zero: this stage has no
precipitation, as the task directs.

At the Case 36877 soil and the declared film (293.15 K, RH 0.4,
`delta = 8.3 mm`):

| | value |
|---|---:|
| initial flux | 3.0153e-5 kg/m2/s (2.61 mm/day) |
| flux when the surface reaches halite saturation | 1.7876e-5 kg/m2/s |
| decline | **40.71 %** |
| time to saturation | 3.72e5 s (103 h) |
| surface saturation `S_w` then | 0.1730 |
| cumulative evaporation | 9.303 kg/m2 |

The flux falls monotonically throughout, which a test asserts over the whole
trajectory. **The surface reaching halite saturation is a prediction of the
model, not an input**: it is what triggers the seed in item 6.

## 5. Item 5: independent budgets under the dynamic boundary

| Study | cells | `dt` [s] | liquid | salt | water |
|---|---:|---:|---:|---:|---:|
| grid | 20 | 250 | 2.88e-12 | 1.86e-12 | 2.89e-12 |
| grid | 40 | 250 | 8.65e-16 | 5.96e-16 | 8.76e-16 |
| grid | 80 | 250 | 5.62e-16 | 5.96e-16 | 5.68e-16 |
| time step | 40 | 1000 | 1.57e-16 | 1.99e-16 | 1.58e-16 |
| time step | 40 | 500 | 1.66e-15 | 1.39e-15 | 1.53e-15 |
| time step | 40 | 125 | 2.34e-15 | 2.58e-15 | 2.22e-15 |
| production | 40 | 500 | 3.39e-15 | **0.00e+00** | 3.35e-15 |

Every row is four to six orders inside the registered 1e-10 gate. The
production run's salt residual is **exactly zero**: with no salt flux at the
top and a no-flux bottom, the column's NaCl inventory is bit-for-bit
unchanged however strongly the surface concentrates.

### 5.1 Two limitations the same tables expose

**The time step is limited by the surface cell, and the limit is measured.**
`dt = 500 s` converges on 40 cells and fails on 80; `dt = 250 s` converges on
80. On 20 cells `dt = 500 s` runs for 205 steps and then fails with a singular
block once the surface cell has dried to `S_w = 0.128`, essentially the
residual saturation. The limit scales with surface-cell thickness and is
reported as a row of the evidence table rather than tuned around.

**The surface concentration is not grid-converged, and this matters.** At a
common 50000 s horizon:

| cells | surface `w/C_eq` | successive change |
|---:|---:|---:|
| 20 | 0.1078 | - |
| 40 | 0.1686 | 36.1 % |
| 80 | 0.2800 | 39.8 % |

It is diverging, not converging, and the reason is physical: evaporation
builds a surface salt boundary layer of thickness `D_eff/v_pore` with
`v_pore = j/(rho phi S_w)`, which here is about **0.33 mm**, while 80 cells
gives `dz = 2.5 mm`. The peak sharpens with every refinement.

The budgets are unaffected - conservation is a property of the scheme, not of
the resolution - but **the quantity that drives the entire coupling is the one
that is unresolved**. The time to saturation (3.72e5 s), and therefore the
seeding time, is a grid-dependent number and is not quoted as a prediction
anywhere in this document. Time-step convergence, by contrast, is clean and
first order: successive changes in the surface mass fraction of 1.34e-4,
6.72e-5 and 3.36e-5 as `dt` halves from 1000 to 125 s.

## 6. Item 6: the seed

The crust is seeded only when the soil surface actually reaches halite
saturation. Its salt and water are **removed from the soil's top cells**, so
the combined initial inventory is exactly the soil's pre-seed inventory and
the coupled budget has one unambiguous origin.

| | value |
|---|---:|
| seed height | 1.0e-5 m |
| crust cells | 156 |
| crust reaction length `delta_r` | 5.1396e-7 m |
| seed height / `delta_r` | 19.5 |
| seed NaCl | 0.0137101 kg/m2 |
| fraction of the column's NaCl | **4.91 %** |
| seed water | 0.0038 kg/m2 |
| **budget-neutrality residual** | **0.00e+00** |
| minimum crust porosity | 0.43, connected |

The seed height is declared at 10 um rather than the 100 um of P2-A and
P2-B, and the reason is recorded in the criteria file before the runs: a
100 um crust at `eps_0 = 0.43` holds 0.123 kg/m2 of NaCl, **44 % of this
column's entire inventory**, which is not a minimal perturbation of the soil.
At 10 um the crust is still 19.5 reaction lengths thick, so it remains a
continuum. An oversized seed is refused rather than silently rescaled, which
a test exercises.

## 7. Item 7: the coupled flux conditions

At the interface `z = L` (soil top = crust base) the two sides exchange one
solution-mass flux and one NaCl-mass flux. Let `w_i` be the NaCl **mass
fraction at the interface**.

**Crust side.** The quasi-steady reduced transport problem is solved with
`C(0) = w_i` instead of Prat's `C(0) = C_eq`. Writing `u = C - C_eq`, the
problem is linear with `u(0) = u_b = w_i - C_eq` and `u(h) = s`, so it
separates as

```
u(z) = u_b psi(z) + s phi(z),     phi(0)=0, phi(h)=1,  psi(0)=1, psi(h)=0
```

which costs one extra tridiagonal solve and reduces to the committed P2-C
computation exactly when `u_b = 0`. The moving-interface balances then give

```
rho U(0) = [rho eps_0 S + rho_cr (1-eps_0) - rho_v] v_s + j + R_int (1 - rho S/rho_cr)
F_s,crust(w_i) = w_i rho U(0) - rho eps_0 S D_eff dC/dz|_0
```

and `j` is itself a function of `w_i`, because the gas film is driven by the
water activity of the **crust-top** concentration that the same solve
produces. That inner dependence is closed by bisection on
`j = film(C_top(j))`, whose bracket `[0, film(0)]` is guaranteed because the
mismatch is positive at zero flux and negative at the free-water flux.

**Soil side.** The same face reconstruction every interior face uses, with the
interface half a cell above the top cell centre:

```
F_s,soil(w_i) = F_tot w_up - rho_face D_face (w_i - w_N)/(dz/2)
```

**One unknown, one equation.** The task forbids prescribing both the interface
concentration and the interface salt flux. Neither is prescribed:

```
g(w_i) = F_s,soil(w_i) - F_s,crust(w_i) = 0
```

is solved for `w_i`, and the salt flux is whatever that `w_i` implies. The
total mass flux is not a second independent condition: it is handed to the
soil as the crust's water demand, so the two sides share one solution flux by
construction.

### 7.1 Which of the two crust fluxes the soil is given

The crust offers two expressions for the same base salt flux: the **transport
route** `w_i rho U(0) + D_base`, and the **inventory route**
`d(M_solid + M_dissolved)/dt`. P2-C established that they differ by the
residual of the audited linearised sink, of order `s L_top / F_s`.

The soil is handed the **inventory route**, because that is what the crust
actually accumulates, so the coupled pair conserves NaCl exactly rather than
to the linearisation error. The transport route is then a genuinely
independent audit of the same number. Handing the soil the transport route
instead leaks exactly the linearisation residual; that was measured before the
design was fixed, and section 8 reports both.

## 8. Item 8: the compatible-flux solve, and the three routes

Sixty coupled steps of 5 s each, from the seeded state:

| | first step | last step |
|---|---:|---:|
| `t` [s] | 5 | 300 |
| `w_i / C_eq` | 0.999997332 | 0.999996836 |
| evaporation `j` [kg/m2/s] | 1.78152e-5 | 1.78152e-5 |
| `v_s` [m/s] | 3.8000e-9 | 3.8775e-9 |
| crust height [m] | 1.00190e-5 | 1.11579e-5 |
| minimum crust porosity | 0.42269 | 0.24749 |

The first result is physical, not numerical: **the crust pins the interface to
halite saturation**, to within 3.2e-6 relative, and holds it there. That is
precisely the condition Prat's reduction assumes, arrived at rather than
imposed. The mechanism is stiffness: a 2 % departure from `C_eq` changes the
crust's solution demand by two orders of magnitude and reverses its sign,
because the reduced sink acts over a half-micrometre reaction length. That
stiffness is why the search is bracketed inside the declared envelope and
widened only as needed.

The minimum porosity falls from 0.4227 towards 0.2475 at 300 s, heading for
the P2-C wake plateau `eps_0 - (1-eps_0)/sqrt(Da*) = 0.219`. It has not
arrived in 300 s and no plateau is claimed here.

### 8.1 The three routes

| Route | what it is | agreement |
|---|---|---:|
| 1, crust transport | `w_i rho U(0) + D_base` from the concentration profile | - |
| 2, crust inventory | `d(M_solid + M_dissolved)/dt`, what the soil is given | - |
| 3, soil component flux | the soil's own top-face reconstruction | - |
| | worst 2 vs 3 over 60 steps | **5.47e-10** |
| | worst 1 vs 2 over 60 steps | 9.38e-6 |
| | predicted 1 vs 2, `s L_top / F_s` | 5.81e-6 |

Routes 2 and 3 agree to 5.5e-10, against a registered gate of 1e-8: that is
the interface equation actually being solved. Route 1 is the independent
audit, and it departs by 9.4e-6 where the closed form derived in P2-C predicts
5.8e-6 - the same quantity to within a factor 1.6, which is as close as that
identity can be expected to hold when the base is not exactly at `C_eq`.
**The departure is the audited linearised sink, not a bookkeeping error.**

### 8.2 The coupled budget, and a registered tolerance missed

| | value |
|---|---:|
| combined NaCl residual | **8.18e-10** |
| combined water residual | 1.59e-11 |
| registered gate | 1e-10 |

The salt budget misses its gate by a factor of eight, and the cause is
identified rather than absorbed. Refining the time step eightfold does not
move it:

| `dt` [s] | steps | salt residual | crust height [m] |
|---:|---:|---:|---:|
| 20 | 5 | 7.550e-10 | 1.0382532e-5 |
| 10 | 10 | 7.533e-10 | 1.0383021e-5 |
| 5 | 20 | 7.525e-10 | 1.0383290e-5 |
| 2.5 | 40 | 7.517e-10 | 1.0383433e-5 |

so it is not a time-integration conservation error. It is a **one-off
inconsistency at activation**: the seed's dissolved inventory is defined with
a halite-equilibrated base, and the first coupled step immediately
re-evaluates it at the solved interface value, which is 3.2e-6 below `C_eq`.
Evaluating that single jump directly gives

```
M_d(base = w_i) - M_d(base = C_eq) = -2.223e-10 kg/m2
                                   =  7.957e-10 of the initial column NaCl
```

against the observed 8.180e-10, i.e. **97.3 % of the residual is this one
term**; the remaining 2.7 % is the accumulated Newton tolerance. The crust
step itself conserves: `dh` and `eps_min` converge cleanly at first order
under the same refinement (successive changes 4.7e-5, 2.6e-5, 1.4e-5 in
height).

This is a defect of the activation protocol, not of the coupling, and the fix
is to activate the seed at the interface concentration the first coupled solve
will find. That needs a predictor and was not attempted here.

## 9. Item 9: the declared stops

The coupled model stops rather than extrapolating when:

1. the crust loses **hydraulic connectivity** (`eps_min <= 0.02`), which is
   the primary stop, or clogs geometrically (`<= 0.01`);
2. no compatible interface concentration exists inside the declared
   halite-equilibrium envelope `C_eq (1 +/- 0.02)`;
3. the solved interface reaches the undersaturated edge of that envelope, so
   Prat's halite-equilibrated base is not defensible;
4. the top closure has no positive supersaturation root (growth arrested);
5. the crust-top or soil-surface **molality** leaves the Pitzer correlation's
   declared range of 8 mol/kg;
6. `Pe` is not small against `Da`, so the audited reduction is invalid;
7. the soil Newton solve fails.

None is a clamp. The growth velocity is always
`beta k_r rho C*_top / rho_cr` and is never capped to match available supply,
which a test asserts.

## 10. Item 3: the C++ module

`model/p2_d_coupling/cpp/pratcrust.hh` is a header-only, dependency-free
C++17 port of the minimum crust functionality: the thermophysical
correlations, the Pitzer water activity, the gas film, Prat equation (17) and
the `sqrt(eps)` tortuosity closure, the graded top-refined grid, the
quasi-steady tridiagonal solve with the generalised base condition, the exact
top closure, equation (2) porosity evolution, the moving top, inventories, the
corrected base fluxes, the gas-film fixed point and the interface root solve.

Its interface with P1-R2 is two numbers:

```cpp
const auto f = prat::solveInterface(crust, params, film, soilSaltFlux, dt);
values[conti0EqIdx] = f.totalMass;   // total liquid mass leaving the soil
values[conti1EqIdx] = f.saltMass;    // dissolved NaCl mass leaving the soil
```

matching `Problem::neumann`'s sign convention, where positive means mass
leaves the domain.

**It has never been compiled.** What replaces that, until a toolchain exists:

* `reference_vectors.txt` holds 64 quantities generated by the Python model
  that P2-A, P2-B and P2-C verified;
* `verify_vectors.cc` recomputes every one from the header and reports the
  worst relative difference against a 1e-11 gate;
* three tests check statically that the two sides compute the **same 64
  names** and share the **same declared configuration**, so they cannot drift
  apart silently, and that the stored vectors are current.

That is the gate. Running it is the first thing to do on a machine that can
build C++.

## 11. What this does not establish

1. **Nothing ran inside P1-R2.** The coupling was demonstrated against a 1D
   mirror of its closures. The mirror's discretisation differs (cell-centred
   1D two-point flux versus 2D vertex-centred box), so no discrete value here
   may be compared digit for digit with P1-R2 output.
2. **The C++ module has never been compiled**, so every number it would
   produce is unverified.
3. **The soil surface concentration is not grid-converged** (section 5.1).
   The surface boundary layer is about 0.33 mm against a 2.5 mm cell at the
   finest resolution run. Seeding time, and therefore the entire coupled
   trajectory's starting point, inherits that.
4. **The coupled window is 300 s**, long enough to establish the interface
   condition and the flux closure, far too short for the crust to reach its
   wake plateau or for the soil to respond.
5. **One registered tolerance is missed** and explained (section 8.2).
6. The crust is a single column under a single patch: no lateral transport,
   no 2D geometry, no instability question addressed.
7. Prat's continuum validity remains doubtful at the crystal scale, by both
   this project's and salt_grow's own assessment: the crust reaction layer is
   0.51 um against an estimated 3.5 um crystal size.
8. A closing budget verifies an implementation, not the physics of any
   closure.

## 12. GO/NO-GO

**GO** for the formulation: the coupled flux conditions are correctly posed
(one unknown, one equation, neither concentration nor flux prescribed), the
compatible interface is found robustly, the three flux routes close to
5.5e-10 where they should and depart by exactly the audited linearisation
where they must, and the combined system conserves NaCl and water to 8e-10
and 1.6e-11 with the one residual term identified in closed form.

**NO-GO** for running the coupling in production, and therefore still NO-GO
for any 2D instability simulation. Four blockers, in order of severity:

1. **The C++ module is uncompiled** and P1-R2 cannot be built in this
   environment. Until `verify_vectors` passes, there is no verified crust code
   for the production solver. *(This is P2-A blocker 2, unchanged.)*
2. **The soil surface boundary layer is unresolved.** At 0.33 mm it needs
   cells an order of magnitude finer than the 2.5 mm used here, i.e. of order
   2000 cells over 0.2 m in 1D, and the time step scales down with them. The
   coupling is driven by the surface concentration, so this is not a cosmetic
   convergence issue. **This is a new blocker that P2-D found.**
3. **The seeding criterion is resolution-dependent.** Because of blocker 2,
   the time at which the surface reaches halite saturation is not converged,
   and the seed is triggered by exactly that.
4. **`d_sub` is unmeasured for this soil.** Item 1 resolved where `a_v` comes
   from, and the paper supplies `k_r` and `eps_s` with experimental
   provenance, but the substrate grain size that sets `a_v` for *this* soil,
   and whether a crust grown on it resembles the measured efflorescences, are
   still open.

Gates 3 and 5 of the P2-A GO/NO-GO are unchanged: the atmospheric control for
Case 36877 is still unresolved, and no coupled budget has been demonstrated
inside the production solver.

What P2-D removes from the blocker list: the coupled flux conditions are no
longer an open design question, and the conservation properties of the
exchange are now measured rather than asserted.

## 13. Reproducing this

```bash
cd model/p2_d_coupling/verification
python -m unittest test_p2d_coupling          # 47 tests
python run_p2d.py                             # evidence tables
python generate_cpp_vectors.py                # regenerate cpp/reference_vectors.txt
```

On a machine with a compiler, the gate that would verify item 3:

```bash
cd model/p2_d_coupling/cpp
g++ -std=c++17 -O2 -o verify_vectors verify_vectors.cc && ./verify_vectors reference_vectors.txt
```

Regression across the whole surface module:

```bash
cd model/p2_a_prat_surface/verification
python -m unittest test_p2a_surface test_p2b_evolving test_p2c_graded   # 137 tests
```

184 tests pass in total. Evidence tables are written to
`model/p2_d_coupling/results/p2d_coupling/`:

| File | Contents |
|---|---|
| `film_evaporation.csv` | grid study, time-step study, the measured step limit, and the production run to solubility |
| `film_trajectory.csv` | the production film trajectory |
| `coupled_summary.csv` | seed, 60 coupled steps, the three routes and the combined budget |
| `coupled_trajectory.csv` | per-step coupled diagnostics |
| `coupled_timestep.csv` | coupled refinement from `dt = 20` to 2.5 s |
| `summary.json` | headline numbers, and the recorded absence of a C++ toolchain |

No frozen P1-R2 source, input or result is modified, and
`model/p2_a_prat_surface` keeps every committed result: the only change there
is an optional base-concentration argument whose default reproduces the
committed model exactly, which the 137 regression tests enforce.
