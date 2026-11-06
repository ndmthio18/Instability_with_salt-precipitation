# P2-A surface subgrid: physical film evaporation and Prat (2026) crust growth

## 1. Scope and source boundary

P2-A adds a **standalone external surface module** next to the P1-R2
conservative Richards-salt solver. It does not modify P1-R2, does not modify
any frozen source, and runs no DuMuX simulation. Its isolated location is
`model/p2_a_prat_surface/`.

The scientific setting is unchanged: unsaturated homogeneous soil, one
volumetric Richards domain, no precipitation inside the original soil pores,
and an external porous salt structure represented by a surface subgrid.

What P2-A delivers:

1. the heterogeneity project's physical gas-film evaporation formulation,
   reused without calibration and without weight normalisation;
2. the Prat (2026) growth law with its finite-`Da` form, its large-`Da` limit,
   crust porosity, crust height, time integration and water/salt inventories;
3. an explicitly declared seed/precipitation activation event;
4. the soil-side water and NaCl fluxes reconstructed from Prat's
   moving-interface balances and the reduced crust transport solution;
5. a quantified comparison of the two evaporation-interface choices.

What P2-A does **not** deliver: any coupled DuMuX run, any density-instability
result, any experimental validation, and any claim that the Prat closure is
valid at the declared parameters. Section 7 lists the open items.

## 2. Reused sources, traced

| Quantity | Source file | Note |
|---|---|---|
| NaCl water activity (Pitzer osmotic coefficient) | heterogeneity `source/appl/water_activity.hh`, `naclWaterActivity` | Rose et al. (2008) Eqs. A15-A16 with Pitzer and Mayorga (1973) parameters at 298.15 K; use at 293.15 K is the upstream file's own stated approximation and is inherited unchanged |
| Gas-film evaporation flux | heterogeneity `source/appl/problem.hh`, `filmWaterFlux_` | mass-fraction-gradient film closure, see section 3 |
| Pure-water saturation pressure | DuMuX `dumux/material/components/iapws/region4.hh` | IAPWS-IF97 Region 4, ten coefficients transcribed |
| H2O-air gas diffusion coefficient | DuMuX `dumux/material/binarycoefficients/h2o_air.hh` | `2.13e-5 (1e5/p) (T/273.15)^1.8` |
| Halite solid density | DuMuX `dumux/material/components/nacl.hh` | 2165 kg/m3 |
| Gas constant, air molar mass | DuMuX `constants.hh`, `air.hh` | 8.314472, 0.02896 |
| H2O/NaCl molar masses, brine density law, solubility mole fraction | frozen Case 36877 input and P1-R2 `main.cc::rhoW` | reproduces the registered P1-R2 conversions |
| Prat equations (1)-(3), (5)-(8), (13)-(15), (20), (28), (29)/(32), (36) | `docs/PRAT_PHYSICS_AUDIT.md` | the audit was written against the publisher PDF; no equation is reconstructed from memory |

**Deliberately not reused**: `normalizedEvaporationWeights` from
`water_activity.hh`. P2-A item 4 forbids normalising local evaporation weights
to a fixed total rate, so the helper is excluded by design, and a test asserts
that it has not reappeared.

## 3. Evaporation: physical inputs only

The water-component mass flux through the gas film is

```
j = rho_ref * D_v / delta * ( X_w(y_s) - X_w(y_ref) )        [kg H2O m^-2 s^-1]
y_s   = a_w * p_sat(T) / p_atm
y_ref = RH  * p_sat(T) / p_atm
X_w(y) = M_H2O y / ( M_H2O y + M_air (1-y) )
```

with `rho_ref` the gas density at the ambient composition `y_ref`. The inputs
are exactly the five physical quantities required by P2-A item 2: temperature,
ambient relative humidity, gas-film thickness, vapour diffusion coefficient and
atmospheric pressure.

Three properties follow and are enforced by tests:

* **No calibration.** `delta` is a declared input. The module can *report* the
  film thickness that would reproduce a given flux
  (`equivalent_thickness_for_flux`), but that value is never fed back. At the
  declared reference the film gives 5.82 times the legacy P1 prescribed flux,
  and reproducing the legacy flux would require `delta = 4.83e-2 m`. This is
  reported as a model difference, not resolved by tuning.
* **No normalisation.** The flux is evaluated pointwise from the local water
  activity. The total evaporation rate is free to change with salinity, and it
  does: between fresh water and halite equilibrium it falls by 41.6 per cent.
* **Signed.** Condensation is returned as a negative flux, following the
  upstream comment, rather than being clipped to zero.

The closure is Fickian in the **mass fraction** with a mass-averaged gas
density. That is the same reference frame as the P1-R2 liquid closure
`j_s = -rho_l D_eff grad(w_s)`, so the gas and liquid sides compose without a
frame change. The gas phase is treated as an ideal binary H2O/air mixture; this
is an explicit assumption of the standalone module, adopted so the film state
can be evaluated without a DuMuX fluid system.

Before any surface salt exists, the water activity is computed from the **local
Richards soil-surface NaCl concentration** (P2-A item 3). The primary variable
is the mole fraction `x_s`; it is converted to a mass fraction with the frozen
P1-R2 closure and then to molality as `m = x_s / ((1-x_s) M_H2O)`, exactly as
in the upstream heterogeneity top-face loop.

## 4. The Prat crust module

### 4.1 Imported verbatim

Dimensionless groups (audit section 3):

```
Pe      = (j/rho) h / (eps_c D_s_eff)
Da      = h^2 (1-C_eq) a_v k_r / (eps_c D_s_eff)
Da_top  = beta (1-eps_c) (1-C_eq) h k_r / (eps_c D_s_eff)
Da*     = beta^2 (1-eps_c)^2 (1-C_eq) k_r / (eps_c D_s_eff a_v)
```

Finite-`Da` growth law, equation (28), and its large-`Da` limit, equations
(29)/(32):

```
v_s = C_eq j / ( rho_cr (1-eps_c)(1-C_eq) ) * f
f      = sqrt(Da*) tanh(sqrt(Da)) / ( 1 + sqrt(Da*) tanh(sqrt(Da)) )
f_inf  = sqrt(Da*) / ( 1 + sqrt(Da*) )
```

Moving-interface balances, equations (6) and (7), including the `-rho_v dh/dt`
sweep term, which is retained and reported rather than hidden in a sign
convention.

Internal reaction: equations (1)-(2) give the solid production rate density
`r_i = a_v rho k_r (C - C_eq)`, while the transport equation for `C`,
equation (5), carries an extra `(1-C)` factor. The two are kept distinct, so
`(1-C_eq)` appears in `Da` but not in `r_i`.

### 4.2 Derived here, and labelled as such

The audit supplies the boundary conditions `C*(0) = 0`, `C*(1) = C*_top`, the
quasi-steady `Pe << Da` reduction, and equation (8). It does not print the
reduced solution. The module therefore derives:

```
C*(z*) = C*_top sinh(sqrt(Da) z*) / sinh(sqrt(Da))
C*_top = C_eq Pe tanh(sqrt(Da)) / ( sqrt(Da) (1 + sqrt(Da*) tanh(sqrt(Da))) )
```

`C*_top` is fixed uniquely by requiring that equation (8) evaluated on this
profile reproduce equation (28). Consistency is then checked three ways, and
all three agree to round-off (relative difference below 1e-12 over
`h` from 1e-5 to 2e-2 m):

1. equation (28) directly;
2. equation (8) with the analytic top gradient;
3. the outer-surface kinetic relation `v_s = beta k_r rho C*_top / rho_cr`,
   which is the form of equation (20) compatible with equation (28). Note that
   the variant carrying an extra `(1-C_eq)` is *not* compatible; the module
   records which form it uses.

A second derived result, which Prat (2026) does not state and which the audit
explicitly warns must not be assumed:

```
internal reaction / top skeleton rate  ->  1 / sqrt(Da*)      (Da >> 1)
hence  f_inf = 1 / ( 1 + internal/top )
```

So in the large-`Da` limit the Prat efficiency **is** the top-growth share of
total precipitation, and `1 - f_inf` is the internal share. This is derived
from the integrated analytic profile, not postulated, and it is verified
numerically to 1e-6 at `Da >= 1e4`. It is *not* a partition of the incoming
soil salt flux, because that flux also feeds dissolved storage; the audit's
prohibition on writing `dM_s,i/dt = (1-f) F_s,soil` still stands.

### 4.3 Activation is separate from growth

`SeedEvent` carries the declared activation rule. Until it fires:

* `h = 0`, and every Prat quantity raises `CrustModelDomainExit` rather than
  returning a number, because at `h = 0` the continuum domain, the
  diffusion-time argument and the boundary-layer argument all fail;
* the patch reports the ordinary soil evaporation boundary with **zero** net
  salt flux, i.e. exactly the P1-R2 top condition.

The trigger is a thermodynamic *eligibility* threshold (soil surface reaching
halite equilibrium). Eligibility time and activation time are recorded
separately so that a solubility crossing is never reported as a nucleation
prediction. The seed height and porosity are declared modelling inputs.

The seed carries real mass. That mass is charged to the soil as one
instantaneous soil-to-crust transfer and recorded in `seed_salt`/`seed_water`,
because creating crust inventory from nothing would reproduce exactly the
defect P2-A item 8 forbids.

## 5. Conservative soil/crust exchange

### 5.1 State and rates

The crust inventories are **states**, advanced by the same rates that are
reported as interface fluxes:

```
dh/dt          = v_s
dM_s,solid/dt  = (1-eps_c) rho_cr v_s + R_int
dM_s,d/dt      = eps_c S_c rho C_top v_s
dM_w,c/dt      = eps_c S_c rho (1-C_top) v_s
R_int          = a_v rho k_r h C*_top tanh(sqrt(Da)/2) / sqrt(Da)
```

so that the two identities of `docs/INTERFACE_CONSERVATION.md` hold by
construction:

```
F_s,soil = d( M_s,solid + M_s,d )/dt
F_w,soil = j - rho_v v_s + dM_w,c/dt
```

No arbitrary salt sink is introduced, and `F_s,soil` is strictly positive
whenever solid accumulates.

### 5.2 Two routes to the base salt flux, and their exact gap

The base flux is also reconstructed independently as

```
rho U(h) = ( rho eps_c S_c + rho_cr (1-eps_c) - rho_v ) v_s + j      (eq 6)
rho U(0) = rho U(h) + R_int                  (quasi-steady solution mass)
F_s,base = rho C_eq U(0) - rho eps_c S_c D_s_eff dC/dz|_0+
```

The two routes do not agree identically, and the gap has a closed form:

```
F_s,base - F_s,balance = -( eps_c S_c rho C*_top + C_eq rho_v ) v_s
```

The leading term is the dissolved salt in the newly created pore volume *in
excess of equilibrium*; the second is the equation (6) vapour sweep. The gap
exists because equation (7) is exact in `C_top` while equation (8), on which
equation (28) is built, sets `C ~ C_eq` at the top. It is the
top-linearisation residual of the published derivation, not a numerical defect,
and it vanishes as `C*_top -> 0`. Tests confirm the closed form to 1e-12 of the
flux scale. The reported `F_s,soil` is the **balance route**, because that is
the one that matches the tracked inventories term by term.

### 5.3 Crust-top concentration is not assumed equal to the soil top

`C_top = C_eq + C*_top` with `C*_top > 0` strictly. P2-A item 8's prohibition
is satisfied: the module never assumes `C_top = C_b`.

### 5.4 What the quasi-steady reduction drops

The reduced solution retains the new-pore-volume part of `dM_s,d/dt` but drops
the change of the dissolved *profile shape*. That omission is reported, not
hidden: every step records `dissolved_salt_quasi_steady_kg_m2` (the profile
diagnostic) next to the integrated state and their difference. Over the
verification trajectories the gap stays between 8.9e-5 and 1.5e-3 of the
dissolved inventory.

### 5.5 Supply limitation is a diagnostic, never a cap

`supply_limitation_diagnostic` compares the reconstructed demand against an
available soil flux and reports the ratio. It never clips `v_s`, and with
`raise_on_shortfall` it raises `CrustModelDomainExit`. This follows
`INTERFACE_CONSERVATION.md` section 7.1: a pointwise
`min(Prat demand, available flux)` changes `v_s` without updating `C_b`,
`C_top`, `U`, the dissolved storage, the new pore water and the internal
reaction, and is therefore not conservation-consistent.

Loss of liquid connectivity is handled the same way: below the declared
`connectivity_threshold` the saturated closure is disabled rather than
extrapolated.

## 6. Declared parameters that the audit could not supply

`a_v` and `k_r` have no value in `PRAT_PHYSICS_AUDIT.md`. They are declared
inputs in `verification/config/acceptance_criteria.json`, registered before the
runs, and swept over three to four decades in the evidence. Absolute growth
rates are therefore **parameter-conditional**. The verification targets are
identities, limits and convergence orders, which are parameter-independent.

`eps_c` and the seed height are likewise declared. `D_s_eff` uses a declared
tortuosity closure `D_mol sqrt(eps_c)` on the frozen Case 36877
`SaltDiffusivity`; it is not calibrated.

## 7. Open items and known tensions

1. **Constant porosity versus equation (2).** Equation (2) forces the crust
   porosity to change wherever `C != C_eq`, yet the reduced solution holds
   `eps_c` fixed (audit section 3, reduction 3). The module reports
   `porosity_consumption_timescale_s = eps_c rho_cr h / R_int` and flags
   `constant_porosity_defensible` when it exceeds ten growth times. **At the
   declared reference parameters the flag is false at every height examined**
   (3 of 61 heights pass). Because internal/top = `1/sqrt(Da*)`, a
   self-consistent constant-porosity crust needs `Da*` of order one, which is
   exactly the regime in which `f_inf -> 1` and the efficiency factor carries
   no information. This tension is a property of the published closure, not of
   this implementation, and it must be resolved before the crust feeds back
   into a soil simulation.
2. **The module is Python.** It is a reference implementation and a verified
   specification, not yet a DuMuX surface subgrid. No C++ compiler was
   available in this environment, so a C++ port is deferred and is a
   prerequisite for any coupled run.
3. **No feedback to the soil.** P2-A item 7 is a one-way diagnostic by
   instruction. Nothing here has been run against a live Richards solution.
4. **Film thickness.** The declared `delta` is inherited from the heterogeneity
   cylinder experiment, whose geometry differs from the Case 36877 column. It
   gives 5.82 times the legacy P1 flux. Choosing the right atmospheric control
   for this problem is an open physical question, listed as gate 1 of the
   P1-R2 GO/NO-GO.
5. **Lateral crust transport** is zero by the independent-patch approximation,
   and the horizontal variability that the density-instability study is about
   is exactly where that approximation is weakest.
