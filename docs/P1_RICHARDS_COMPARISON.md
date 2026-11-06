# P1 comparison with the frozen Richards--salt baseline

## Matched reference

Case 36877 is the comparison case. Both formulations use the 0.60 m by 0.20 m domain, homogeneous material parameters, water-level depth 0.6307756333 m, initial NaCl mole fraction 0.0036, evaporation velocity \(3.065544359\times10^{-9}\) m s\(^{-1}\), salt molecular diffusivity \(1.5\times10^{-9}\) m\(^2\) s\(^{-1}\), gravity 9.8 m s\(^{-2}\), and the frozen mode-28 eigenprofile at amplitude \(10^{-6}\).

The reference timing labels remain distinct:

- global LSA neutral crossing: 31409.158 s, critical mode 5;
- observable mode-28 LSA crossing: 45474.355 s;
- frozen Richards numerical detection: 50394.334 s.

The new-model amplitude minimum reported below is another numerical diagnostic. It is not relabelled as any of those three times.

## Quantitative comparison

| Observable | Frozen Richards reference | P1 compositional result |
|---|---:|---:|
| Domain / production grid | 0.60 m × 0.20 m / 600 × 200 | same formal input; executed benchmark resolution reported below |
| Conservation variable | normalized \(\phi S_wx_s\) | physical H2O/air/NaCl molar inventories |
| Surface NaCl flux | zero | exactly zero prescribed |
| Global LSA crossing | 31409.158 s, mode 5 | not recomputed by P1 |
| Mode-28 LSA crossing | 45474.355 s | independent reference line only |
| Numerical detection | 50394.334 s | not detected by 60000 s; mode-28 minimum remained at the endpoint |
| Post-minimum log-amplitude slope | reference LSA value not recovered in P1 audit inputs | undefined because no post-minimum interval exists |
| Final surface salinity | reference field series exists | B0 0.008018127631; B1 0.008018127628 |
| Final liquid saturation at top | reference field series exists | B0 0.258141086768; B1 0.258141086768 |
| Maximum absolute NaCl budget residual | normalized audit | B0 \(2.11\times10^{-12}\); B1 \(1.32\times10^{-12}\) mol m\(^{-2}\) |

## Constitutive differences that prevent bitwise matching

1. The Richards code advances a water equation and a normalized dissolved-salt variable. P1 advances physical molar H2O, air and NaCl component balances across both phases.
2. Richards uses its reference viscosity and quadratic density law. P1 retains upstream `BrineAir` water activity, salinity/pressure-dependent liquid density and viscosity, gas properties, Henry/Raoult equilibrium and dissolved-air transport.
3. P1 uses upstream mass-averaged Fick diffusion with conversion to molar component equations. Its effective coefficient contains Millington--Quirk porosity/saturation scaling.
4. P1 phase-state switching changes the primary-variable set. Richards has no gas-component conservation or phase-composition equilibrium.
5. TPFA fluxes, nonlinear scaling and phase upwinding differ from the frozen solver.

Accordingly, direct comparison of a Richards \(\phi S_wx_s\) integral with P1 mol m\(^{-2}\) is dimensionally invalid. A justified transformation would require the local phase molar density and the same control volume at each saved state.

## Evidence limit

The executed comparison used the full 0.60 m × 0.20 m physical domain but an
exploratory 280 × 100 grid, fixed 227.372 s maximum time step and a 0--60000 s
window. The grid provides ten cells per imposed mode-28 wavelength; it is not
the formal 600 × 200 production grid. Both B0 and B1 ended normally after 264
accepted steps. B0 remained symmetric with final mode amplitude
\(3.56\times10^{-16}\). B1 decreased monotonically from
\(7.61\times10^{-7}\) to \(4.70\times10^{-7}\), including after the LSA and
Richards reference times. Therefore the compositional benchmark did not
reproduce the frozen Richards instability within the registered window.

The registered one-factor checks further show that the exploratory result is
time-step insensitive in its terminal top salinity (0.00773% relative change)
and mode amplitude (0.0282%), and horizontally insensitive for those terminal
observables (\(1.34\times10^{-11}\) and 1.53%). It is not vertically insensitive:
reducing 100 to 50 vertical cells changes terminal top salinity by 26.96% and
mode amplitude by 48.61%. Moreover, every mode-amplitude minimum remains at
the endpoint, so none of the pairs has a post-minimum growth rate to compare.
These are bounded sensitivity findings, not an asymptotic convergence claim.

The comparison tests whether the same pre-precipitation mechanism is present under matched primary inputs. It does not validate either model against an experiment and does not license surface precipitation or Prat growth.
