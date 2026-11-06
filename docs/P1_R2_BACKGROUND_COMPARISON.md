# P1-R2 Case 36877 background comparison

## Comparison design

The frozen normalized Richards solver and the conservative solver were run as one-horizontal-cell, 200-cell vertically graded columns over 136423.0637723484 s. Permeability, porosity, initial pressure and salinity, retention and relative-permeability parameters, density and viscosity laws, gravity, height and bottom reservoir were held fixed. No parameter was tuned to reproduce the frozen profile.

The conservative solver was run as:

- Mode A: fixed H2O mass demand referenced to 998.203 kg m^-3;
- Mode B: fixed equivalent liquid volumetric flux using the current face density;
- Mode A with half the base time step.

Profiles were compared at the 13 physical output times shared exactly by all four runs. The older solution is a normalized-PDE reference, not a physical-component-conservation reference.

## Final matched-time differences

At the final shared output, 136423 s:

| Comparison | \(p_l\) relative L2 | \(S_l\) relative L2 | \(x_s\) relative L2 | \(w_s\) relative L2 | \(\rho\) relative L2 | \(q_z\) relative L2 | max absolute \(q_z\) difference (m s^-1) |
|---|---:|---:|---:|---:|---:|---:|---:|
| Mode A vs frozen normalized | 3.78e-9 | 3.84e-8 | 8.10e-4 | 7.67e-4 | 8.91e-6 | 6.78e-3 | 1.16e-10 |
| Mode B vs Mode A | 1.64e-6 | 1.73e-5 | 2.43e-2 | 2.30e-2 | 2.67e-4 | 4.37e-2 | 1.35e-10 |
| Mode A \(\Delta t/2\) vs base | 9.51e-12 | 9.52e-11 | 6.12e-6 | 5.87e-6 | 6.75e-8 | 8.30e-8 | 1.01e-15 |

The largest relative Mode A/frozen discrepancy is in vertical velocity, but the absolute difference is only \(1.16\times10^{-10}\) m s^-1 because the reference velocity is small. Pressure and saturation remain nearly identical. Composition differs at approximately \(8\times10^{-4}\) relative, which is expected because the conservative salt storage and diffusion law are not algebraically equivalent to \(\theta x_s\) and \(D\nabla x_s\).

Mode A/B differences are much larger than time-step effects: the final mass-fraction profile differs by 2.30% relative L2, versus \(5.87\times10^{-6}\) for time-step halving. Mode B removes mass in proportion to the salinity-dependent face density, whereas Mode A imposes a constant pure-water-referenced mass demand. As surface density increases, the two boundary conditions diverge. This is a model-definition effect and must not be absorbed by retuning diffusivity or evaporation.

## Inventories and budgets

The conservative runs independently close total, salt and reconstructed-water budgets. Maximum full-window absolute residuals are below \(1.38\times10^{-13}\) kg for total mass and \(1.60\times10^{-15}\) kg for salt. The bottom reservoir can supply liquid while the top evaporates, so total in-domain mass is not required to decrease monotonically; closure is assessed using both independently reconstructed boundary transfers.

Physical inventories reconstructed from the frozen normalized fields are diagnostic conversions only. They do not turn the old residual into a physical mass balance, so no physical-budget PASS is assigned to the frozen run.

## Numerical interpretation

Time-step halving shows that the 1D conservative background is insensitive to this temporal refinement at the reported precision. The acceptance file pre-registered component-budget tolerances but did not pre-register a separate profile-error threshold for V5. Therefore the profile comparisons are reported quantitatively rather than retroactively labelled as passing a formal convergence tolerance.

The Stage C evidence was sufficient to authorize the single bounded Stage D comparison because:

1. all component budgets passed their registered gates;
2. the forced retry was correctly represented in accepted-step accounting;
3. the half-step solution was close to the base step;
4. Mode A/frozen differences were bounded and mechanistically attributable to the changed balance law;
5. no parameter tuning or production campaign was introduced.

This authorization is not a claim of spatial convergence or physical validation.

## Figures and data

- `model/p1_r2_conservative_richards/results/stage_c_20261106/BACKGROUND_FINAL_PROFILES.png`
- `model/p1_r2_conservative_richards/results/stage_c_20261106/BACKGROUND_DIFFERENCES.png`
- `model/p1_r2_conservative_richards/results/stage_c_20261106/BACKGROUND_TRAJECTORIES.png`
- Machine-readable profiles, velocities, inventories, comparisons and summary are in the same result directory.
