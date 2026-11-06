# P1-R2 bounded instability comparison

## Authorized scope

Exactly one new 2D run was executed after Stages B and C passed their dependent gates. It is a one-wavelength representation of Case 36877 mode 28:

- domain: 0.021428571428571432 m by 0.2 m;
- grid: 24 by 200 cells, vertical grading -0.98;
- frozen eigenprofile SHA-256: `d721e605d126e555f92b43200aa73b999ed4e2b61397edd1e98c0d9fddb36575`;
- perturbation amplitude: \(10^{-6}\);
- Mode A evaporation;
- 600 accepted steps to 136423.0637723484 s;
- no rejected steps and no solubility stop.

The comparison baseline is the frozen normalized `R24-M1-full` trajectory at the identical grid, time step, eigenprofile and amplitude. This makes the balance-law change identifiable at a fixed discretization. R24 had previously been qualified as a reduced fixed-grid reference, but it is not continuum truth. No new vertical-grid refinement was authorized or executed in Stage D.

## Mode evolution

| Diagnostic | Conservative Mode A | Frozen normalized R24 | Relative interpretation |
|---|---:|---:|---|
| Final projected \(A(t)/A(0)\) | 1.132055 | 1.129188 | new is 0.254% higher |
| Final-quarter fitted \(d\ln A/dt\) | 3.3308e-6 s^-1 | 3.3094e-6 s^-1 | new is 0.647% higher |
| Final mean surface \(x_s\) | 0.01945864 | 0.01942558 | new is 0.170% higher |
| Surface-sigma minimum time | 50021.79 s | 50703.91 s | new minimum is 682 s earlier |

At the global LSA reference time (nearest actual time 31377.3 s), both projected amplitudes are still decaying and the conservative smoothed growth rate is \(-8.27\times10^{-7}\) s^-1. At the exact mode-28 LSA reference time 45474.35 s, the conservative growth rate remains slightly negative, \(-1.08\times10^{-7}\) s^-1. Near the old numerical-detection reference (actual sampled time 50476.53 s), it is positive, \(1.38\times10^{-7}\) s^-1. The directly observed conservative surface-sigma minimum at 50021.79 s is close to, but must not be conflated with, the frozen 50394.33 s detection diagnostic.

The fixed-grid evidence therefore supports a shallow decay-to-growth transition and a slightly faster late-window growth branch after conservative reformulation. The difference is small relative to the Mode A/B background-model difference, and no experimental evidence identifies either trajectory as physically correct.

## Density and horizontal perturbation

The density field develops an increasingly sharp near-surface high-density layer while the horizontal perturbation remains small. The surface horizontal standard deviation initially decreases from \(7.2\times10^{-7}\), reaches \(6.77\times10^{-7}\) near 50022 s, and rises to \(8.06\times10^{-7}\) at the final time. This is consistent with the projected-mode trajectory and is not merely a change in the horizontal mean.

The density-field figure uses a common colour scale at initial, global-LSA, mode-28-LSA and final states. Horizontal density contrasts are visually small on the full density scale; the mode and surface-sigma diagnostics provide the more sensitive horizontal-perturbation measures.

## Conservative downward salt flux

`salt_down_kg_s` is defined positive downward. At the final state, the net flux is negative (upward) at every registered cut. At 1.03 mm below the surface:

\[
F_{s,adv}^{down}=-1.7824\times10^{-9}\ \mathrm{kg\,s^{-1}},
\]

\[
F_{s,diff}^{down}=+1.3178\times10^{-9}\ \mathrm{kg\,s^{-1}},
\]

giving net \(-4.6460\times10^{-10}\) kg s^-1. Thus downward diffusion partially offsets upward advective salt transport near the evaporating surface. At depths of approximately 20-100 mm, the final net is nearly uniform at \(-7.6552\times10^{-10}\) kg s^-1 and is almost entirely advective.

These are area-integrated values per unit out-of-plane thickness. They are conservative solver outputs, not a precipitation sink or crust-growth rate.

## Numerical status and limitations

The 2D component budgets pass: maximum total, salt and water residuals are \(9.05\times10^{-14}\), \(8.06\times10^{-16}\), and \(8.98\times10^{-14}\) kg. This establishes discrete closure on the executed grid and time-step sequence.

It does not establish:

1. spatial or temporal convergence of the mode growth rate;
2. sufficiency of R24 for all near-surface horizontal structure;
3. physical validity of the mixture-velocity or diffusion closure;
4. experimental validation of the onset time;
5. any post-precipitation behaviour.

The task required one bounded run rather than a refinement campaign. Consequently, the evidence-conclusion boundary is: **fixed-grid instability compatibility is supported; converged or physically validated instability prediction is not.**

## Figures and provenance

- `model/p1_r2_conservative_richards/results/stage_d_20261106/INSTABILITY_COMPARISON.png`
- `model/p1_r2_conservative_richards/results/stage_d_20261106/DENSITY_FIELDS.png`
- `model/p1_r2_conservative_richards/results/stage_d_20261106/DOWNWARD_SALT_FLUX.png`
- Machine-readable comparison, flux and budget files are in the same result directory.
- Raw VTU/PVD fields are retained in the remote archive recorded in `model/p1_r2_conservative_richards/provenance/source_manifest.json` and are excluded from Git.
