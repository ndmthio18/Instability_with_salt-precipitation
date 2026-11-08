# P2-H reactive-area sensitivity

Branch: `research/p2-h-av-sensitivity`; base: `301097e3d2de3541bf63741cd729bd5cda22329a`.

The latest verified surface-precipitation closure is ported as an explicit mode:

```text
P2H.PrecipitationMode = implicit_surface
```

Historical `projection` behavior remains available. The sensitivity matrix uses
the corrected mode so that pre-activation salt storage does not inherit the known
top-BOX splitting error.

## Parameters and code boundary

`P2H.ReactionCoefficient` defaults to the historical `2.3e-3 m/s` and controls
both the M3 internal reaction and the M3 direct top-growth kinetic law. The task
matrix explicitly sets `1e-3 m/s`.

`P2H.AvMultiplier` defaults to one and enters only

```text
a_v(epsilon) = AvMultiplier * 36 * (1 - epsilon) / d_sub.
```

It does not enter the top root or top kinetic velocity directly. It can change
top growth indirectly through the coupled concentration and porosity fields.
`d_sub=25e-6 m`, `epsilon_birth=0.43`, the B2 diffusivity, M3 conservative
operator, gas film, donor rules and all physical bounds are unchanged.

The physical comparisons use `h_m` and `top_growth_velocity_m_s`, both evaluated
by the ALE model with `P2H.ReactionCoefficient`. The historical `h_prat_m` column
is an older P2-G reference diagnostic and is excluded from this sensitivity
analysis; it does not drive the coupled state.

The standalone `av_parameter_test.cc` checks linear internal-rate scaling,
direct top-law independence and signed dissolution. The frozen M3 regression is
run with the historical defaults. `av_mechanism.cc` is a separate prescribed-state
mechanism check; it is not classified as a coupled DuMux result.

## Cumulative reaction convention

The conservative BDF net internal-solid ledger remains authoritative. Positive
precipitation and dissolution are reported separately as optional BDF state
ledgers integrated with the same accepted internal steps:

```text
P(t) = integral sum_i max(R_i, 0) V_i dt
D(t) = integral sum_i max(-R_i, 0) V_i dt.
```

`P-D-net_internal_solid` is retained as an explicit ledger-closure diagnostic.
The optional two-state extension is disabled by default, so the frozen historical
M3 state dimension and regression remain unchanged. Coupled sensitivity cases
enable it explicitly. Rejected coupled trials never update these accepted-state
accumulators.

`L95_bottom` integrates physical cell volumes from the soil-crust interface.
`L95_shortest` searches mass-coordinate interval endpoints and linearly resolves
the fractional boundary cells. Both are undefined when the positive rate is
numerically negligible. Cell indices are never interpreted as material points.

## Registered stages

- `build_av.slurm`: corrected DuMux build, parameter unit tests, frozen M3 full
  regression, standalone matched-state mechanism check and short coupled smoke.
- `matrix_av.slurm primary_exact`: 80/160 soil cells, 64 ALE cells, three multipliers.
- `matrix_av.slurm dt_half_exact`: the same six cases with halved soil coupling step.
- `matrix_av.slurm ale_refine_exact`: selected 80-cell cases with 128 ALE cells, used
  only after the 64-cell reaction-zone resolution audit.
- `matrix_av.slurm ale256_exact`: the selected 80-cell, `10av` case with 256 ALE
  cells after the 128-cell local-porosity gate remained unresolved.

All coupled cases start with zero crust and zero solid inventory, activate at
3.5 micrometres from actual surface precipitation, and use conservative adaptive
donation with a 5 mm search depth and retained saturation of 0.99.
