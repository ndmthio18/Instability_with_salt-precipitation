# Figure descriptions

## Figure 1 — Timing and preregistered alignment

Nonlinear harmonic signal `H23` for amplitudes `1e-6`, `1e-5`, and `1e-4` on
R24/R48/R96.  The top row uses identical physical time; the bottom row uses the
preregistered coordinate `t-t_nl`, where `t_nl` is the first of 50 consecutive
accepted steps with `H23 >= 1e-5`.  No curve-specific fitted shift is used.

## Figure 2 — `A=1e-4` three-grid transport

Surface mean concentration, conservative SCVF advective salt flux at 2 mm,
peak positive advective flux, and integrated positive advective flux for
R24/R48/R96.  The first two panels retain the signed flux trajectory; the lower
panels summarize positive downward transport.  R48 and R96 approach one another
for most 2--10 mm quantities, while the 1 mm quantity remains less consistent.

## Figure 3 — Modal energy at physical and aligned times

Energy fractions of modes m=2, 4, 6, and 10 at the common physical time 6T and
at one T after each trajectory's registered nonlinear milestone.  The figure
separates acceleration from a genuine change in modal selection.

## Figure 4 — Symmetry and timing

Maximum odd-mode energy fraction, maximum horizontal mirror-symmetry error,
registered nonlinear milestone, and resolved local-downflow onset.  Logarithmic
symmetry panels use an explicit positive floor only for display; source values
are preserved in the comparison tables.

## Figure 5 — `A=1e-5` three-grid transport

The same conservative transport and surface metrics as Figure 2, evaluated at
the intermediate amplitude without pooling it with `A=1e-4`.  Timing and surface
peaks become closer from R48 to R96, while positive advective-flux magnitudes at
1--5 mm remain grid sensitive and 10 mm has not developed positive advection by
6T.

All panels show deterministic trajectories; error bars and significance tests
would be inappropriate.  The source data are the nine verified CSV/JSON
reductions listed in `../REPRODUCE.md`.
