# P2-E DuMuX integration overlay

This directory is an isolated patch series against the frozen P1-R2 source
tree. It does not modify P1-R2 or P2-D. Apply, in order, `p2e_main.patch`,
`p2e_initial_diagnostic_fix.patch`, `p2e_post_accept_dt_fix.patch`, and
`p2e_seed_transfer.patch` at the P1-R2 `source/` root. Place
`pratcoupling.hh` and the frozen `pratcrust.hh` beside
`appl/simpleRichardsnc/main.cc`, then build with the pinned stack recorded in
`model/p1_r2_conservative_richards/provenance/source_manifest.json`.
In that staged DuMuX copy only, change the include in `pratcoupling.hh` from
the repository-relative P2-D path to `#include "pratcrust.hh"`. The tracked
header keeps the repository-relative include so its standalone CMake target
continues to build without duplicating the frozen P2-D header.

The Box/CVFE top vertex is the soil/crust interface unknown.  Therefore the
actual DuMuX path evaluates the crust at that Newton iterate and returns its
component Neumann fluxes.  It does not create a duplicate scalar interface
unknown.  The persistent crust is read-only in `Problem::neumann`; after the
global Newton solve accepts a time step, `acceptPratStep` repeats the pure
evaluation with the accepted state and commits once.

`p2e_active_no_flux.input` is only the bounded, initially active crust
verification case.  It is not evidence for seed activation, surface-grid
convergence, reservoir supply, or production readiness.  Those remain
separate gates in the preregistered acceptance file.
