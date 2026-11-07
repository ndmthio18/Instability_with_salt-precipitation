# P2-F isolated DuMuX overlay

Base commit: `02c352f8c0f64f9dff952121b78ae7a61259a2e6`.

Apply `p2f_main.patch` at the root of the frozen P1-R2 `source/` tree and
place `../cpp/analyticalgrowth.hh` beside `appl/simpleRichardsnc/main.cc`.
The overlay adds the `p2f_analytical` boundary mode. It does not modify the
frozen P1-R2 tree in this repository and has no dependency on P2-D/P2-E,
`solveProfile`, a crust mesh, or an internal reaction source.

The benchmark is geometrically one-dimensional: one horizontal cell and a
graded vertical grid. Grid, step, seed, and supply variants are set with DuMuX
command-line overrides by `hpcc/run_coupled.slurm`.
