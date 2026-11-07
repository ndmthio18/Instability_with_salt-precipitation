# P2-F analytical surface growth

This isolated module implements the reaction-free, constant-porosity surface
growth law derived in `docs/P2_F_MODEL_FORMULATION.md`. It contains Python and
C++ implementations, a P1-R2 DuMuX integration patch, HPCC Slurm scripts,
pre-registered acceptance criteria and lightweight verification results.

No frozen P1-R2 or P2-A--E source is modified. The production path has no
crust mesh, internal reaction solve, dynamic porosity or `solveProfile`
dependency.

Scientific decision: see `GO_NO_GO.md` (NO-GO for 2D coupling).
