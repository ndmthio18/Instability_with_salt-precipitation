# Crystal-size activation with corrected top faces

Use `build_crystal.slurm`, `matrix_crystal.slurm` and
`postprocess_crystal.slurm` for the registered revision in `CRYSTAL_PLAN.md`.
The old runners explicitly specify historical thresholds for regression;
they do not define the current application default.

The application default is `P2H.ActivationHeight=3.5e-6` m. Coupled state
starts empty. Surface precipitation accumulates actual solid; activation uses
its mass-derived height, which may slightly exceed the threshold at the first
accepted crossing. There is no rounding to exactly 3.5 micrometres and no
injection of seed salt. `P2H.DonorMode=adaptive_multilayer` is explicit in the
new campaign, with 5 mm maximum depth and retained saturation >=0.99.

`epsilon_birth=0.43` and 64 ALE cells are unchanged. Computational widths are
`h * (faces[i+1]-faces[i])`, with the existing graded ALE coordinate. A cell is
not assigned a 3.5-micrometre thickness. The standalone M3 initial-height
parameter is not used as a coupled seed: coupling calls `initial(h,mixedC)`
only after conservative activation succeeds.

The frozen soil source is copied and patched by `prepare_dumux.py`. Only its
P2-H Neumann top selection changes: position AND upward face normal. On this
rectangular one-horizontal-cell mesh, lateral faces have vertical normal zero.
Boundary SCVFs are logged to `boundary_faces.csv`; postprocessing independently
checks upward normals, exact top coordinates, selected area and legacy false
positives. No evaporation law, water ledger or soil constitutive law is changed.

The 3.5-micrometre estimate is an initialization assumption supplied for this
revision. It does not calibrate the unchanged ALE microstructure/kinetic
parameters (including its existing 25-micrometre grain parameter), nucleation
time or capillary-filling timescale. Complete mixing and instantaneous filling
remain the documented numerical approximation.
