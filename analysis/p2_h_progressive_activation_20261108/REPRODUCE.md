# Reproduction instructions

## Repository and HPCC scope

- Branch: `research/p2-h-progressive-crust-activation`
- Base commit: `77e9da996d6d784743ed2974ab64be49e5b26f19`
- HPCC root: `/WORK/hulm_work/zhixin/p2_h_progressive_activation_20261108`
- All compilation, numerical simulation and plotting run through Slurm compute nodes.
- Login-node actions are limited to staging, submission and read-only monitoring.

## Registered sequence

```bash
cd /WORK/hulm_work/zhixin/p2_h_progressive_activation_20261108
build=$(sbatch --parsable progressive_build.slurm)
pilot=$(sbatch --parsable --dependency=afterok:${build} progressive_pilot.slurm)
# Inspect the pilot initialization and accepted ALE segment before continuing.
sbatch progressive_matrix.slurm
sbatch progressive_ale.slurm
# After all arrays terminate:
sbatch progressive_postprocess.slurm
```

The run registry is `model/p2_h_ale_coupling/run_progressive.py`.  It freezes
`x_NaCl=0.08` (mole fraction), reservoir bottom, `WaterLevel=0.01 m`,
`kr=1e-3 m/s`, `alpha_av=1`, `epsilon_birth=0.43`, activation height
`3.5e-6 m`, and ALE 64 cells.  The progressive matrix is the Cartesian product
of soil cells `80/160/320` and maximum time step `0.5/0.25/0.125 s`; old
`adaptive_multilayer` and `top_only` controls use `Nsoil=80`, `dt=0.5 s`.

ALE 128/256 checks are registered separately and are interpreted only after
the initializer's mass, volume and state gates pass.  Raw histories remain on
HPCC; Git receives compact tables, figures, logs, checksums and failure records.
