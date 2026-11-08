# Registered xNaCl=0.08, 3.5-micrometre activation test

- Branch: `research/p2-h-reservoir-x008-h3p5`.
- Parent: `b2e1797a54a68f56ea6e8c49c1ad72360f847ba9`.
- Controlled change: `P2H.ActivationHeight` changes from `1e-7 m` to
  `3.5e-6 m`.
- Initial soil and bottom reservoir remain `x_NaCl=0.08`.
- Surface solid and crust height start at zero. ALE activation must wait for
  actual conservatively accumulated surface solid to cross 3.5 micrometres.
- Unchanged: 80 soil cells, 64 ALE cells, 0.5 s maximum timestep, 50 h
  requested horizon, adaptive 5 mm donor depth, `Sw_min=0.99`,
  `epsilon_birth=0.43`, gas film, Prat/M3/ALE equations, reaction and
  diffusivity closures, boundary reservoir and conservation tolerances.

Primary causal comparison: the earlier 0.1-micrometre run stopped at
20580.510687 s because its soil-side ALE cell reached `epsilon=0.95` after net
internal dissolution. The new run tests whether waiting for a crystal-scale
actual solid inventory prevents or delays that failure.

Acceptance gates:

1. Registered and accepted activation height is at least `3.5e-6 m`.
2. Activation transfer residual is at most `1e-12` relative.
3. Global water and NaCl residuals are at most `1e-8` relative.
4. The run reaches 180000 s with a nontrivial accepted ALE trajectory.

A completed Slurm job without ALE activation is not a coupled-model pass. A
failed run remains censored at its last accepted physical time and will be
diagnosed; no porosity clipping or tolerance relaxation is permitted.
