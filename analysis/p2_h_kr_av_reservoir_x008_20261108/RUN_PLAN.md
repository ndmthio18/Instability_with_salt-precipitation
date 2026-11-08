# P2-H reservoir xNaCl=0.08 registered run

- Branch: `research/p2-h-kr-av-reservoir-x008`.
- Parent: `b6fed0087531aa7015b5ba8b435076a7b99332d3`.
- Only declared physical change: both the accepted initial soil field and the
  bottom reservoir Dirichlet condition use `Problem.SaltMoleFraction=0.08`.
- Unchanged: reservoir water level, soil/ALE grids, gas film, activation
  threshold, zero initial solid/crust, donor depth, `Sw_min=0.99`, kinetics,
  constitutive laws, timestep and 50 h requested horizon.
- The actual accepted concentration is verified from `reservoir_config.csv`;
  the configuration is not inferred from the command line alone.

The registered dependency gate is:

1. Repeat the short hydrostatic/evaporating/no-flux verification at x=0.08.
2. Run one reference 80-soil/64-ALE, `dt_max=0.5 s`, 50 h pilot.
3. Continue only if the pilot records real surface solid, a conservative
   accepted ALE activation, and at least 100 s of nontrivial ALE evolution.
4. If that gate passes, run ALE128 and `dt/2` diagnostics. Only a conservative
   diagnostic pass permits the five unique `k_r/a_v` cases.

This concentration change is a controlled scenario, not a correction to the
x=0.05 result and not evidence that x=0.08 is experimentally representative.
