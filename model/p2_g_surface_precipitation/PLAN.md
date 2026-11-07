# P2-G registered benchmark

Base: 3a84f863eae63da8a09aef1ca2e78c5053daaaeb. Separate branch/module.
Seven one-dimensional cases, 1000 s maximum, 40/80/160 vertical cells;
dt=.5/.25 s; continuum switch .5/1/2 micrometres; no-precipitation control.
Gas film uses frozen Pitzer activity and diffusion correlation, held at the
previous accepted surface composition for each step (first-order splitting).
Only surface BOX vertices undergo a post-transport conservative projection.
Remove salt, preserve water and invert density/saturation to recover p,x.
Solid inventory is per surface area; h_mass=solid/[2165*(1-.43)].
Prat finite-Da diagnostic uses d_sub=25 um, k_r=.0023, eps=.43, beta=2,
D=sqrt(eps)*1.3e-9 and a_v=36*(1-eps)/d_sub. Its implicit internal reaction
is a comparison hypothesis, not an implemented physical process.
The diagnostic starts at the actual mass height crossing the chosen switch;
later predicted increments never debit soil or control evaporation.
No crust pore-water inventory or hydraulic resistance is modelled. Actual
height is a solid-equivalent geometric diagnostic, not a verified wet crust.
Gates: relative total salt/water budgets <1e-8; projection water unchanged;
accepted top w<=Ceq+1e-12; all nominal runs must activate precipitation and
cross continuum threshold. Compare grid/dt differences against 2%/1%; failure
limits prediction claims, never permits threshold tuning. Independent Prat
demand compatibility band is 10% (exploratory), not physical validation.
Retain failed runs and full outputs on HPCC; retrieve compact accepted records.
