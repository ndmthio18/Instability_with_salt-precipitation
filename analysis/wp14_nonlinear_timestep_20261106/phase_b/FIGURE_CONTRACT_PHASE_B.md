# WP14 Phase B figure contract

## Figure 5 — R48 dt versus dt/2 surface concentration and modal evolution

- Population: Case 36877, R48 x 200, `A=1e-4`, original A28 eigenprofile, 0--6T.
- Comparison: frozen WP12 nominal `dt=113.68588647695695 s` versus the single WP14 refinement `dt/2=56.842943238478475 s`.
- Alignment: identical physical times only; original accepted step `k` is paired with refined accepted step `2k`.
- Panels: surface mean concentration, normalized nonseed/seed harmonic ratio, and `m=2,4,6,10` energy fractions for each timestep.
- Scope: a finite-timestep branch-stability test, not a proof of time convergence or physical validation.

## Figure 6 — R48 dt versus dt/2 conservative transport

- Population and pairing are identical to Figure 5.
- Observation planes: approximately 1, 2, 5 and 10 mm below the evaporating surface.
- Top row: conservative SCVF net advective salt flux, downward positive.
- Bottom row: time integral of the resolved positive part of the same advective flux over the frozen 0--6T window.
- Values below the registered `1e-15 m/s` resolvability threshold remain in the signed time series but are not reclassified as a transport event.
- These are transport-pulse diagnostics, not a validated physical `t_eff`.

## Figure 7 — Spatial versus temporal finite-resolution differences

- Horizontal difference: `Delta_x = |Q_R96,dt - Q_R48,dt|` from frozen WP12 runs.
- Timestep difference: `Delta_t = |Q_R48,dt/2 - Q_R48,dt|` using the single WP14 run.
- The scatter reports each difference relative to the corresponding finer-resolution value so heterogeneous QoIs can share one plot.
- The bar panel reports the absolute ratio `Delta_t/Delta_x` for preregistered representative QoIs.
- The shaded range `0.25 <= Delta_t/Delta_x <= 4` means that the two observed finite-resolution differences are comparable by the registered interpretation, not that errors are additive or statistically independent.
- No Richardson order or rigorous truncation-error bound is inferred, especially if nonlinear modal branches differ.
