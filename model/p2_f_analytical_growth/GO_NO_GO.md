# P2-F GO/NO-GO

## Decision: NO-GO for two-dimensional coupling

GO is supported only for the isolated reaction-free analytical law, its C++
implementation, transactional accepted-step update, seed conservation and
combined component bookkeeping.

NO-GO applies to production Richards--crust coupling and any 2D instability
calculation because:

1. activation time is strongly grid dependent (52.675, 12.350 and 1.575 s on
   40/80/160 vertical layers);
2. the 80/160-layer runs leave the declared equilibrium interface envelope;
3. seed height changes the valid coupled duration by more than a factor of
   four;
4. P2-F has no crust hydraulic-pressure closure, so adequate saturated supply
   is monitored but not solved.

The failures are not conservation or advection errors: combined salt/water
residuals are below `4.46e-14`, activation residuals below `1.92e-16`, and the
maximum Peclet measures remain below `6.84e-5`.

The next stage requires a grid-converged soil-side surface boundary layer and
a defensible hydraulic supply closure. It must not proceed by relaxing the 2%
interface gate or clipping the growth velocity.
