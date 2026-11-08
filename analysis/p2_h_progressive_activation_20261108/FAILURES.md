# Failure register

## Pilot

The progressive 80-soil/64-ALE, `dt_max=0.5 s` pilot accepted conservative
activation at approximately 23.136 ks and then evolved for approximately
501 s.  It subsequently reached the unchanged M3 ALE physical-bound screen;
successive time-step reductions did not produce another admissible state.

Classification: `INITIALIZATION_IMPROVED_BUT_INTERFACE_FAILURE_REMAINS`.

This is not an initialization mass failure.  The initial ALE pore liquid came
from the progressive top-boundary ledger at the actual saturated interface
composition.  Positive internal precipitation and larger cumulative internal
dissolution subsequently coexisted, consistent with the previously audited
post-activation interface closure problem.  Final exact values and all matrix
failures will be populated from the retained accepted and rejected histories.

No failed run is extrapolated to the requested horizon, and scheduler
completion is not classified as scientific success.
