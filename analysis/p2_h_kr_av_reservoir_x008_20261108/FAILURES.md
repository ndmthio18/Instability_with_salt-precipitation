# Failure record

## Registered pilot failure

- Requested horizon: 180000 s (50 h).
- Last accepted time: 20580.510687092206 s.
- Solver exit: 1 (`ALE accepted state outside M3 physical bounds`).
- Scientific status: `FAIL_ALE_UPPER_POROSITY_BOUND`.
- Matrix status: `NOT_RUN_PILOT_GATE_FAILED`.

## Reproduced bound violation

An independent rerun on another HPCC compute node reproduced the same accepted
prefix and failure. The last accepted ALE field had no bound violations. The
first rejected trial had exactly one violation:

- cell: 0 (nearest the soil interface);
- height above soil: 2.00271947e-8 m;
- porosity: 0.950089672149217;
- upper bound: 0.95;
- concentration violations: 0;
- rejected trials recorded: 23.

The accepted minimum porosity was 0.3875437744268013, so the stop was not an
`epsilon_min = 0.02` clogging event. Net internal solid was negative because
cumulative internal dissolution exceeded positive internal precipitation.

## Actions deliberately not taken

- No porosity clipping.
- No relaxation of the upper bound.
- No salt or water inventory correction.
- No artificial seed.
- No matrix, ALE-grid, or timestep extension after the registered pilot gate
  failed.

