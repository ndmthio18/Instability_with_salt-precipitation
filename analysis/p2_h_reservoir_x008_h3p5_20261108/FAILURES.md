# Failure records

## Numerical pilot

- Job: 16438585.
- Requested horizon: 180000 s.
- Last accepted time: 23613.815040957943 s.
- Accepted activation: 23116.0 s at 3.500127024 µm actual solid-equivalent height.
- Solver exit: 1, `ALE accepted state outside M3 physical bounds`.
- Scientific status: `FAIL_ALE_UPPER_POROSITY_BOUND`.

## Independent spatial reproduction

- Job: 16438590.
- Underlying solver failure reproduced on the same compute node with full ALE
  accepted/rejected fields.
- Last accepted `epsilon_max`: 0.9499999999125653, cell 0.
- First rejected `epsilon_max`: 0.9500526589230762, cell 0.
- Last accepted `epsilon_min`: 0.19093606251501768, cell 34.
- Concentration bound violations: zero.
- Rejected trial states: 25.

## Postprocessing failures retained

- Job 16438601: compute node Python lacked matplotlib; failed before analysis.
- Job 16438603: categorical `status` field was incorrectly parsed as float.
- Job 16438605: corrected postprocessing completed.
- Job 16438615: final donor-mixture provenance update completed.

No numerical simulation was repeated or altered for either postprocessing
failure.

## Deliberately not done

- No porosity or concentration clipping.
- No relaxation of the 0.95 upper bound.
- No artificial 3.5 µm seed.
- No forced equilibrium pore liquid without an explicit mass ledger.
- No ALE-grid, timestep, `k_r`, or `a_v` extension after the pilot failed.
