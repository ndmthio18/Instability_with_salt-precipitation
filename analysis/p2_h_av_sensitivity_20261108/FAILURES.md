# Failure and exclusion record

No failure listed here is silently discarded or relabelled as a successful 2500 s physical run.

## Development and audit failures

| Event | Evidence | Classification | Resolution/use |
|---|---|---|---|
| Initial build job 16405764 stopped after the prescribed-state `av=1` mechanism case reached the physical porosity bound | `logs/build-16405764.out`, raw HPCC log | driver `set -e` classification error; frozen M3 regression had already passed | mechanism failure retained; smoke resumed with outcome-aware driver |
| Exact build job 16405925 failed `ALE state size` | `logs/build-exact-16405925.out` | unit-test fixture omitted the optional two ledger states | fixture corrected; exact binary rebuilt and unit test passed in 16405951 |
| Preliminary accepted-endpoint trapezoid integration did not close positive minus dissolution against the conservative net ledger for A3 | preliminary `primary`, `dt_half`, `ale_refine` retained in full archive | invalid diagnostic integration method | excluded from final tables; replaced by positive/dissolution BDF state ledgers |
| Audit job 16406099 failed converting CSV `reference` text to float | `logs/audit-16406099.out` | post-processing parser error | parser corrected; audits 16406144/16406274/16406416 completed |

For the preliminary A3 80-layer run, endpoint quadrature gave positive `0.241373 g m-2` and
dissolution `0.024790 g m-2` while net was `0.241750 g m-2`; this contradiction exposed the method.
The exact BDF ledger gives positive `0.242402059 g m-2`, dissolution `0.000652363 g m-2`, net
`0.241749697 g m-2`, with closure `-2.44e-19 kg m-2`.

## Registered physical/numerical stops

| Case | Last accepted (s) | First rejected (s) | accepted epsilon_min | rejected epsilon_min | Outcome |
|---|---:|---:|---:|---:|---|
| n80_av10, 64 ALE | 741.693716 | 741.693718 | 0.0200000026 | 0.0192916941 | STOPPED |
| n160_av10, 64 ALE | 741.693562 | 741.693563 | 0.0200000026 | 0.0192912656 | STOPPED |
| n80_av10_dt2 | 741.693787 | 741.693788 | 0.0200000006 | 0.0196742095 | STOPPED |
| n160_av10_dt2 | 741.444159 | 741.444160 | 0.0200000018 | 0.0196751395 | STOPPED |
| n80_av10, 128 ALE | 734.536174 | 734.536175 | 0.0200000020 | 0.0182056395 | STOPPED |
| n80_av10, 256 ALE | 730.647424 | 730.647425 | 0.0200000002 | 0.0182993063 | STOPPED |

The first rejected fields, local concentration, `a_v`, reaction, effective diffusivity and bottom interface
fluxes are retained in `results/reaction_profiles/*first_rejected.csv` and the compact case records.
These are physical-bound rejections of converged trial states, not generic Newton failures. Their grid-dependent
timing prevents a quantitative pore-clogging claim.

## Standalone mechanism outcomes

The independent prescribed-state test intentionally uses an aggressive matched boundary condition. `0.1av`
completed 5 s; `av` and `10av` reached the porosity bound after 4.4 s and 0.8 s. These results demonstrate the
direct kinetic multiplier and its feedback but are not DuMuX--ALE coupled predictions.

