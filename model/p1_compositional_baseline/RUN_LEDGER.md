# P1 run ledger

All paths below are under `/temp-hd/zhixin/p1_compositional_20261106` on `sal`. No shared DUNE/DuMuX source was edited.

| Record | Outcome | Purpose |
|---|---|---|
| `build_attempt_01` | failed at configure | missing `p1compositional.pc.in`; retained |
| `build_attempt_02` | failed at compile | exposed SuperLU 7/DUNE 2.9 header incompatibility and local-residual constructor error |
| `build_attempt_03` | compiled | first real smoke run; default MPI/UCX initialization deadlocked |
| `smoke_v1_attempt02` | completed | TCP/`ob1` runtime workaround; V1 smoke residual at roundoff scale |
| `verification_run_02` | mixed/failed gates | diagnosed inactive-phase audit and ill-posed V5 setup |
| `verification_run_03` | all executions completed | corrected inactive-phase handling and two-phase hydrostatic gas pressure |
| `verification_run_04` | V3/V4/V6 gate failures | retained unchanged predeclared thresholds |
| `verification_run_05` | V3 failure only | V4/V6 corrections passed; constant-density gas hydrostatics still above velocity gate |
| `build_attempt_06` | compiled | second-order gas hydrostatic initialization, executable hashes recorded |
| `verification_run_06` | V1--V7 PASS | accepted verification evidence |
| `benchmark_cost_probe_01` | completed | B0, 280 × 100, 10 steps; 119.15 s wall, 67,996 KiB maximum RSS |
| `benchmark_run_01` | running/then final state documented in P1 reports | bounded B0/B1 Case 36877 comparison to 60,000 s |

Failed or superseded records are provenance, not accepted results. Acceptance criteria were not weakened after inspection.
