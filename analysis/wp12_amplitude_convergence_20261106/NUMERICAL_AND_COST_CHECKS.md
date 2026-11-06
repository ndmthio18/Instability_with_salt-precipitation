# Numerical, provenance and cost checks

## Scope

This record covers the six newly authorized 6T simulations at amplitudes
`1e-4` and `1e-5`.  The three existing `1e-6` trajectories were reused read-only.
Scientific convergence is assessed separately in `REPORT_ZH.md`; passing this
record does not validate a solubility-crossing time.

## Frozen identities

| Item | SHA-256 |
|---|---|
| DuMuX executable | `112076759fdff0942000307f729c2bd3c7e96963229b43c1b67c1f12bbf9df81` |
| `main.cc` | `dc7c629220fd0f3b23a6fcb6f156537d0d1d38b376550606235aeda7729dd2b4` |
| A28 eigenprofile | `d721e605d126e555f92b43200aa73b999ed4e2b61397edd1e98c0d9fddb36575` |

Input SHA-256 values:

| Amplitude | R24 | R48 | R96 |
|---|---|---|---|
| `1e-4` | `ebcab330653ad76495be4ac960e2e8c622402d0c60fb1d55e58e7d0d7101f0db` | `ea7993e1b65fc5329fe1da9d3725641bfcdf6d0a715c8b227e469fcb9f986ca2` | `582c35ad9a6d7b57b447f40d746231abe3cecc84c3c31a1db854d66aa6b29424` |
| `1e-5` | `c1c7ecee27ed49e10f70bed0aa78e12077fe2220acfb20cf5df00cb1b2842e88` | `9eb09a60c355eb81845985877cbde2b5843513651102b55189504d9f98d06d32` | `120a22fe9c53380cfd114396ffa7adafa64736fe85d7806b94555d74fe0e4680` |

The input-difference audit found only `Problem.Amplitude`, `Problem.Name` and
the authorized 6T `TimeLoop.TEnd` difference relative to the reference template.

## Per-run acceptance

| Amplitude | Grid | Simulation job | Reduction job | Acceptance job | Stop | Acceptance |
|---|---|---:|---:|---:|---|---|
| `1e-4` | R24 | 16309362 | 16309390 | 16309403 | 7200 steps / requested horizon | PASS |
| `1e-4` | R48 | 16309363 | 16309391 | 16309399 | 7200 steps / requested horizon | PASS |
| `1e-4` | R96 | 16309364 | 16309392 | 16309400 | 7200 steps / requested horizon | PASS |
| `1e-5` | R24 | 16309618 | 16309621 | 16309690 | 7200 steps / requested horizon | PASS |
| `1e-5` | R48 | 16309619 | 16309623 | 16309624 | 7200 steps / requested horizon | PASS |
| `1e-5` | R96 | 16309620 | 16309625 | 16309626 | 7200 steps / requested horizon | PASS |

Each accepted run satisfies:

- registered job, grid, four-thread binding, executable/input/profile hashes;
- solver exit zero and node-local scratch removal after durable archive creation;
- 7201 consecutive scalar records with strictly increasing accepted time;
- fixed nominal timestep, except the registered final horizon clamp of
  `8.4703501102e-8 s`;
- 145 complete VTU frames and one PVD mapping on the intended grid;
- independently reconstructed salt/water storage and conservative-flux budgets;
- measured initial concentration range equal to `x0 ± amplitude`.

## Conservation and online-versus-field agreement

Across the six new runs:

- maximum independent salt-inventory reconstruction difference:
  `2.1833121248426846e-17 m`;
- maximum independent water-inventory reconstruction difference:
  `5.665606872540252e-15 m`;
- maximum online-versus-field weighted diagnostic difference:
  `5.880712583561376e-16`;
- cumulative salt and water flux reintegration differences: exactly `0` in the
  recorded floating-point representation.

The detailed values for every observation volume are stored in each
`archive_acceptance.json` and `wp12_scalar_metrics.json`.

## Runtime, allocation and storage

| Amplitude | Grid | Solver wall time (s) | Allocated core-hours (4 cores) | Peak RSS (KiB) | Field archive bytes | Product archive bytes |
|---|---|---:|---:|---:|---:|---:|
| `1e-4` | R24 | 901.36 | 1.0015 | 31,176 | 33,493,798 | 12,171,976 |
| `1e-4` | R48 | 2,343.32 | 2.6037 | 48,212 | 64,195,704 | 12,283,322 |
| `1e-4` | R96 | 6,309.75 | 7.0108 | 82,424 | 125,507,841 | 12,403,412 |
| `1e-5` | R24 | 901.37 | 1.0015 | 31,252 | 32,764,856 | 12,053,580 |
| `1e-5` | R48 | 2,283.46 | 2.5372 | 48,148 | 62,327,648 | 12,147,383 |
| `1e-5` | R96 | 5,708.98 | 6.3433 | 85,344 | 121,933,419 | 12,288,958 |

Allocated core-hours are `4*solver_wall_time/3600`; they are a reproducible
allocation cost, not measured CPU utilization.  R96 live scheduler sampling
showed about 1.6 actively used cores and about 86 MiB RSS, consistent with the
known partially serial solver structure.

Before launch, `/WORK` had about 231 GiB headroom and the `cnmix` partition had
ample idle CPUs.  The account was above its soft file-count quota, so all durable
output inodes were pre-created, VTUs were written on node-local scratch, and each
run was collapsed into two durable archives plus fixed evidence files.  No
unbounded file creation or login-node computation was used.

## Preserved failure evidence

- Initial `1e-4` jobs 16309354--16309356 failed before the solver because fixed
  evidence files were staged in the wrong directory.  The logs were preserved;
  the same frozen inputs were submitted as 16309362--16309364.
- The first R96 reference extraction (16309369) failed only on an over-strict tar
  member-order assertion; jobs 16309376/16309386 produced the verified 0--6T
  prefix.  No production trajectory was modified.
- R24 `1e-4` acceptance 16309398 initially rejected the legitimate terminal
  horizon clamp.  The gate was corrected and 16309403 passed without rerunning
  the simulation.
- R24 `1e-5` acceptance 16309622 initially used the hard-coded `1e-4` initial
  range.  The gate now derives `x0 ± amplitude` from the archived input;
  replacement acceptance 16309690 passed without rerunning the simulation.

These failures are workflow/gate evidence, not hidden solver retries or favorable
trajectory selection.
