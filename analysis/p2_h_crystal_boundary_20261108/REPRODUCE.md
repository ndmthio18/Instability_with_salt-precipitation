# Reproduction and evidence

Base `c1fb2e56ccee892206486283e5262cd41354a814`; branch
`research/p2-h-crystal-boundary`. The original branches are unchanged.

HPCC root: `/WORK/hulm_work/zhixin/p2_h_crystal_boundary_20261108`.
Use the `hpcc` SSH alias. All compilation, simulations, numerical auditing,
plotting and archive assembly execute in Slurm allocations, not on login nodes.

1. Read `model/p2_h_ale_coupling/CRYSTAL_PLAN.md`. Stage the module plus frozen
   dependencies exactly as recorded in `source_r02.tar.gz`: P1-R2 soil source,
   P2-G surface code, P2-D C++ Prat code and P2-F integration input.
2. Provide the exact-stack archive identified in `build_crystal.slurm`
   (SHA-256 `2a29d942b87c157eb0bd58dd35577c1088932a048411e069837bf36b80cad684`)
   and the existing SUNDIALS6.7 static install, read-only through `dependencies`.
   GCC12.2/OpenMPI4.1.4/CMake3.25.2 are loaded by the build script. Create `build/`.
3. Submit `build_crystal.slurm`. It extracts a fresh stack, applies the isolated
   P2-H overlay, compiles the application and donor unit test, runs four controlled
   default-threshold allocations and the real 320-cell120-s smoke test.
   Compiled `main.cc`, ELF and full compiler log remain in `build/` on HPCC.
4. Submit `matrix_crystal.slurm` with `--dependency=afterok:<build-job>`.
   It verifies the binary and smoke/geometry/conservation gates first, then runs
   twelve cases. Each case records the exact CLI array and input. Directories
   must be new; the runner refuses to overwrite them. Each child is bounded
   to3600 wall seconds; no numerical failure is automatically retried.
5. After all cases terminate, submit `postprocess_crystal.slurm`. It audits full
   accepted histories, compares frozen source hashes and old-mode histories,
   checks actual SCVFs and mass-weighted birth fields, then plots and packages.
   A separate optional `thread_check.slurm` reproduces the original120-s smoke
   with `DUMUX_NUM_THREADS=1`; this run's exact-history comparison is retained.
6. `package_crystal.slurm` only rebuilds deliverables after a packaging change.
   It does not run the solver or alter histories. The source manifest hashes
   source archives and raw outputs. Derived archive hashes are exclusively in
   `transfer.sha256`, avoiding stale self-referential archive entries.

For another workspace, update hardcoded roots and historical reference paths
in a separately recorded operational revision. The verification uses prior
frozen P2-H source and prior adaptive-donor old-mode histories as references.
The source archive predates later analysis additions; `remote_manifest.json`
identifies the final source files actually staged, and `provenance.json`
records the compiled solver and archive identities.

Actual jobs:

| Job | Role | Terminal outcome |
|---|---|---|
| 16374128 | Build, C++ units, controlled and smoke | COMPLETED |
| 16374135 | Twelve-case matrix | COMPLETED; two child solvers exited1 |
| 16374139 | Initial audit/plots | COMPLETED |
| 16378274 | One-DuMux-thread120-s check | COMPLETED, identical history |
| 16379345 | Final audit with component budgets and plots | COMPLETED |
| 16379358 | Packaging correction | COMPLETED |

See `final_accounting.txt` for node names, exit codes, allocated CPUs and wall
times. Eight allocated CPUs do not imply eight truly single-thread solvers:
the production environment constrained OMP/BLAS but retained TBB auto workers.
This execution deviation and the independent thread check are in the plan and
report; the production cases were not restarted or silently retuned.

`deliverables/audit/verification.json` separates implementation integrity from
convergence and physical validation. The two320-cell coupled failures are
accepted-prefix outcomes at1417.629079 and1377.976931 s. Their final2500-s
heights and fields do not exist. Do not treat scheduler exit0, successful
activation, or prefix mass closure as successful full coupling.

The compact archive includes sampled histories, full donor and boundary
profiles, available soil/crust fields, commands, inputs, compiler logs and
audit/figure products. Full accepted histories and raw failure records are
in the separately hashed full-history archive on HPCC; ELF, VTK and full
solver stdout remain there too. Do not integrate sampled history rows as if
they were all accepted steps. No failed state or missing profile was invented.

Initial source archive r01 was rejected by Slurm for CRLF before any job ran;
r02 fixes line endings. Both source hashes remain recorded. The packaging
correction removes prior derived archives from the source manifest; the final
`transfer.sha256` is the authoritative archive checksum record.
