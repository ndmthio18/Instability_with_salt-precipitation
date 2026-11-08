# Reproduction and evidence layout

Baseline: `25a41470aac2fa52ae772d313f9ae2b35a49339f` on
`research/p2-h-ale-coupling`; development branch `research/p2-h-adaptive-donor`.
Read the registered [matrix and gates](../../model/p2_h_ale_coupling/ADAPTIVE_PLAN.md)
and [initialization equations](../../model/p2_h_ale_coupling/ADAPTIVE_FORMULATION.md).

HPCC root: `/WORK/hulm_work/zhixin/p2_h_adaptive_donor_20261108`.
All builds, numerical audits and plotting execute on allocated compute nodes.

1. Export the module and frozen dependencies to `repo/model/` as listed in
   `source_r01.tar.gz`: `p2_h_ale_coupling`,
   `p1_r2_conservative_richards/source`, `p2_g_surface_precipitation`,
   `p2_d_coupling/cpp`, `p2_f_analytical_growth/integration`.
2. Provide the exact-stack archive identified by `build_adaptive.slurm` and
   the existing task-local SUNDIALS 6.7.0 static install. This run reads the
   unchanged P2-H install through a symlink; it never rebuilds inside the old
   workspace. The standalone M3 operator and BDF are byte-identical to the
   base, whose numerical regression remains documented in the prior report.
3. Submit `repo/model/p2_h_ale_coupling/build_adaptive.slurm`. It builds the
   actual application, executes C++ donor unit tests, 51 controlled allocation
   cases and two 400-s DuMuX smoke tests. Logs and binary SHA are retained.
4. Submit `matrix_adaptive.slurm` only after successful smoke execution.
   It audits smoke/controlled histories before starting 63 registered runs
   with eight independent processes, one thread each. The exact executable
   argument array and input file are saved per case. Existing case directories
   are never overwritten by the runner. Individual timeout is 1800 wall seconds;
   a timeout or solver failure is an outcome, never a successful coupled run.
5. After the matrix ends, submit `postprocess_adaptive.slurm` using the final
   tracked analysis scripts. It checks complete accepted histories, independent
   donor mass reconstruction, interface budgets, old-mode regression and
   exploratory convergence gates, then produces figures and archives.
   Independently run `refusal_adaptive.slurm` before this audit: its two real
   400-s histories must be exactly identical, proving repeated refused adaptive
   attempts do not alter the precrust transport state. The final audit also
   compares frozen source hashes, including the complete original soil application source.
6. Verify `transfer.sha256` on retrieval. The compact archive contains
   sampled histories, full allocation profiles, case commands and inputs,
   verification tables/JSON, figures, compiler logs, job logs and source/output
   checksums. Full accepted histories remain in a separately hashed archive
   outside Git; raw VTK and ELF stay on HPCC. Do not integrate sampled CSV rows
   as if they were complete time-step records.

The optional `resume_adaptive.slurm` is a scheduling safeguard, not a solver
restart scheme. Under a primary Slurm TIMEOUT it moves interrupted prefixes
to `failures/16371535/` and replays only missing/interrupted cases from the same
frozen initial conditions. It never overwrites prefixes or automatically retries
numerical failures. `resume_metadata.json` records whether this path was used.

Elapsed Slurm time is wall time, and allocated CPUs are not measured CPU use.
An exit code zero can mean a completed **precrust** case or a controlled
refusal. Scientific status requires the separate activation and convergence
fields in `audit/verification.json`. Source archive r01 predates added analysis
scripts; the final remote source manifest records the files actually audited.

The read-only attribution is reproduced with `boundary_audit.slurm` after the
matrix: it uses actual BOX coordinates, the frozen boundary predicate and
the exact upstream BOX SCVF geometry headers. It does not repair histories.
The visualization-only `plot_revision.slurm` was executed after rendered QA
to improve conservation tick labels, align donor-mass axes and remove crowded
initial-saturation labels. `figure_revision.tar.gz` supersedes the compact
archive's figures; numerical tables and histories are unchanged. The initial
remote manifest predates this plot-only change; `provenance.json` and
`local_delivery.sha256` record the final delivered revision.

Actual jobs, failures and retrieval scope are in [FAILURES.md](FAILURES.md)
and `final_accounting.txt`. The raw/full accepted-history archive is retained
on HPCC; only the compact evidence and final figure revision were retrieved
and verified locally. Compiler logs are compressed under
`deliverables/compiler_logs/`, including the executed C++ unit test log.

## Diagnostic controls

- `P2H.DonorMode`: `top_only` (compatibility default) or `adaptive_multilayer`.
- `P2H.MaxDonorDepth`: physical metres, adaptive default `.005`.
- `P2H.MinimumSupplySaturation`: unchanged `.99`.
- `P2H.ControlledDonorTest`: prescribed-mass allocation harness only; exits
  before time integration and does not create a production salt seed.
- `P2H.PrecrustOnly`: continues actual surface precipitation without activation.
- Existing injected-rejection controls test activation and post-activation
  rollback independently. Attempt rows are not necessarily accepted states.
