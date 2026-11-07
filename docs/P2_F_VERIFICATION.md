# P2-F verification

## Scope and execution identity

All compilation and numerical testing was executed on HPCC Slurm compute
nodes. No test was run on the HPCC login node or on the local workstation.

- frozen base commit: `02c352f8c0f64f9dff952121b78ae7a61259a2e6`
- branch: `research/p2-f-analytical-growth`
- pinned source-stack archive SHA-256:
  `2a29d942b87c157eb0bd58dd35577c1088932a048411e069837bf36b80cad684`
- built solver SHA-256:
  `e4effd7b257fc661b553d1b1d3825f5579fbb68a0fa514785cd36afae9591711`
- compiler/build: GCC 12.2.0, CMake 3.25.2
- isolated verification job: `16367196`
- successful DuMuX build job: `16367199`
- final coupled array: `16367255`
- final acceptance analysis: `16367281` (exit 2 is the registered NO-GO result)

The local result bundle contains only CSV/JSON/text evidence. VTK fields and
binaries remain on HPCC at
`/WORK/hulm_work/zhixin/p2_f_analytical_growth_20261107`.

## A. Analytical growth law

Twelve Python tests and an independent C++17 implementation passed. The
cross-language values of equilibrium mass fraction, density, effective
diffusivity, growth velocity and implicit constant-flux reference height were
identical at printed double precision; maximum relative difference was zero.

The tests establish:

- zero evaporation gives exactly zero growth;
- growth is nonnegative for nonnegative prescribed evaporation;
- the joint `a_v -> 0` limit of the finite-reaction expression gives
  `Da_top/(1+Da_top)`;
- an independently solved reaction-free advective boundary-value problem
  approaches the closed law when `Pe_U` is small and differs when advection is
  deliberately large;
- the large-`Da_top` limit is recovered;
- porosity is immutable and neither `solveProfile` nor internal specific area
  appears in the production model.

The appropriate reaction-free validity measures are `Pe_U=|U|h/(eps D_eff)`
and `Pe_move=vh/D_eff`. The old condition `Pe << Da` is not used because it is
singular as `Da -> 0`.

## B. Isolated crust

For `h0=1e-5 m`, `j=1.8e-5 kg m-2 s-1` and 2000 s, the independent implicit
reference is `2.0379285095813205e-5 m`. Explicit Euler errors for time steps
40, 20, 10 and 5 s were respectively `9.9996e-10`, `4.9874e-10`,
`2.4906e-10` and `1.2445e-10 m`, giving the expected first-order decrease.
Solid salt, dissolved pore salt and pore water increments exactly match the
flux-defined inventory changes within the registered floating-point tests.

Pure trial evaluation is repeatable; rollback leaves persistent state
unchanged; discarded, stale and already committed trials are rejected.

## C. One-dimensional DuMuX coupling

The integration overlay retains the P1-R2 total-liquid/NaCl mass equations.
Before activation it applies prescribed water evaporation and zero salt flux.
At the first accepted solubility crossing, a declared seed is debited from
soil salt and water. After activation, Newton evaluations are pure trials and
height is committed once after acceptance.

The final benchmark used one horizontal cell, 40/80/160 graded vertical
layers, `x_NaCl,0=0.09999 < 0.1`, `j=1e-6 kg m-2 s-1`, `dt=0.025 s`, and a
closed bottom. The registered interface envelope was unchanged:
`|C_i/C_eq-1| <= 0.02`, `S_w >= 0.99`, `Pe_U <= 0.1`, `Pe_move <= 0.1`.

| Case | activation (s) | last accepted (s) | outcome | salt residual | water residual |
|---|---:|---:|---|---:|---:|
| 40 layers, 10 um | 52.675 | 60.000 | completed | 4.45e-14 | 3.93e-14 |
| 80 layers, 10 um | 12.350 | 24.075 | interface gate at 24.100 s | 1.08e-14 | 7.09e-15 |
| 160 layers, 10 um | 1.575 | 1.600 | interface gate at 1.625 s | 5.31e-16 | 1.72e-15 |
| 80 layers, 10 um, dt/2 | 12.3375 | 24.025 | interface gate at 24.0375 s | 2.09e-14 | 1.72e-14 |

Residuals are relative combined soil-plus-crust component residuals. Across
all nominal cases, activation-ledger relative residuals were at most
`1.92e-16`; combined salt and water residuals remained below `4.46e-14`.
Maximum `Pe_U` was `6.84e-5` and maximum `Pe_move` was `6.79e-6`, so the
rejections were not caused by the neglected-advection conditions.

The grid result is a failure of the required convergence gate, not evidence
for a converged activation prediction. Activation changed by factors far
larger than the registered 2% tolerance, and the finer grids left the
equilibrium interface envelope earlier. A shared endpoint comparison is
therefore undefined. The half-step 80-layer run reproduced nearly the same
activation and domain-exit times, showing that reducing `dt` did not remove
the supply incompatibility.

## D. Supply and seed sensitivity

At 80 layers, changing seed height altered the time at which the interface
left the validity envelope:

- 5 um: activation 12.350 s, domain exit 58.225 s;
- 10 um: activation 12.350 s, domain exit 24.100 s;
- 20 um: activation 12.350 s, domain exit 15.475 s.

Thus the seed inventory is conservative but the sustained equilibrium-supply
assumption is strongly seed dependent. A deliberately dry/high-flux case
activated at 0.825 s and was stopped at 0.850 s because interface saturation
was 0.715, far below 0.99. Growth was not clipped.

## Investigated failed attempts

Failures were retained rather than relabelled:

- build `16367197` stopped before compilation because the script referenced a
  nonexistent optional P1-R2 `source/dumux` overlay; build `16367199` removed
  that packaging error and passed;
- array `16367214` stopped before solver initialization because this OpenMPI
  build cannot be direct-launched through `srun`; direct single-process launch
  was used subsequently inside the Slurm allocation;
- exact `p_c=0` initialization produced a NaN Jacobian in the inherited
  saturation endpoint. Moving to `WaterLevel=0.01 m` kept `S_w>0.99` while
  avoiding that singular endpoint;
- high prescribed flux `1.08e-5 kg m-2 s-1` caused the interface-equilibrium
  gate to fail even more rapidly and is retained as insufficient-supply
  evidence, not tuned away.

## Assessment

The analytical model and conservative bookkeeping are verified. The present
Richards boundary coupling does not provide a grid-converged, hydraulically
supported equilibrium interface. This is an evidence--conclusion mismatch for
any claim that P2-F is ready for two-dimensional density-instability studies.
