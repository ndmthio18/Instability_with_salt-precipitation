# Verification and numerical provenance

## Accepted physical gate

HPCC compute-node job **16572081**, GCC 12.2.0 / SUNDIALS 6.7.0.
The original `ale.hh` and `bdf.hh` are unchanged. The new solver integrates
crystal volume, dissolved-salt mass, two positions, and independent cumulative
flux/reaction ledgers. Gate values below are from actual job output.

| Test | Result | Acceptance criterion |
|---|---:|---:|
| Fixed-bottom recovery, maximum relative top-speed difference | 2.39773e-7 | <2e-6 |
| Fixed-bottom displacement | <2.8e-33 m; speed exactly zero | <1e-14 m |
| Recovery maximum mass residual | 2.79e-18 kg/m2 | <2e-9 |
| Pure dissolution recession vs independent constant-bath limit | 5.82766e-5 relative | <1e-4 |
| Pure dissolution maximum mass residual | 5.96e-19 kg/m2 | <2e-9 |
| Simultaneous motion, kappa=1: top / bottom speed at 100 s | 3.39833e-9 / 1.54850e-9 m/s | both positive |
| Simultaneous motion, kappa=2: top / bottom speed at 100 s | 3.40496e-9 / 2.11733e-9 m/s | bottom response >5% |
| Two-region maximum mass residual | 2.39e-18 kg/m2 | <2e-9 |
| Connected exact-zero transfer | passed | NaCl, water, crystal differences <1e-16 |

The pure dissolution problem disables top and volume reactions, uses a bath
at C=0.26 and D=1e-6 m2/s, and compares recession at 0.1 s against
rho/rho_cr * k_r * (Ceq-0.26) * t. Its finite diffusion resistance explains
the small nonzero asymptotic error. Crystal mass loss is checked independently
as rho_cr*(1-eps_birth)*recession. This benchmark diffusion coefficient is not
the production coefficient 1.3e-9 m2/s.

Recovery samples are 0.01, 0.1, 1, 10 and 100 s, with both solvers at relative
tolerance 1e-10. At 1e-9, one experimental Jacobian build missed the 2e-6
recovery threshold at the earliest sample; the threshold was retained and both
integrations were tightened. A later exact-equality test for bottom position
was replaced by the stated positional numerical tolerance after observing
2.3e-33 m linear-solver roundoff; its physical derivative remains exactly zero.

## Conservation interpretation

Each accepted trajectory compares reconstructed inventories with BDF-integrated
flux ledgers, including separate crust/liquid salt and water budgets. Interface
algebraic residuals are recorded separately: they do not constitute an
independent temporal conservation test. The geometric-rate identity is an
implementation check, not evidence for the corrected analytical growth law.

The zero-solid test transfers two connected bottom cells, retains a separate
interior zero cell, leaves the top fixed and preserves dissolved salt and water.
It does not claim that the original area-proportional-to-solid reaction law
reaches exact zero at finite time. Python overlap-remap tests independently
check mass preservation and decreasing reconstruction error at 32/64/128 cells.
No threshold deletion is allowed, even for a positive solid volume of 1e-30.

## Preserved development failures

| Job(s) | Formulation / result |
|---|---|
| 16572033, 16572036 | Original full-Q physical verification passed |
| 16572037, gate 16572044 | Five 1 h pilots passed; 0.1 nm seed failed at internal t=13.6836 s |
| 16572052 | Affine excess inventory, T-relative weights: pure dissolution failed at 0.0220593 s |
| 16572056 | Dependent pilot cancelled after its prerequisite failed; no computed result |
| 16572059 | Excess inventory with Q-scale weights: pure dissolution failed at 0.0397035 s |
| 16572064 | Bath-referenced excess inventory: extreme diffusion failed at 0.0619107 s |
| 16572069 | Less stiff pure-dissolution benchmark passed; kappa=2 failed at 24.4999 s |
| 16572074, 16572077 | Physical-scale Jacobian, excess inventory: early recovery threshold missed |
| 16572078 | Tight recovery passed; excess-inventory kappa=2 failed at 3.80311 s |
| 16572080 | Full Q plus physical-scale Jacobian: exact-zero assertion detected 2.3e-33 m |
| 16572081 | Full Q plus physical-scale Jacobian: all physical gates passed |

Failures are not renamed as successful registered production cases. Source
snapshots, hashes and stderr are retained. The evidence shows improved
robustness for the delivered representation; it does not isolate a unique
floating-point mechanism for every failed run.

The one-hour Jacobian pilots, production grid/time sensitivity, and final
trajectory bounds are summarized in VALIDATION_REPORT.md once available.

## Production and temporal controls

Six Jacobian pilots (16572082, gate16572088) all reached3600 s, including the
previously failing0.1 nm seed. The pilot maximum mass residual was4.10e-16 kg/m2.
At1 h, the0.1/10 nm seed net-speed differences from1 nm were about7.6e-9/1.7e-7
relative;64->128 crust-cell net-speed difference was0.58%.

The23-case main array16572089 produced22 full supported trajectories and one
failed temporal control: u0.6_strict at last output2 s (internal failure2.5201 s).
Both512-cell additions16572116 reached50 h. Additional strict-tighter job
16572125_2, using the SAME executable and rtol=1e-10, stopped at last output5 s
(internal failure7.26963 s). Thus24 of26 registered runs reached their endpoints;
two failed prefixes remain. Full-duration u0 and x008 strict controls passed.
The short u0.6 prefixes cannot establish50 h time convergence. Spatial error
reduction and mass conservation alone do not close this verification gap.
