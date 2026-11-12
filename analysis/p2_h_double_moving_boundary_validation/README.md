# Task 27: independently computed geometric salt-crust growth


**Final scientific status: NO-GO for final quantitative net-growth validation.**
Physical/conservation gates passed. Twenty-four of 26 registered trajectories
reached their supported endpoints (50 h controlled forcing; 44.70875 h x008).
The u=0.6 rtol=1e-9 and 1e-10 temporal controls stopped after 2 and 5 s of saved
output. All baseline/spatial-refinement runs completed, but this unresolved
numerical robustness issue prevents claiming full temporal convergence.

Start with [the Chinese validation report](VALIDATION_REPORT.md). Eight figures
are in figures/ (PNG, PDF, Times New Roman text SVG). Full computed histories
and profiles are in results/*.csv.gz; failed cases are explicitly labeled.
The results separate actual geometry from crystal-inventory-equivalent height.

Baseline c7da956041b7e42f85278603a068c3915394ec20; branch
research/p2-h-double-moving-boundary-validation. Original ALE/BDF and all
previous studies remain unchanged. New code is isolated in this directory.

Read GOVERNING_EQUATIONS.md and BOTTOM_FRONT_CLOSURE.md first. The selected
bottom law assumes exposed reactive area kappa_b phi_b and uses local
undersaturation with the existing intrinsic k_r. It is an explicitly added
physical closure, not a prediction uniquely determined by conservation.
No theoretical core fraction or corrected growth formula enters the solver.

Source: two_front.hh, two_bdf.hh, run_case.cc, verify.cc. The first advances
both regions and both interface positions; the second is an independent CVODE
wrapper. The original ale.hh is a read-only source for the top kinetic root;
the original bdf.hh is used only for recovery verification.

Production uses HPCC Slurm cnmix at
/WORK/hulm_work/zhixin/p2_h_double_moving_boundary_20261112.
Compiler GCC 12.2.0 and SUNDIALS 6.7.0 are pinned as in task 26.
build_verify.sbatch runs the physical verification gate before building the
production executable. run_pilot.sbatch and check_pilot.sbatch verify one-hour
transport and seed sensitivity before run_production.sbatch is eligible.

Reproduce the registered matrices with `python prepare.py`. After retrieving
raw histories, profiles and status files into raw/, run `python analyze.py`.
The x008 forcing file is copied verbatim from the preceding comparison study.
No extrapolation beyond its 44.70875 h support is allowed.

`python -m unittest -v test_remap.py` tests conservative overlap transfer and
the exact-zero connected-prefix rule. Full results-level tests are separate
from physical validity. Intermediate failures are retained by distinct case
IDs/source snapshots rather than overwritten by retries.

Geometric thickness is always zt-zb. Inventory-equivalent thickness and the
integrated analytical prediction are labeled separately. A contact-flux version
of the top law is a diagnostic conditional on the computed boundary flux,
not a standalone analytical prediction of that flux.

After raw data retrieval, run `python analyze.py`, `python audit_fluxes.py`, `python write_report.py`,
and `python audit_delivery.py`. The delivery audit checks evidence integrity;
it does not override the scientific NO-GO. Download archives from the HPCC
work directory named above: task27_evidence.tar.gz, task27_refine.tar.gz and
task27_temporal.tar.gz; their SHA256 receipts are retained in provenance/.
Source snapshots match all 59 archived source-hash entries. The source in
provenance/source/ is the exact compiled production source; root source files
must match it byte-for-byte. Old original ALE/BDF operators are unchanged.

Early contact-diffusion spikes are also unconverged; see results/flux_resolution.csv. Figure 8 shows startup data separately, without deleting the spikes or interpreting them as established physical events.
