# P2-H staged implementation plan, 2026-11-07

User instruction: task/9.md. Base P2-G a3dae1374705f55a08d9d146c55567a3bb1becd0;
external salt_grow reference 8c0df6eddeea95a8a76c61f245cef2f891255c86.
Isolated branch research/p2-h-ale-coupling. Frozen branches/files remain intact.

A: Port the exact M3 finite-volume operator to C++ and use SUNDIALS 6.7.0
CVODE BDF/Newton, not explicit integration. Reproduce B2, 640 cells, power 3,
21600 s, max step 60 s, rtol 1e-8, field/inventory atol 1e-12, height atol 1e-14.
Target h=1.64919695 mm, eps_min=.23612005, internal fraction=.25378548,
V=3.00578141e-8 m/s. Cross-integrator scalar tolerance 5e-5 relative;
independent global salt/water residual <=1e-8 relative. Also compare complete
frozen output histories and sampled spatial operators, including reverse mesh
velocity. First run 80-cell/300-s smoke, then 640-cell/6-h regression.
Source-only SUNDIALS is built/installed within this isolated HPCC workspace.
Original release archive SHA256 d434d00a68d0f18dc2c8f4c306e4becf62ac6d37a259253ec9ea616152babba4.
Minimal source archive SHA256 d5ef5a9227d1ecb002c784def3975a6f099c6bd88fd0d6f802f124f43420e3ae.
The external M3 Python implementation is used only as a verification oracle.

B: Keep P2-G pre-crust transport/projection. Begin with zero external solid.
At 1/2/5/10 um solid-equivalent height, transfer existing solid once and
withdraw the required pore solution only from surface BOX control volumes,
at their actual accepted concentration. No distant donor or invented liquid.
Recover soil primary variables preserving removed water and dissolved salt.
Refuse activation if donor inventories or retained surface saturation fail.
Thresholds 1/2 um are below the estimated 3.5-um crystal size.

C (conditional on A/B): Trial ALE at the candidate soil-interface concentration,
apply its time-integrated total-liquid/salt flux to DuMuX as equal and opposite
Neumann exchange, and solve the interface concentration consistency equation.
Gas-film flux is computed from ALE exposed-top concentration. Remove the old
surface precipitation sink after activation. Commit both domains together;
discard both on rejection. No assumed crust permeability is introduced.
Demonstration restricted to continuous saturated crust supply with accepted
soil interface Sw>=.99; this is an explicit screening condition, not a derived
hydraulic potential or unsaturated crust closure. Detect failure rather than
silently clipping growth or continuing with a prescribed reservoir.

D/E: If coupled states are admissible, compare P2-G and Prat over actual common
accepted times. Initial bounded cases: soil 80/ALE 64/dt .5 s, soil 160,
ALE 128, dt .25 s, activation 1/2/5/10 um; 2500 s cap so the pre-crust control
can reach the 10-um threshold. Predeclare relative budgets 1e-8, source/transfer
balances 1e-12 (scaled by initial domain inventories), interface concentration
residual 1e-10 absolute. Report grid/time/threshold differences without tuning
away failures; 2% grid and 1% time gates for height/solid/evaporation.
These are short 1D demonstrations, not new 2D production or physical validation.

Implementation audit amendment (before B/C): the same-tolerance primitive
CVODE regression matched all four requested endpoint values, but its nonlinear
inventory reconstruction failed the stricter global budget (salt 4.264e-8,
water 1.425e-7 relative). Preserve job 16371399 as failed. The next implementation
integrates E and Q directly with the exact M3 flux operator. Tolerance coefficients
remain 1e-8/1e-12/1e-14, mapped to E/Q; the error norm differs from the primitive
reference. No residual is redistributed or corrected. Repeat all A gates before B.

Interface audit amendment: job 16371422 exposed a legacy final assignment
zeroing salt flux after the new boundary branch. Preserve the failed C histories
and binary; correct the overlay exclusion, explicitly log cumulative interface
liquid/salt exchanges and audit them independently against soil inventories.
Rebuild and repeat B/C. Add a 1e-6 s minimum retry-step stop to prevent numerical
stall prefixes from being extended to meaningless tiny timesteps. No mass-budget
threshold or physical parameter is relaxed.

If a required physical closure prevents C, deliver verified A/B and label
subsequent calculations diagnostic-only, as explicitly permitted by task/9.md.
Output four-panel vertical porosity maps and separate P2-G comparisons from
accepted records only. Preserve compiler/test logs and all failed attempts.
