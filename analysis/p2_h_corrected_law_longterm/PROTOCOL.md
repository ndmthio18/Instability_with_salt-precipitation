# Registered analysis and adaptive refinement record

P0 matrix: u=0,.3,.6,.9 and x008; n=32,64,128 at power=1.5;
64-cell power=1,2 sensitivity; rtol=1e-8/maxstep=30 s and
rtol=1e-9/maxstep=15 s paired temporal sensitivity. A paired change bounds
combined temporal sensitivity; it does not identify separate tolerance and
maximum-step contributions.

Target times are 1,5,10,20,50 h (x008 support ends at 44.70875 h).
Reports keep threshold-stopped intervals and do not extrapolate them to targets.
Output times are logarithmically spaced during start-up and 60 s thereafter;
profiles every 600 s and at the final structural threshold. Events are bracketed,
not claimed to be root-located.

Acceptance: reconstructed global and local inventory errors <2e-9 kg/m2;
top-speed spatial/temporal difference <1% on the common interval. Concentration,
reaction, solid and dissolved inventories are reported as additional convergence
observables. Passing total mass conservation alone cannot establish convergence.

Adaptive refinement was triggered by the measured 8.86% late top-speed difference
between n64 and n128 at power=1.5. Added n128,256,512 at power=2, with the same
physical parameters. This changes numerical resolution, not the analytical law.
The original failed resolution is retained. These controls remain mathematical
continuum calculations: h0=3.5 um and ell0=0.779 um are smaller than the 25 um
grain parameter. Subgrain mesh sizes do not independently validate the continuum
closure or mechanical integrity of the skeleton.

P1 changes beta and a_v separately. It is not authorized by a solver completion
alone: P0 budget, spatial resolution and structural validity must first be
assessed. P2 similarly remains an interface-consistency diagnostic unless a
hydraulic closure is available; no correction factor is introduced.

P1 registered 24 runs with n=128/256, u=0/.6. Eight low-Da* runs reached an
accepted-state bound between sparse output requests. Their valid prefixes and
failure logs are preserved; eight distinct `_event` reruns add 1 s output
requests to detect the epsilon=.05 screen before that bound. No physical
coefficient was altered. A measured >1% 128/256 top-speed difference then
triggered five n=512 checks: both beta=.5 and avmult=16 cases, plus jfactor=1000,
u=0. This finite adaptive check does not promise convergence of the local
porosity extreme or threshold time; unresolved cases remain labeled as such.

The quasi-steady BVP C freezes a measured material field and shares its mesh.
It is not independent structural prediction. Agreement with C cannot replace
the separate ALE mesh checks. High-Pe cases have a nominal clock-scaled output
schedule; this does not shorten the physical diffusion relaxation time by fiat.
