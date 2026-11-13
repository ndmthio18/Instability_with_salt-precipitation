# Core solid fraction: definitions and identifiable evidence

The analytical prediction is phi_core=(1−epsilon_t)(1+X)/X, with X=sqrt(Da*) and all transport and internal-area coefficients frozen at the birth porosity. Its numerical value is recomputed from the registered material inputs in `results/constants.json`; it is not used by the ALE.

Three different quantities are reported:

- Birth solid fraction 1−epsilon_t=0.57, a boundary material parameter.
- Volume mean phi_bar=Mcr/(rho_cr hc), obtained from independently evolved solid inventory and actual geometric thickness.
- Core mean and spatial standard deviation over the overlap of finite-volume cells with 3 ell < z−zb < hc−3 ell, using the same explicitly fixed reference reaction length ell for all cases. A nonexistent core is NaN, not zero.

A saturated core is not automatically a spatially uniform, stationary-density core. The analytic derivation additionally requires stable translating reaction layers and complete conversion of the entering reference flux in the upper zone, negligible dissolved-storage corrections, and compatible top/bottom core densities. Material created at earlier times can carry its own solid-fraction history down to the bottom front.

The ALE evolves av=36 phi/d and eps-dependent diffusivity eps^(3/2)D. These differ substantially from coefficients frozen at eps_t. For example the transport coefficient at phi=.84 relative to that at birth is (.16/.43)^(3/2), while internal area is .84/.57 times the birth value. This algebra establishes a violation of the frozen-coefficient approximation; it does not isolate the causal share of each mechanism in the observed core error.

Finite-thickness H=hc/ell, profile resolution, κ-dependent recession, core spread/drift and volume-average evolution are separate diagnostics. A coefficient-freezing counterfactual would be needed to quantify a unique mechanism; an unexplained residual cannot be called caused by variable coefficients merely because those coefficients vary. If numerical gates remain unresolved, neither recovery nor rejection of the theoretical constant is definitive. No replacement constant is fitted to the simulations.
