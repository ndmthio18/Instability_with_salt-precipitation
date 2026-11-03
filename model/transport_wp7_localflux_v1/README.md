# Isolated WP7 local-SCVF transport diagnostics

Copied byte-for-byte from `model/transport_wp6_v1` before the WP7 edit. The
physical residual calculation, constitutive laws, boundary conditions, initial
condition and original detector remain unchanged. `model.hh` only exposes the
already-selected upwind mole fraction alongside the existing conservative flux
components; it does not alter any component returned to the residual.

The diagnostic requires seven absolute `Transport.CutCoordinates` values and
rejects any coordinate not matching a BOX conservative face to 1e-12 m. Original
deep region names stay unchanged. Four shallow regions are added; `free_domain`
remains last for the prescribed-bottom reaction accounting. Budget schema is 2.

`weighted_online.csv` adds accepted-step surface weighted standard deviation,
physical-wavenumber Fourier amplitude, volume-weighted eigenprojection and RMS,
and second/third top harmonics. Background eigenprojection is blank/undefined;
the original unweighted detector is untouched. Independent per-frame checks use
VTU tensor quadrature and reintegrate recorded face fluxes. They are independent
storage checks, not a separately implemented physical flux solver.

When `Transport.WriteLocalFluxes = true`, `local_vertical_flux.csv` records the
signed water, advective-salt, diffusive-salt and total-salt contribution of each
SCVF crossing one of the locked planes, together with face area, horizontal
coordinate and the residual kernel's upwind mole fraction. For every accepted
step and plane, summing an integrated component and dividing by domain width
must reproduce the corresponding unchanged `transport_budget.csv` value. This
identity is the bridge gate before interpreting down-flow fractions,
concentration-flow covariance or penetration. VTK velocity is not used.

Build/runtime/source hashes, threshold fixtures, bridge and 2T results belong in
`analysis/wp6_execution_20261103`. Raw fields remain in ignored local cache and
the dedicated remote execution directory. No old source or run is overwritten.
