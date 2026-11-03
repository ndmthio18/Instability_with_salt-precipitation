# Isolated WP6 seven-plane transport diagnostics

Copied from `model/transport_v1/source`; the only source modification is
`appl/simpleRichardsnc/transportdiagnostics.hh`. The physical residual,
constitutive laws, boundary conditions, initial condition and original detector
in `main.cc`/`model.hh` are byte-identical to the validated version.

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

Build/runtime/source hashes, threshold fixtures, bridge and 2T results belong in
`analysis/wp6_execution_20261103`. Raw fields remain in ignored local cache and
the dedicated remote execution directory. No old source or run is overwritten.
