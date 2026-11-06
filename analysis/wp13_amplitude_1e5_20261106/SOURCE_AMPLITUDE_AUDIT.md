# `Problem.Amplitude` source audit

- Source: `model/transport_v1/source/appl/simpleRichardsnc/main.cc`
- Source SHA-256: `dc7c629220fd0f3b23a6fcb6f156537d0d1d38b376550606235aeda7729dd2b4`
- A28 SHA-256: `d721e605d126e555f92b43200aa73b999ed4e2b61397edd1e98c0d9fddb36575`

For `PerturbationType = Eigenprofile`, the initial salt mole fraction is
`x(x,y) = SaltMoleFraction + amplitude * A28(y) * cos(2*pi*x/wavelength)`.
`Problem.Amplitude` is therefore an additive mole-fraction coefficient; it does
not change A28 or the wavelength.

With `SaltMoleFraction = 0.0036`, `Problem.Amplitude = 1e-5`, normalized A28
maximum `1`, and the cosine extrema represented on every even horizontal grid,
the initial concentration range is `0.00359` to `0.00361`. The generated-input
audit confirms that each new input differs from its own `1e-6` reference only
in `Problem.Amplitude`, `Problem.Name`, and the authorized 6T `TimeLoop.TEnd`.
