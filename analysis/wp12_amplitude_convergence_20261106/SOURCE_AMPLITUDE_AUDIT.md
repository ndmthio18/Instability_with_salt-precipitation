# `Problem.Amplitude` source audit

- Source: `model/transport_v1/source/appl/simpleRichardsnc/main.cc`
- Source SHA-256: `dc7c629220fd0f3b23a6fcb6f156537d0d1d38b376550606235aeda7729dd2b4`
- A28 SHA-256: `d721e605d126e555f92b43200aa73b999ed4e2b61397edd1e98c0d9fddb36575`

The constructor reads `Problem.Amplitude` directly into `amplitude_`. For an
`Eigenprofile` perturbation, the initial salt mole fraction is

`x(x,y) = SaltMoleFraction + amplitude_ * A28(y) * cos(2*pi*x/wavelength)`.

The wavelength is set independently as the physical domain width divided by
`Problem.WaveNumber`; changing `Problem.Amplitude` therefore does not change the
lateral wavelength or the A28 vertical profile.

The frozen A28 table has minimum `-8.642409378252e-11` and maximum `1`. All
three even horizontal grids contain both cosine extrema. With
`SaltMoleFraction = 0.0036` and `Problem.Amplitude = 1e-4`, the evaluated initial
range is exactly `0.0035` to `0.0037` to floating-point precision. This is
positive and far below the registered solubility threshold `0.1`.

The generated input audit confirms that each new input differs from its own
`1e-6` reference only in `Problem.Amplitude`, `Problem.Name` and the authorized
6T `TimeLoop.TEnd`.
