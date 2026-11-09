# Saturated soil-crust interface without dissolution

This isolated control asks whether the existing xNaCl=0.09, 80-soil-cell,
64-ALE-cell system can satisfy all of the following at once after activation:

1. the soil-side and crust-side interface concentration equals `C_eq`;
2. internal dissolution is suppressed;
3. the interface phase-change source is exactly zero;
4. the existing gas-film, soil equations, finite-rate top growth, ALE operator,
   activation transfer, and timestep remain unchanged.

`P2H.InterfaceMode=saturated_no_source_strict` fixes the ALE boundary face at `C_eq`,
advances ALE with no added phase inventory, applies the resulting conservative
liquid and salt fluxes to DuMux, and checks the reconstructed soil surface
concentration. The accepted step requires
`abs(C_soil-C_eq)<5e-13` and the existing `Sw>=0.99` screen.

`P2H.InterfaceMode=crust_saturated_soil_undersaturated` instead keeps only the
crust-side face at `C_eq`. It accepts a soil-side concentration jump satisfying
`-P2H.MaxSoilSideUndersaturation <= C_soil-C_eq <= 5e-13`, while retaining the
same conservative liquid/salt flux and zero phase source. The first exploratory
bound is `1e-4` in NaCl mass fraction. This jump is a declared numerical closure,
not a calibrated interfacial mass-transfer law.

An incompatible step is rejected transactionally and terminates without
timestep subdivision. This is deliberate: shrinking `dt` until the endpoint
change falls below the concentration tolerance would hide a continuous-time
compatibility failure. The rejected trial is logged in
`interface_compatibility.csv`; it never enters accepted inventories.

This mode does not project concentrations, clip a salt deficit, add a reservoir,
reduce evaporation, dissolve interface solid, or settle the skeleton. A failed
pilot is therefore a scientific result, not a solver result to be tuned away.

