# P1 implementation and reproducibility

## Architecture

The isolated module is `model/p1_compositional_baseline/source`. Both executables inherit upstream `TwoPNC`, use cell-centred TPFA, tensor-product YaspGrid meshes and `UseMoles=true`. One-dimensional and two-dimensional targets share the same problem, fluid system and spatial parameters.

`P1FluidSystem` inherits upstream `BrineAir` and overrides only the liquid H2O--NaCl binary molecular diffusion coefficient. `P1SpatialParams` supplies homogeneous permeability, porosity, temperature, gravity and Van Genuchten--Mualem relations. `P1Problem` supplies initial states, phase switching, boundaries and diagnostics.

No `2pncmin` tag, mineral volume variable, precipitation source, solid saturation or crust state is compiled. The archived v12 model was used only as an API/physics reference.

## Conservation diagnostics

At every accepted step, the code integrates

\[
M_k=\int_\Omega \phi\sum_\alpha S_\alpha c_\alpha x_\alpha^k\,dV
\]

and all boundary component fluxes. Dirichlet boundary fluxes are evaluated with the same face operator at the accepted state; Neumann fluxes are integrated from the prescribed density. The primary residual is \(M_k(t)-M_k(0)+\sum_n\Delta t_nF_{k,n}\).

For V7, `WriteCellStates=true` writes initial and final control-volume geometry and thermodynamic state. A separate Python implementation reconstructs all three inventories and combines them with the analytic V1/V2 boundary flux. It does not read the solver's inventory columns.

The audit also reports phase-state counts, saturation/composition bounds, top salt concentration and density, mode-28 amplitude, and reconstructed liquid/gas Darcy velocities.

## Environment evidence

The exact module commits, compiler, source hashes and executable hashes are in `model/p1_compositional_baseline/manifests/dependencies_20261106.json`. Shared DUNE/DuMuX sources were not modified. SuperLU was disabled only for this module because the host's SuperLU 7 headers are incompatible with the DUNE 2.9 wrapper and the selected iterative solver does not require it.

The first runtime attempt exposed an OpenMPI/UCX allocator deadlock during `MPI_Init`. The recorded single-rank TCP/`ob1` environment removed the deadlock; the failed attempt remains under the remote versioned work root.

## Case 36877 mapping

The production inputs preserve the 0.60 m by 0.20 m domain, 600 by 200 mesh, evaporation rate, water-table depth, porosity, permeability, capillarity, residual saturation, diffusion coefficient, initial salt fraction, gravity, mode 28 and amplitude \(10^{-6}\). The frozen eigenprofile is referenced by repository path and SHA-256.

The following differences are intentional and must remain visible: physical molar component balances instead of normalized \(\phi S_wx_s\); mobile gas and dissolved air; upstream `BrineAir` density, viscosity and phase equilibrium; mass-averaged Fick diffusion converted by upstream code to molar equation units; and TPFA/nonlinear-solver differences.
