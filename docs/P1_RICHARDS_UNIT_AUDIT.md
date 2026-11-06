# P1-R1 Richards–Salt Unit Audit

## Audit decision

**The frozen Richards–salt implementation conserves two normalized quantities, not physical H2O/NaCl mass or mole components.** Its discrete balances are internally useful and can be exactly conservative in their own units, but their inventories and fluxes must not be labelled `kg`, `mol`, `kg/s`, or `mol/s`.

This is an evidence–claim mismatch, not merely a metadata problem. Multiplying the existing diagnostic output by density or molar density does not turn the governing residual into a physical component balance when composition-dependent density and mean molar mass are retained.

No frozen source, material law, boundary value, mesh, initial perturbation, or Case 36877 setting was changed during this audit. No production simulation was started.

## Provenance and scope

The primary audited source is the immutable Paper-1 snapshot:

- `model/transport_wp8_background_v2/source/appl/simpleRichardsnc/model.hh`
- `model/transport_wp8_background_v2/source/appl/simpleRichardsnc/main.cc`
- `model/transport_wp8_background_v2/source/appl/simpleRichardsnc/transportdiagnostics.hh`
- Paper-1 commit: `158395f3d25a29680dea5bcd0d3a24e52a28b4cd`

The local source manifest was checked without finding a hash mismatch. Exact file hashes and the read-only upstream references are recorded in `analysis/p1_r1_unit_audit_20261106/source_manifest.json`.

The upstream comparison used DuMuX commit `63568e6ad0016124fa89ecd932adb613118527bf`, in particular:

- `dumux/porousmediumflow/compositional/localresidual.hh`
- `dumux/porousmediumflow/richardsnc/model.hh`
- `dumux/porousmediumflow/richardsnc/volumevariables.hh`
- `dumux/flux/{box/fickslaw.hh,fickslaw_fwd.hh,referencesystemformulation.hh}`

## What the current code actually solves

Define liquid saturation (S_w), porosity \(\phi\), salt mole fraction (x_s), Darcy phase velocity \(\mathbf q\), and

\[
\theta=\phi S_w.
\]

The custom local residual implements, in continuum notation,

\[
\frac{\partial \theta}{\partial t}+\nabla\cdot\mathbf q=0,
\]

\[
\frac{\partial (\theta x_s)}{\partial t}
+\nabla\cdot\left(x_s\mathbf q-D_{\mathrm{eff}}\nabla x_s\right)=0.
\]

The source evidence is direct:

- `normalizedStorage()` returns \(\phi S_w\) and \(\phi S_w x_s\).
- `transportFlux()` returns area-integrated \(q_n\) and \(q_n x_{s,\mathrm{up}}-D_{\mathrm{eff}}\nabla x_s\cdot\mathbf n\).
- Neither storage nor flux contains liquid mass density, liquid molar density, or component molar mass.
- `molarRhoW(x)` is defined in `main.cc`, but is not used by the residual or the transport diagnostics.
- `transportdiagnostics.hh` correctly describes its salt inventory as `phi*Sw*x` and its flux unit as `m/s`; the bottom Dirichlet contribution is inferred as a reaction term rather than independently reconstructed.

### Dimensional meaning

| Quantity | Current expression | Local unit | Domain-integrated meaning |
|---|---:|---:|---:|
| Equation 0 storage | \(\theta\) | 1 | liquid-volume-equivalent, \(\mathrm{m^3}\) |
| Equation 0 flux | \(q_n\) | \(\mathrm{m\,s^{-1}}\) | \(\mathrm{m^3\,s^{-1}}\) after face integration |
| Equation 1 storage | \(\theta x_s\) | 1 | mole-fraction-weighted liquid volume, \(\mathrm{m^3}\) |
| Equation 1 flux | \(q_nx_s-D_{\mathrm{eff}}\nabla x_s\cdot\mathbf n\) | \(\mathrm{m\,s^{-1}}\) | weighted volume rate after face integration |

For a column diagnostic divided by horizontal area, both inventories have equivalent-depth units of metres. Equation 1 is **not** salt moles, salt mass, or salt volume.

## Initial and boundary conditions

- Bottom: hydrostatic pressure and fixed salt mole fraction are imposed as Dirichlet data for both equations.
- Top and lateral boundaries: Neumann conditions are used.
- Top evaporation: `EvaporationRate` in `m/s` is inserted directly into equation 0; equation 1 is zero.
- Initial concentration perturbations and eigenprofiles are prescribed in mole fraction and are not themselves unit-inconsistent.

Consequently, the present top condition is consistent with the normalized equation as an equivalent liquid-volume demand. It is not yet a declared H2O mass or molar flux.

## Material-property audit

The following closures are present and dimensionally interpretable:

- Darcy velocity: permeability, relative permeability, viscosity, pressure gradient, density and gravity produce \(\mathbf q\) in \(\mathrm{m\,s^{-1}}\).
- Salt-dependent density: mole fraction is converted to mass fraction before evaluating \(\rho(x_s)\).
- Liquid molar density: \(c=\rho/\bar M\) is implemented by `molarRhoW(x)`, with \(\bar M=(1-x_s)M_w+x_sM_s\), but is unused in the balance law.
- Effective diffusion: \(D_{\mathrm{eff}}=\phi^{1.5}S_w^{3.5}D_{\mathrm{mol}}\), in \(\mathrm{m^2\,s^{-1}}\).

The material-property functions therefore do not repair the missing density factors in storage and transport.

## Comparison with pinned upstream DuMuX

The upstream `CompositionalLocalResidual` consistently selects one physical basis:

- with `UseMoles=true`, storage is molar density times mole fraction and advection is molar density times mole fraction times phase mobility/flux;
- with `UseMoles=false`, storage is mass density times mass fraction and advection is the corresponding mass flux;
- diffusive flux is converted according to its declared mass-averaged or molar-averaged reference system;
- `ReplaceCompEqIdx` can replace one component equation with a total balance.

Upstream RichardsNC defaults to a compositional local residual and exposes these properties. The frozen custom model, however, defines its own model traits, volume variables, storage and flux. It does not expose the complete fluid-system interface expected by the upstream compositional residual. **Setting `UseMoles=true` in the current custom model would therefore not fix the equations.** Importing the full RichardsNC stack would be a framework replacement and is outside P1-R1 scope.

## Deterministic conversion check

`analysis/p1_r1_unit_audit_20261106/conversion_audit.py` reproduces the binary H2O–NaCl conversions and asserts component mass identities. Its committed JSON output gives, for the Case 36877 initial state \(x_s=0.0036\):

- \(\bar M=0.01814559408\ \mathrm{kg\,mol^{-1}}\),
- \(w_s=0.011594775\),
- \(\rho=1007.8280\ \mathrm{kg\,m^{-3}}\),
- \(c=55541.1984\ \mathrm{mol\,m^{-3}}\).

For the explicitly conditional interpretation of \(e=3.065544359\times10^{-9}\ \mathrm{m\,s^{-1}}\) as a pure-water volume demand referenced to \(1000\ \mathrm{kg\,m^{-3}}\):

\[
F_w^m=3.065544359\times10^{-6}\ \mathrm{kg\,m^{-2}\,s^{-1}},
\]

\[
F_w^n=1.7030801994\times10^{-4}\ \mathrm{mol\,m^{-2}\,s^{-1}}.
\]

These are boundary conversions under a declared reference-density assumption; they are not evidence that the interior normalized residual conserves physical components.

## Direct answers from the audit

1. **What is currently conserved?** Normalized liquid content \(\theta\) and normalized salt-weighted content \(\theta x_s\), with matching area-integrated normalized fluxes.
2. **Are the current diagnostics wrong?** Their arithmetic is appropriate for the implemented normalized residual, but any interpretation as physical mass/moles is wrong.
3. **Can a diagnostic-only conversion create exact conservation?** No, except under a constant-density, constant-mean-molar-mass approximation with consistently scaled sources and boundaries.
4. **Can upstream `UseMoles` be enabled as a one-line fix?** No; the custom residual bypasses the upstream compositional machinery.
