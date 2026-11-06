# P1-R1 Conservation Derivation

## 1. Definitions and exact binary conversions

Let (x_s) be the NaCl mole fraction, (x_w=1-x_s), (M_s) and (M_w) the component molar masses, and

\[
\bar M=x_wM_w+x_sM_s.
\]

The NaCl mass fraction is

\[
w_s=\frac{x_sM_s}{\bar M}, \qquad w_w=1-w_s.
\]

For liquid density \(\rho(x_s)\), the liquid molar density is

\[
c=\frac{\rho}{\bar M}.
\]

With liquid content \(\theta=\phi S_w\), the exact bulk-volume inventories are:

\[
n_s=\theta c x_s,\qquad n_w=\theta c(1-x_s)
\quad [\mathrm{mol\,m^{-3}_{bulk}}],
\]

\[
m_s=\theta\rho w_s,\qquad m_w=\theta\rho(1-w_s)
\quad [\mathrm{kg\,m^{-3}_{bulk}}].
\]

They obey \(m_i=M_i n_i\) exactly.

## 2. Why post-processing multiplication is not an exact conversion

The implemented water-like residual is

\[
R_0=\partial_t\theta+\nabla\cdot\mathbf q.
\]

The physical total-liquid mass residual would be

\[
R_m=\partial_t(\theta\rho)+\nabla\cdot(\rho\mathbf q)
=\rho R_0+\theta\partial_t\rho+\mathbf q\cdot\nabla\rho.
\]

Thus \(R_m\ne\rho R_0\) whenever density changes materially along the transient or spatial flow.

Likewise, the implemented salt-like residual is

\[
R_s^{\mathrm{norm}}=\partial_t(\theta x_s)
+\nabla\cdot(x_s\mathbf q-D_{\mathrm{eff}}\nabla x_s).
\]

Even if a molar salt flux is provisionally written as

\[
\mathbf N_s=c x_s\mathbf q-cD_{\mathrm{eff}}\nabla x_s,
\]

the physical residual \(\partial_t(\theta c x_s)+\nabla\cdot\mathbf N_s\) contains additional derivatives of \(c(x_s)\). It is not \(cR_s^{\mathrm{norm}}\). The deterministic audit quantifies the storage correction \(x_s c'(x_s)/c\): approximately \(-2.69\times10^{-4}\) at the initial Case 36877 concentration and \(-2.72\times10^{-2}\) at \(x_s=0.1\). Small at the initial state does not mean identically zero, and the discrepancy grows with concentration.

Therefore a diagnostic-only rescaling is admissible only as an explicitly approximate, constant-property reporting convention. It cannot close a precipitation mass balance.

## 3. Three candidate formulations

### A. Keep the normalized equations; convert diagnostics only

**Benefit:** zero change to the existing solution, perturbation, linear-stability comparison and Case 36877 history.

**Limitation:** the converted outputs do not satisfy an exact physical component balance when \(\rho\) or \(\bar M\) varies. This option cannot support conservative coupling to a growing salt crust.

**Permitted claim:** “normalized-equation conservation, with approximate reference-state conversion.”

**Forbidden claim:** “H2O/NaCl mass or molar conservation.”

### B. Molar component balances

Retain pressure and salt mole fraction as primary variables, but use physical storage:

\[
\partial_t(\theta cx_w)+\nabla\cdot\mathbf N_w=0,
\]

\[
\partial_t(\theta cx_s)+\nabla\cdot\mathbf N_s=0.
\]

Advective component fluxes are \(cx_i\mathbf q\). The diffusive terms require an explicit reference-system choice:

- **Molar-averaged Fick:** \(\mathbf J_s^n=-cD_{\mathrm{eff}}\nabla x_s\), \(\mathbf J_w^n=-\mathbf J_s^n\).
- **Mass-averaged Fick:** \(\mathbf j_s^m=-\rho D_{\mathrm{eff}}\nabla w_s\), \(\mathbf j_w^m=-\mathbf j_s^m\), followed by \(\mathbf J_i^n=\mathbf j_i^m/M_i\).

These formulations are not interchangeable: in a mass-averaged frame the two component molar diffusive fluxes do not generally sum to zero.

The candidate \(\mathbf N_s=cx_s\mathbf q-cD_{\mathrm{eff}}\nabla x_s\) is a consistent molar component flux only if (i) the advective velocity represented by \(\mathbf q\) is the molar-average liquid velocity, or an incompressible/equal-partial-molar-volume approximation makes it equivalent; (ii) diffusion is defined in the molar-average reference frame; (iii) the water diffusive molar flux is the equal and opposite flux; and (iv) the spatially varying \(c\) remains inside the divergence operator. It is not obtained by multiplying the divergence of the old normalized flux by a local \(c\).

### C. Total-mass plus salt-mass balances — recommended

Retain pressure and (x_s) as numerical primary variables, but convert (x_s\mapsto w_s) inside the residual:

\[
\partial_t(\theta\rho)+\nabla\cdot(\rho\mathbf q)=0,
\]

\[
\partial_t(\theta\rho w_s)
+\nabla\cdot(\rho w_s\mathbf q+\mathbf j_s^m)=0,
\]

with the standard mass-averaged closure

\[
\mathbf j_s^m=-\rho D_{\mathrm{eff}}\nabla w_s,
\qquad \mathbf j_w^m=-\mathbf j_s^m.
\]

Water mass follows by subtracting salt mass from total liquid mass. This basis matches the existing P0 Prat interface, where crust storage and atmospheric loss are already expressed as component masses.

**Scientific consequence:** changing \(\nabla x_s\) to \(\nabla w_s\), and inserting density in storage/advection, changes the PDE. The original Case 36877 solution and growth rate are therefore comparison baselines, not expected regression identities. Preserving \(-D\nabla x_s\) while labelling it a salt mass flux would instead require a new empirical-flux definition and cannot be called standard mass-averaged Fick diffusion.

The scalar constitutive function \(D_{\mathrm{eff}}(\phi,S_w)\) can remain unchanged. The required change from a mole-fraction gradient to a mass-fraction gradient is a change in reference-system representation, not parameter tuning; it must nevertheless be treated as a governing-equation change in the Case 36877 comparison.

### Side-by-side assessment

| Criterion | A: diagnostics only | B: molar components | C: total mass + salt mass |
|---|---|---|---|
| Mathematical consistency | Exact only for normalized equations; physical conversion approximate | Exact if velocity and molar diffusion reference are declared consistently | Exact with mass-average flux and component mass storage |
| Physical conservation | No, for variable \(\rho,c\) | Yes | Yes |
| C++ change | Diagnostics only | Custom storage, advection, diffusion, BC and diagnostics | Custom storage, advection, diffusion, BC and diagnostics |
| Existing LSA compatibility | Same PDE and solution | Same primary perturbation variable, different PDE/operator | Same primary perturbation variable, different PDE/operator |
| Unsaturated hydraulics | Identical | Darcy/capillary laws retained; composition coupling changes balance | Darcy/capillary laws retained; composition coupling changes balance |
| Prat coupling | Insufficient for exact closure | Possible after component-wise mass/mole conversion | Directly compatible with mass-based crust storage |
| Decision | Baseline/reporting only | Defensible but reference-frame sensitive | Recommended minimum physical basis |

## 4. Top boundary before precipitation

Use outward/upward flux as positive. The present normalized condition is

\[
F_0=e\ [\mathrm{m\,s^{-1}}],\qquad F_s^{\mathrm{norm}}=0.
\]

A physical boundary condition must declare what (e) means. If it is an equivalent pure-water volume demand referenced to \(\rho_{w,\mathrm{ref}}\), then

\[
F_w^m=\rho_{w,\mathrm{ref}}e,qquad
F_w^n=\frac{\rho_{w,\mathrm{ref}}e}{M_w},qquad F_s=0.
\]

For formulation C, the equation-0 total-mass flux is \(F_w^m\) and the equation-1 salt-mass flux is zero. Evaporation must not be implemented as an unqualified bulk-brine outflow \(\rho\mathbf q\), because that would advect salt through a boundary intended to be salt-impermeable.

## 5. Future conservative Prat interface

Let \(M_s^c\) and \(M_w^c\) be salt and water stored per interfacial area in the crust, and let \(j_{atm}\ge0\) be water mass transferred from the crust to the atmosphere. Conservation requires

\[
F_{s,soil}^m=\frac{dM_s^c}{dt},
\]

\[
F_{w,soil}^m=j_{atm}+\frac{dM_w^c}{dt}.
\]

For formulation C, with equation 0 = total liquid mass and equation 1 = salt mass,

\[
F_0^m=F_{w,soil}^m+F_{s,soil}^m,
\qquad F_1^m=F_{s,soil}^m.
\]

For a molar formulation, convert the two component mass fluxes separately:

\[
F_i^n=F_i^m/M_i.
\]

There is no generally exact single factor that maps these future component fluxes back to the two current normalized equations when density and composition vary.

For comparison only, a frozen reference state \((\rho_r,c_r,x_r,w_{s,r})\) would give the approximate original-unit mapping

\[
F_0^{\mathrm{norm}}\approx
\frac{F_w^m+F_s^m}{\rho_r},
\qquad
F_1^{\mathrm{norm}}\approx
\frac{F_s^m}{c_rM_s}
=\frac{x_r}{\rho_rw_{s,r}}F_s^m.
\]

Before precipitation this reduces to \(F_0^{\mathrm{norm}}=e\), \(F_1^{\mathrm{norm}}=0\) under the stated pure-water reference-density convention. After precipitation it is not an exact conservative boundary law and must not be used for Prat coupling.

## 6. Scientific conclusion

The minimal defensible target is **formulation C implemented inside the existing custom Richards architecture**, not a new RichardsNC solver and not a diagnostic-only relabelling. Pressure and mole fraction can remain primary variables, so the existing initial perturbation and LSA input representation are preserved. What changes is the conservation law evaluated from them.

This recommendation is conditional on two decisions being recorded before code changes:

1. confirm mass-averaged Fick diffusion in \(w_s\) as the physical closure;
2. confirm whether the imposed evaporation number is a pure-water volume demand, water mass flux, or another experimentally derived quantity.

Until those choices are explicit, implementation is scientifically ambiguous and P1-R1 should stop at the audited change plan.
