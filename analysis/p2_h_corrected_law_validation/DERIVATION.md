# Derivation Package

## Target

Independently derive and verify the corrected finite-\(Da\) top-growth law for
a prescribed downward diffusive salt flux at the crust bottom.

## Status

**COHERENT AS STATED**, subject to the explicitly listed linearization,
quasi-steady, constant-property and boundary-condition assumptions.

## Invariant Object

The upward total salt flux

\[
J=S_{\mathrm{ref}}-\rho_l\mathcal D\,\frac{dc}{dz},
\qquad
S_{\mathrm{ref}}=\frac{C_{\mathrm{eq}}j}{1-C_{\mathrm{eq}}}.
\]

## Assumptions

- \(c=C-C_{\mathrm{eq}}\) is small.
- \(Pe_\ell=j\ell/(\rho_l\mathcal D)\ll1\); perturbation advection is neglected.
- The reaction–diffusion profile is quasi-steady.
- \(\epsilon\), \(D_{\mathrm{eff}}\), \(a_v\), and \(k_r\) are constant.
- The bottom prescribes the diffusive flux only:
  \(\rho_l\epsilon D_{\mathrm{eff}}C'(0)=B_{\mathrm{soil}}\).
- The top uses the same finite-rate precipitation and salt balance as the ALE
  surface closure after linearization.

## Notation

\[
\mathcal D=\frac{\epsilon D_{\mathrm{eff}}}{1-C_{\mathrm{eq}}},
\quad
\ell=\sqrt{\frac{\mathcal D}{a_vk_r}},
\quad
H=\frac{h}{\ell}=\sqrt{Da},
\]

\[
k_s=\beta(1-\epsilon_t)k_r,
\quad
X=\frac{k_s\ell}{\mathcal D}=\sqrt{Da^*},
\quad
u=\frac{B_{\mathrm{soil}}}{C_{\mathrm{eq}}j}.
\]

All \(Da\), \(Da^*\), \(X\), \(H\), and \(u\) are dimensionless.
\(j\), \(B_{\mathrm{soil}}\), \(J\), and \(G_{\mathrm{top}}\) have units
\(\mathrm{kg\,m^{-2}\,s^{-1}}\). Both \(k_s\) and \(k_r\) have units
\(\mathrm{m\,s^{-1}}\).

## Derivation Strategy

Solve the linear reaction–diffusion equation, impose the bottom Neumann
condition, then impose the top Robin balance and eliminate the two integration
constants.

## Derivation Map

1. Linearized salt transport gives \(J=S_{\mathrm{ref}}-\rho_l\mathcal D c'\).
2. Salt conservation with \(R=\rho_la_vk_rc\) gives \(c''=c/\ell^2\).
3. The bottom condition fixes the odd integration constant.
4. The top balance fixes the even integration constant and \(c(h)\).
5. \(G_{\mathrm{top}}=\rho_lk_sc(h)\) gives the velocity law.

## Main Derivation

With \(\zeta=z/\ell\),

\[
c=a\cosh\zeta+b\sinh\zeta.
\]

The bottom condition gives

\[
\rho_l\epsilon D_{\mathrm{eff}}c'(0)=B_{\mathrm{soil}}
\quad\Longrightarrow\quad
b=\frac{uS_{\mathrm{ref}}}{K},
\qquad
K=\frac{\rho_l\mathcal D}{\ell}.
\]

At the top,

\[
S_{\mathrm{ref}}-Kc_\zeta(H)=XKc(H).
\]

Solving for the top excess concentration,

\[
c(h)=\frac{S_{\mathrm{ref}}}{K}
\frac{\cosh H-u}{\sinh H+X\cosh H}.
\]

Therefore

\[
G_{\mathrm{top}}
=XKc(h)
=S_{\mathrm{ref}}
\frac{X[1-u\,\mathrm{sech}(H)]}{X+\tanh H}.
\]

Since \(G_{\mathrm{top}}=\rho_{\mathrm{cr}}(1-\epsilon_t)v_t\),

\[
\boxed{
v_{\mathrm{corr}}=
\frac{C_{\mathrm{eq}}j}
{\rho_{\mathrm{cr}}(1-\epsilon_t)(1-C_{\mathrm{eq}})}
\frac{X[1-u\,\mathrm{sech}(\sqrt{Da})]}
{X+\tanh(\sqrt{Da})}}.
\]

At fixed \(Da\), division by the \(u=0\) solution yields the exact prediction
within this boundary-value problem:

\[
\boxed{
\frac{v_{\mathrm{corr}}(u)}{v_{\mathrm{corr}}(0)}
=1-u\,\mathrm{sech}(\sqrt{Da})}.
\]

## Remarks and Interpretation

- Large \(Da\) screens the top from the lower boundary exponentially in
  \(\sqrt{Da}\).
- \(u=0\) is a zero diffusive-gradient bottom, not the finite-\(Da\) Prat
  Dirichlet condition \(C(0)=C_{\mathrm{eq}}\).
- The independent numerical BVP uses a zero initial guess and does not use the
  closed-form solution for calibration.

## Boundaries and Non-Claims

- The formula does not include transient dissolved-salt storage, perturbation
  advection, evolving porosity, or changing \(Da\).
- ALE agreement is not implied by the mathematical BVP verification.
- A diffusive Neumann condition does not by itself prescribe the advective
  inflow concentration; the ALE boundary implementation must report this
  separately.

## Open Risks

- For small \(Da\), quasi-steady relaxation and moving-boundary time scales may
  not separate.
- For high \(u\), the lower part of the crust may become undersaturated and
  dissolve, invalidating constant-property assumptions.
- The continuum reaction layer can be smaller than the crystal scale.
