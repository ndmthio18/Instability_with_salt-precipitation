"""Prat (2026) porous salt growth law as a standalone surface subgrid module.

Reference: Marc Prat, *The Evaporation-Driven Growth Law of Porous Salt
Structures*, Langmuir 42 (2026) 21269-21282, doi:10.1021/acs.langmuir.6c02211.
Every equation used here is taken from ``docs/PRAT_PHYSICS_AUDIT.md``, which was
written directly against the publisher PDF (SHA-256
82E4F450481F7185AFAF3F1B1C2A7A60C25C16FF04B0A570FCB2C989EB23AB3A). Equation
numbers below refer to that paper. Nothing is reconstructed from memory.

Imported verbatim from the audit
--------------------------------
Dimensionless groups, audit section 3::

    Pe      = U_evap h / (eps_c D_s_eff),      U_evap = j / rho
    Da      = h^2 (1 - C_eq) a_v k_r / (eps_c D_s_eff)
    Da_top  = beta (1 - eps_c) (1 - C_eq) h k_r / (eps_c D_s_eff)
    Da_star = beta^2 (1 - eps_c)^2 (1 - C_eq) k_r / (eps_c D_s_eff a_v)

General finite-Da growth law, equation (28)::

    v_s = C_eq j / (rho_cr (1-eps_c)(1-C_eq))
          * sqrt(Da*) tanh(sqrt(Da)) / (1 + sqrt(Da*) tanh(sqrt(Da)))

Large-Da limit, equations (29)/(32), obtained by tanh(sqrt(Da)) -> 1.

Local surface generalisation, equation (36): the same scalar factor multiplies
the normal component of the evaporation flux.

Moving-interface balances, equations (6) and (7), audit section 2.2.

Internal reaction, equations (1)-(2): the solid production rate density is
``r_i = a_v rho k_r (C - C_eq)`` while the transport equation for C carries the
extra ``(1 - C)`` factor of equation (5). Both are kept distinct here; the
``(1 - C_eq)`` factor therefore appears in ``Da`` but not in ``r_i``.

Derived in this module (labelled, not attributed to the paper)
--------------------------------------------------------------
The audit gives the boundary conditions (14)-(15) ``C*(0) = 0``,
``C*(1) = C*_top`` with ``C* = C - C_eq``, and the quasi-steady, Pe << Da,
uniform-property reduction. That reduction is the linear two-point problem

    d^2 C* / dz*^2 = Da C*,   C*(0) = 0,

whose solution is ``C*(z*) = C*_top sinh(sqrt(Da) z*) / sinh(sqrt(Da))``.

Combining equations (6) and (7) with ``C ~ C_eq`` (audit section 2.2) gives the
paper's equation (8),

    rho eps_c D_s_eff dC/dz|_h = C_eq j - rho_cr (1-eps_c)(1-C_eq) v_s,

up to the negligible ``C_eq rho_v v_s`` term. Substituting the analytic profile
and requiring that the result equal equation (28) fixes the top supersaturation
uniquely::

    C*_top = C_eq Pe tanh(sqrt(Da)) / ( sqrt(Da) (1 + sqrt(Da*) tanh(sqrt(Da))) )

This expression is **derived here**, not quoted. ``verification/`` checks that
it is simultaneously consistent with (a) equation (28), (b) equation (8)
evaluated on the analytic profile, and (c) the outer-surface kinetic relation
``v_s = beta k_r rho C*_top / rho_cr``, which is the form of the paper's
equation (20) compatible with equation (28). The three routes agreeing to
round-off is the module's primary internal-consistency gate.

What this module deliberately does NOT do
-----------------------------------------
* It is not a nucleation model. ``CrustState`` stays inert until an explicit
  seed/precipitation event is declared (P2-A item 6).
* It never applies ``min(Prat demand, available soil flux)`` as a model.
  ``supply_limitation_diagnostic`` reports the ratio and can raise a declared
  model-domain exit; it never silently clips (audit section 7).
* It does not partition solid salt into "top" and "internal" using ``1 - f``.
  Total solid is the tracked state; the top and internal contributions are
  reported separately only because both are computed from the same reduced
  solution (audit section 6).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field, replace
from enum import Enum

from .thermo import BrineDensityLaw


class GrowthLaw(str, Enum):
    """Which form of the Prat growth law to evaluate."""

    FINITE_DA = "finite_da"   # equation (28)
    LARGE_DA = "large_da"     # equations (29)/(32)


class CrustModelDomainExit(RuntimeError):
    """Raised when the Prat closure leaves its declared validity envelope."""


# --- numerically safe elementary helpers -------------------------------------

_SMALL = 1.0e-8


def _tanh_over_u(u: float) -> float:
    """tanh(u)/u with the correct u -> 0 limit."""
    if u < _SMALL:
        return 1.0 - u * u / 3.0
    return math.tanh(u) / u


def _tanh_half_over_u(u: float) -> float:
    """tanh(u/2)/u = (cosh u - 1)/(u sinh u), with the u -> 0 limit 1/2."""
    if u < _SMALL:
        return 0.5 - u * u / 24.0
    return math.tanh(0.5 * u) / u


def _u_over_sinh(u: float) -> float:
    """u / sinh(u), overflow-safe, with the u -> 0 limit 1."""
    if u < _SMALL:
        return 1.0 - u * u / 6.0
    if u > 700.0:
        return 0.0
    return u / math.sinh(u)


def _sinh_ratio(u: float, zstar: float) -> float:
    """sinh(u z*) / sinh(u), overflow-safe, with the u -> 0 limit z*."""
    if u < _SMALL:
        return zstar
    if u > 50.0:
        # Divide numerator and denominator of (e^a - e^-a)/(e^u - e^-u) by e^u.
        return math.exp(u * (zstar - 1.0)) - math.exp(-u * (zstar + 1.0))
    return math.sinh(u * zstar) / math.sinh(u)


def _require_positive_height(height: float) -> None:
    if not math.isfinite(height) or height <= 0.0:
        raise CrustModelDomainExit(
            "the Prat continuum crust problem requires h > 0; at h = 0 the "
            "domain, the diffusion-time argument and the boundary-layer "
            "argument all fail. Declare a seed event first (P2-A item 6)."
        )


# --- parameters ---------------------------------------------------------------


@dataclass(frozen=True)
class CrustParameters:
    """Declared physical parameters of the external porous salt structure.

    All parameters are inputs. None of them is fitted here.

    Attributes
    ----------
    equilibrium_mass_fraction:
        ``C_eq``, dissolved NaCl *mass* fraction at halite equilibrium [-].
    crust_porosity:
        ``eps_c`` [-], assumed uniform and constant (audit section 3,
        reduction 3).
    crystal_density:
        ``rho_cr`` [kg/m^3], halite density.
    specific_area:
        ``a_v`` [1/m], internal crystal area per bulk crust volume.
    reaction_coefficient:
        ``k_r`` [m/s], precipitation/dissolution coefficient.
    effective_salt_diffusivity:
        ``D_s_eff`` [m^2/s].
    beta:
        Upscaled outer-surface geometric factor of equation (20). The paper's
        ideal cubic construction gives ``beta = 2`` (audit section 3).
    solution_density:
        ``rho`` [kg/m^3] of the pore solution, evaluated once at ``C_eq``
        (audit section 3, reduction 1). ``None`` uses the frozen P1 brine law.
    vapour_density:
        ``rho_v`` [kg/m^3], the swept-volume vapour term of equation (6). The
        audit requires it to be stated rather than hidden in a sign convention.
    crust_saturation:
        ``S_c`` [-]. The Prat closure assumes a capillary-supplied, liquid
        saturated crust; any value below ``connectivity_threshold`` triggers a
        declared model-domain exit instead of continued growth.
    connectivity_threshold:
        Declared hydraulic-connectivity limit (audit section 5.4).
    """

    equilibrium_mass_fraction: float
    crust_porosity: float
    crystal_density: float
    specific_area: float
    reaction_coefficient: float
    effective_salt_diffusivity: float
    beta: float = 2.0
    solution_density: float | None = None
    vapour_density: float = 0.0
    crust_saturation: float = 1.0
    connectivity_threshold: float = 0.99

    def __post_init__(self) -> None:
        if not (0.0 < self.equilibrium_mass_fraction < 1.0):
            raise ValueError("C_eq must lie strictly in (0, 1)")
        if not (0.0 < self.crust_porosity < 1.0):
            raise ValueError("crust porosity must lie strictly in (0, 1)")
        for name in (
            "crystal_density",
            "specific_area",
            "reaction_coefficient",
            "effective_salt_diffusivity",
            "beta",
        ):
            value = getattr(self, name)
            if not math.isfinite(value) or value <= 0.0:
                raise ValueError(f"{name} must be positive and finite")
        if self.vapour_density < 0.0:
            raise ValueError("vapour density must be non-negative")
        if not (0.0 <= self.crust_saturation <= 1.0):
            raise ValueError("crust saturation must lie in [0, 1]")

    @property
    def rho(self) -> float:
        """Pore-solution density, evaluated once at ``C_eq``."""
        if self.solution_density is not None:
            return self.solution_density
        return BrineDensityLaw().density(self.equilibrium_mass_fraction)

    @property
    def damkohler_star(self) -> float:
        """``Da*``; height-independent by construction."""
        return (
            self.beta**2
            * (1.0 - self.crust_porosity) ** 2
            * (1.0 - self.equilibrium_mass_fraction)
            * self.reaction_coefficient
            / (
                self.crust_porosity
                * self.effective_salt_diffusivity
                * self.specific_area
            )
        )

    def damkohler(self, height: float) -> float:
        """``Da`` at the given crust height."""
        _require_positive_height(height)
        return (
            height**2
            * (1.0 - self.equilibrium_mass_fraction)
            * self.specific_area
            * self.reaction_coefficient
            / (self.crust_porosity * self.effective_salt_diffusivity)
        )

    def damkohler_top(self, height: float) -> float:
        """``Da_top`` at the given crust height."""
        _require_positive_height(height)
        return (
            self.beta
            * (1.0 - self.crust_porosity)
            * (1.0 - self.equilibrium_mass_fraction)
            * height
            * self.reaction_coefficient
            / (self.crust_porosity * self.effective_salt_diffusivity)
        )

    def peclet(self, height: float, evaporation_flux: float) -> float:
        """``Pe`` from the evaporation mass flux ``j`` [kg/m^2/s]."""
        _require_positive_height(height)
        u_evap = evaporation_flux / self.rho
        return u_evap * height / (
            self.crust_porosity * self.effective_salt_diffusivity
        )


# --- closed growth law --------------------------------------------------------


def efficiency(damkohler: float, damkohler_star: float) -> float:
    """``f(Da, Da*)`` of equation (28)."""
    if damkohler < 0.0 or damkohler_star < 0.0:
        raise ValueError("Damkohler numbers must be non-negative")
    s = math.sqrt(damkohler_star) * math.tanh(math.sqrt(damkohler))
    return s / (1.0 + s)


def efficiency_large_damkohler(damkohler_star: float) -> float:
    """``f_inf(Da*)`` of equations (29)/(32), i.e. tanh(sqrt(Da)) -> 1."""
    if damkohler_star < 0.0:
        raise ValueError("Da* must be non-negative")
    s = math.sqrt(damkohler_star)
    return s / (1.0 + s)


def growth_prefactor(params: CrustParameters, evaporation_flux: float) -> float:
    """``G = C_eq j / (rho_cr (1-eps_c)(1-C_eq))``, the f = 1 growth speed."""
    return (
        params.equilibrium_mass_fraction
        * evaporation_flux
        / (
            params.crystal_density
            * (1.0 - params.crust_porosity)
            * (1.0 - params.equilibrium_mass_fraction)
        )
    )


def growth_velocity(
    params: CrustParameters,
    height: float,
    evaporation_flux: float,
    law: GrowthLaw = GrowthLaw.FINITE_DA,
) -> float:
    """``v_s = dh/dt`` [m/s] from equation (28) or (29)/(32)."""
    da_star = params.damkohler_star
    if law is GrowthLaw.LARGE_DA:
        _require_positive_height(height)
        f = efficiency_large_damkohler(da_star)
    else:
        f = efficiency(params.damkohler(height), da_star)
    return growth_prefactor(params, evaporation_flux) * f


# --- reduced crust transport solution ----------------------------------------


def top_supersaturation(
    params: CrustParameters,
    height: float,
    evaporation_flux: float,
    law: GrowthLaw = GrowthLaw.FINITE_DA,
) -> float:
    """``C*_top = C_top - C_eq`` [-]; derived, see the module docstring."""
    da = params.damkohler(height)
    u = math.sqrt(da)
    da_star = params.damkohler_star
    pe = params.peclet(height, evaporation_flux)
    if law is GrowthLaw.LARGE_DA:
        # tanh(sqrt(Da)) -> 1 consistently in numerator and denominator.
        inv_u = 1.0 / u if u > _SMALL else 1.0
        return (
            params.equilibrium_mass_fraction
            * pe
            * inv_u
            / (1.0 + math.sqrt(da_star))
        )
    return (
        params.equilibrium_mass_fraction
        * pe
        * _tanh_over_u(u)
        / (1.0 + math.sqrt(da_star) * math.tanh(u))
    )


def top_mass_fraction(
    params: CrustParameters,
    height: float,
    evaporation_flux: float,
    law: GrowthLaw = GrowthLaw.FINITE_DA,
) -> float:
    """``C_top`` at the evaporating crust surface."""
    return params.equilibrium_mass_fraction + top_supersaturation(
        params, height, evaporation_flux, law
    )


def concentration_profile(
    params: CrustParameters,
    height: float,
    evaporation_flux: float,
    zstar: float,
    law: GrowthLaw = GrowthLaw.FINITE_DA,
) -> float:
    """``C(z*)`` on the reduced crust problem, ``z* = z/h`` in [0, 1]."""
    if not (0.0 <= zstar <= 1.0):
        raise ValueError("z* must lie in [0, 1]")
    u = math.sqrt(params.damkohler(height))
    c_star_top = top_supersaturation(params, height, evaporation_flux, law)
    return params.equilibrium_mass_fraction + c_star_top * _sinh_ratio(u, zstar)


def internal_reaction_rate(
    params: CrustParameters,
    height: float,
    evaporation_flux: float,
    law: GrowthLaw = GrowthLaw.FINITE_DA,
) -> float:
    """Depth-integrated internal precipitation ``R_int`` [kg NaCl/m^2/s].

    ``R_int = int_0^h a_v rho k_r C* dz``, which on the analytic profile is
    ``a_v rho k_r h C*_top tanh(sqrt(Da)/2)/sqrt(Da)``.
    """
    u = math.sqrt(params.damkohler(height))
    c_star_top = top_supersaturation(params, height, evaporation_flux, law)
    return (
        params.specific_area
        * params.rho
        * params.reaction_coefficient
        * height
        * c_star_top
        * _tanh_half_over_u(u)
    )


def porosity_consumption_timescale(
    params: CrustParameters,
    height: float,
    evaporation_flux: float,
    law: GrowthLaw = GrowthLaw.FINITE_DA,
) -> float:
    """Time for the internal reaction to fill the existing pore space [s].

    Equation (2) states ``-rho_cr d(eps_s)/dt = a_v rho k_r (C - C_eq)``, so a
    non-zero internal reaction necessarily changes the crust porosity. The
    reduced solution used here nevertheless holds ``eps_c`` constant (audit
    section 3, reduction 3). That is only defensible while

        t_eps = eps_c rho_cr / (R_int / h)

    is long compared with the time the model is actually run for, and long
    compared with the growth time ``h / v_s``. This function returns ``t_eps``
    so the assumption can be checked rather than assumed. ``inf`` is returned
    when the internal reaction vanishes.
    """
    r_int = internal_reaction_rate(params, height, evaporation_flux, law)
    if r_int <= 0.0:
        return math.inf
    rate_density = r_int / height
    return params.crust_porosity * params.crystal_density / rate_density


def growth_timescale(
    params: CrustParameters,
    height: float,
    evaporation_flux: float,
    law: GrowthLaw = GrowthLaw.FINITE_DA,
) -> float:
    """``h / v_s`` [s], the time to double the current crust height."""
    v_s = growth_velocity(params, height, evaporation_flux, law)
    if v_s <= 0.0:
        return math.inf
    return height / v_s


def mean_mass_fraction(
    params: CrustParameters,
    height: float,
    evaporation_flux: float,
    law: GrowthLaw = GrowthLaw.FINITE_DA,
) -> float:
    """Depth-averaged ``C`` on the analytic profile."""
    u = math.sqrt(params.damkohler(height))
    c_star_top = top_supersaturation(params, height, evaporation_flux, law)
    return params.equilibrium_mass_fraction + c_star_top * _tanh_half_over_u(u)


def dissolved_salt_inventory(
    params: CrustParameters,
    height: float,
    evaporation_flux: float,
    law: GrowthLaw = GrowthLaw.FINITE_DA,
) -> float:
    """``M_s,d = int_0^h eps_c S_c rho C dz`` [kg NaCl/m^2]."""
    return (
        params.crust_porosity
        * params.crust_saturation
        * params.rho
        * height
        * mean_mass_fraction(params, height, evaporation_flux, law)
    )


def water_inventory(
    params: CrustParameters,
    height: float,
    evaporation_flux: float,
    law: GrowthLaw = GrowthLaw.FINITE_DA,
) -> float:
    """``M_w,c = int_0^h eps_c S_c rho (1-C) dz`` [kg H2O/m^2]."""
    pore_liquid = (
        params.crust_porosity * params.crust_saturation * params.rho * height
    )
    return pore_liquid - dissolved_salt_inventory(
        params, height, evaporation_flux, law
    )


# --- moving-interface balances and base fluxes --------------------------------


@dataclass(frozen=True)
class InterfaceFluxes:
    """Reconstructed soil/crust base exchange for one surface patch.

    Sign convention: positive is upward, soil -> crust, per unit *projected*
    soil area, as in ``docs/INTERFACE_CONSERVATION.md`` section 1.
    """

    height: float
    evaporation_flux: float          # j   [kg H2O/m^2/s]
    growth_velocity: float           # v_s [m/s]
    top_mass_fraction: float         # C_top  [-]
    top_supersaturation: float       # C*_top [-]
    damkohler: float
    damkohler_top: float
    damkohler_star: float
    peclet: float
    efficiency: float
    internal_reaction: float         # R_int [kg NaCl/m^2/s]
    solution_mass_flux_top: float    # rho U(h) [kg solution/m^2/s]
    solution_mass_flux_base: float   # rho U(0) [kg solution/m^2/s]
    salt_flux_balance_route: float   # F_s,soil from equations (6)-(7) + R_int
    salt_flux_base_route: float      # F_s,soil from rho C U - rho eps D dC/dz
    water_flux: float                # F_w,soil [kg H2O/m^2/s]
    top_skeleton_rate: float         # dM_s,t/dt = (1-eps_c) rho_cr v_s
    new_pore_salt_rate: float        # eps_c rho C_top v_s
    new_pore_water_rate: float       # eps_c rho (1-C_top) v_s
    vapour_sweep_rate: float         # -rho_v v_s
    route_difference: float          # base route minus balance route
    # closed form: -(eps_c S_c rho C*_top + C_eq rho_v) v_s
    route_difference_identity: float

    @property
    def salt_flux(self) -> float:
        """The reported ``F_s,soil``; the moving-interface balance route."""
        return self.salt_flux_balance_route

    @property
    def dissolved_salt_rate(self) -> float:
        """``dM_s,d/dt`` retained by the reduced model: the new pore liquid."""
        return self.new_pore_salt_rate

    @property
    def water_storage_rate(self) -> float:
        """``dM_w,c/dt`` retained by the reduced model: the new pore liquid."""
        return self.new_pore_water_rate


def interface_fluxes(
    params: CrustParameters,
    height: float,
    evaporation_flux: float,
    law: GrowthLaw = GrowthLaw.FINITE_DA,
) -> InterfaceFluxes:
    """Reconstruct ``F_w,soil`` and ``F_s,soil`` demanded by the crust.

    Two independent routes to ``F_s,soil`` are evaluated:

    balance route
        ``F_s = eps_c rho C_top v_s + (1-eps_c) rho_cr v_s + R_int``, i.e. the
        new pore liquid of equation (7), plus the new skeleton, plus the
        depth-integrated internal reaction.
    base-flux route
        ``F_s = rho C_b U(0) - rho eps_c D_s_eff dC/dz|_0+``, the boxed
        expression of ``INTERFACE_CONSERVATION.md`` section 3, with ``rho U(0)``
        reconstructed from equation (6) plus the quasi-steady solution-mass
        balance ``d(rho U)/dz = -a_v rho k_r C*``.

    The two do **not** agree identically, and the gap has a closed form. The
    balance route uses the exact moving-top salt balance, equation (7), whose
    right-hand side carries ``C_top``. The base route inherits equation (8),
    which was obtained by setting ``C ~ C_eq`` at the top. Carrying the algebra
    through gives exactly

        base route - balance route
            = - ( eps_c S_c rho C*_top + C_eq rho_v ) v_s,

    whose leading term is the dissolved salt in the newly created pore volume
    *in excess of equilibrium*, and whose second term is the equation (6)
    vapour sweep that equation (8) drops. ``route_difference`` reports it, and
    ``route_difference_identity`` is that closed form; `verification/` checks
    the two against each other to round-off. The gap is the top-linearisation
    residual of equation (8), not a numerical defect, and it vanishes as
    ``C*_top -> 0``.

    The reported ``F_s,soil`` is the **balance route**, because that is the one
    that matches the tracked crust inventories term by term and therefore keeps
    the coupling conservative.
    """
    if params.crust_saturation < params.connectivity_threshold:
        raise CrustModelDomainExit(
            "crust liquid saturation "
            f"{params.crust_saturation!r} is below the declared connectivity "
            f"threshold {params.connectivity_threshold!r}; the saturated, "
            "capillary-supplied Prat closure is disabled rather than "
            "extrapolated (audit section 5.4)."
        )

    rho = params.rho
    eps = params.crust_porosity
    sat = params.crust_saturation
    da = params.damkohler(height)
    u = math.sqrt(da)
    da_star = params.damkohler_star
    pe = params.peclet(height, evaporation_flux)

    f = (
        efficiency_large_damkohler(da_star)
        if law is GrowthLaw.LARGE_DA
        else efficiency(da, da_star)
    )
    v_s = growth_prefactor(params, evaporation_flux) * f
    c_star_top = top_supersaturation(params, height, evaporation_flux, law)
    c_top = params.equilibrium_mass_fraction + c_star_top
    r_int = internal_reaction_rate(params, height, evaporation_flux, law)

    # Equation (6): total solution/solid/vapour balance at the moving top.
    rho_u_top = (
        rho * eps * sat
        + params.crystal_density * (1.0 - eps)
        - params.vapour_density
    ) * v_s + evaporation_flux
    # Quasi-steady solution-mass balance integrated over the crust.
    rho_u_base = rho_u_top + r_int

    top_skeleton_rate = (1.0 - eps) * params.crystal_density * v_s
    new_pore_salt_rate = eps * sat * rho * c_top * v_s
    new_pore_water_rate = eps * sat * rho * (1.0 - c_top) * v_s

    salt_balance = new_pore_salt_rate + top_skeleton_rate + r_int

    # dC/dz|_0+ = (1/h) dC*/dz*|_0 = C*_top sqrt(Da) / (h sinh(sqrt(Da)))
    base_gradient = c_star_top * _u_over_sinh(u) / height
    salt_base = (
        params.equilibrium_mass_fraction * rho_u_base
        - rho * eps * sat * params.effective_salt_diffusivity * base_gradient
    )

    water_flux = rho_u_base - salt_balance

    return InterfaceFluxes(
        height=height,
        evaporation_flux=evaporation_flux,
        growth_velocity=v_s,
        top_mass_fraction=c_top,
        top_supersaturation=c_star_top,
        damkohler=da,
        damkohler_top=params.damkohler_top(height),
        damkohler_star=da_star,
        peclet=pe,
        efficiency=f,
        internal_reaction=r_int,
        solution_mass_flux_top=rho_u_top,
        solution_mass_flux_base=rho_u_base,
        salt_flux_balance_route=salt_balance,
        salt_flux_base_route=salt_base,
        water_flux=water_flux,
        top_skeleton_rate=top_skeleton_rate,
        new_pore_salt_rate=new_pore_salt_rate,
        new_pore_water_rate=new_pore_water_rate,
        vapour_sweep_rate=-params.vapour_density * v_s,
        route_difference=salt_base - salt_balance,
        route_difference_identity=-(
            eps * sat * rho * c_star_top
            + params.equilibrium_mass_fraction * params.vapour_density
        )
        * v_s,
    )


# --- validity envelope --------------------------------------------------------


def validity_flags(
    params: CrustParameters,
    height: float,
    evaporation_flux: float,
    large_da_threshold: float = 10.0,
    pe_da_ratio_limit: float = 0.1,
) -> dict:
    """Report, without clipping, whether the closure assumptions hold.

    ``large_da_threshold`` is the paper's own ``Da > 10`` statement for using
    ``tanh(sqrt(Da)) ~ 1``. ``pe_da_ratio_limit`` is a *declared* numerical
    reading of the paper's qualitative ``Pe << Da``; the paper gives no number.
    """
    da = params.damkohler(height)
    pe = params.peclet(height, evaporation_flux)
    ratio = abs(pe) / da if da > 0.0 else math.inf
    t_eps = porosity_consumption_timescale(params, height, evaporation_flux)
    t_growth = growth_timescale(params, height, evaporation_flux)
    separation = t_eps / t_growth if t_growth > 0.0 else math.inf
    r_int = internal_reaction_rate(params, height, evaporation_flux)
    top_rate = (
        (1.0 - params.crust_porosity)
        * params.crystal_density
        * growth_velocity(params, height, evaporation_flux)
    )
    return {
        "damkohler": da,
        "damkohler_star": params.damkohler_star,
        "damkohler_top": params.damkohler_top(height),
        "peclet": pe,
        "pe_over_da": ratio,
        "porosity_consumption_timescale_s": t_eps,
        "growth_timescale_s": t_growth,
        "timescale_separation": separation,
        "internal_over_top_precipitation": (
            r_int / top_rate if top_rate > 0.0 else math.inf
        ),
        "large_da_valid": bool(da >= large_da_threshold),
        "pe_much_smaller_than_da": bool(ratio <= pe_da_ratio_limit),
        "constant_porosity_defensible": bool(separation >= 10.0),
        "connected": bool(params.crust_saturation >= params.connectivity_threshold),
        "positive_height": bool(height > 0.0),
        "evaporating": bool(evaporation_flux > 0.0),
    }


# --- explicit activation ------------------------------------------------------


@dataclass(frozen=True)
class SeedEvent:
    """Declared surface precipitation / seed-formation event (P2-A item 6).

    The Prat growth law is a continuum law for an *existing* structure. It is
    not a nucleation model. This class carries the separate, explicitly
    declared rule that creates the structure in the first place.

    ``trigger_mass_fraction`` is a thermodynamic *eligibility* threshold, not a
    validated nucleation time. Crossing it is reported as an eligibility event;
    growth starts only because ``seed_height`` is then declared, by choice.
    """

    trigger_mass_fraction: float
    seed_height: float
    seed_porosity: float
    declared_by: str = "operator"
    note: str = (
        "Solubility crossing is a thermodynamic eligibility condition only. "
        "The seed height and porosity are declared modelling inputs, not "
        "predictions of a nucleation model."
    )

    def __post_init__(self) -> None:
        if not (0.0 < self.trigger_mass_fraction < 1.0):
            raise ValueError("trigger mass fraction must lie in (0, 1)")
        if not math.isfinite(self.seed_height) or self.seed_height <= 0.0:
            raise ValueError("seed height must be positive; h = 0 is not a crust")
        if not (0.0 < self.seed_porosity < 1.0):
            raise ValueError("seed porosity must lie in (0, 1)")

    def eligible(self, surface_mass_fraction: float) -> bool:
        return surface_mass_fraction >= self.trigger_mass_fraction


@dataclass
class CrustState:
    """Mutable inventory state of one surface patch, per projected area."""

    height: float = 0.0
    solid_salt: float = 0.0       # M_s,solid [kg NaCl/m^2]
    dissolved_salt: float = 0.0   # M_s,d     [kg NaCl/m^2]
    water: float = 0.0            # M_w,c     [kg H2O/m^2]
    time: float = 0.0
    active: bool = False
    activation_time: float | None = None
    eligibility_time: float | None = None
    cumulative_salt_from_soil: float = 0.0   # integral of F_s,soil [kg/m^2]
    cumulative_water_from_soil: float = 0.0  # integral of F_w,soil [kg/m^2]
    cumulative_evaporation: float = 0.0      # integral of j        [kg/m^2]
    cumulative_vapour_sweep: float = 0.0     # integral of rho_v v_s [kg/m^2]
    seed_salt: float = 0.0                   # one-off seed transfer [kg/m^2]
    seed_water: float = 0.0                  # one-off seed transfer [kg/m^2]
    history: list = field(default_factory=list)

    def copy(self) -> "CrustState":
        return replace(self, history=list(self.history))

    @property
    def total_salt(self) -> float:
        """``M_s,solid + M_s,d``; the inventory that ``F_s,soil`` must match."""
        return self.solid_salt + self.dissolved_salt


def activate(
    state: CrustState,
    params: CrustParameters,
    seed: SeedEvent,
    evaporation_flux: float,
    time: float,
    law: GrowthLaw = GrowthLaw.FINITE_DA,
) -> CrustState:
    """Create the crust from a declared seed event. A no-op if already active.

    The seed carries a finite solid, dissolved and pore-water inventory. That
    mass is **charged to the soil** as a single instantaneous transfer and
    recorded in ``seed_salt`` / ``seed_water``, so that the component budgets
    still close across the activation event. Creating crust inventory out of
    nothing would reproduce exactly the defect P2-A item 8 forbids: solid salt
    accumulating in the crust with no corresponding soil salt flux. The seed
    *size* remains a declared modelling input, not a prediction.
    """
    if state.active:
        return state
    if abs(seed.seed_porosity - params.crust_porosity) > 1.0e-12:
        raise ValueError(
            "seed porosity must equal the declared crust porosity; a different "
            "seed porosity would require a second, separately justified "
            "porosity state"
        )
    state.height = seed.seed_height
    state.active = True
    state.activation_time = time
    state.solid_salt = (
        (1.0 - params.crust_porosity) * params.crystal_density * seed.seed_height
    )
    state.dissolved_salt = dissolved_salt_inventory(
        params, state.height, evaporation_flux, law
    )
    state.water = water_inventory(params, state.height, evaporation_flux, law)
    state.seed_salt = state.solid_salt + state.dissolved_salt
    state.seed_water = state.water
    state.cumulative_salt_from_soil += state.seed_salt
    state.cumulative_water_from_soil += state.seed_water
    return state


# --- supply limitation (diagnostic only) --------------------------------------


@dataclass(frozen=True)
class SupplyDiagnostic:
    """Diagnostic comparison of crust demand against an available soil flux."""

    salt_demand: float
    salt_available: float
    water_demand: float
    water_available: float
    salt_ratio: float
    water_ratio: float
    salt_limited: bool
    water_limited: bool
    note: str = (
        "Upper-bound diagnostic only. Capping v_s by the available salt flux "
        "is NOT conservation-consistent: it would change v_s without updating "
        "C_b, C_top, U, the dissolved storage, the new pore water and the "
        "internal reaction (INTERFACE_CONSERVATION.md section 7.1)."
    )


def supply_limitation_diagnostic(
    fluxes: InterfaceFluxes,
    available_salt_flux: float,
    available_water_flux: float,
    raise_on_shortfall: bool = False,
) -> SupplyDiagnostic:
    """Report whether the soil can meet the reconstructed crust demand.

    Never clips. With ``raise_on_shortfall=True`` it raises
    ``CrustModelDomainExit``, which is the conservative behaviour: the closure
    is outside its validity envelope and must not be silently continued.
    """
    salt_demand = fluxes.salt_flux
    water_demand = fluxes.water_flux
    salt_ratio = (
        salt_demand / available_salt_flux
        if available_salt_flux != 0.0
        else math.inf
    )
    water_ratio = (
        water_demand / available_water_flux
        if available_water_flux != 0.0
        else math.inf
    )
    salt_limited = salt_demand > available_salt_flux
    water_limited = water_demand > available_water_flux
    if raise_on_shortfall and (salt_limited or water_limited):
        raise CrustModelDomainExit(
            "soil supply cannot meet the reconstructed crust demand "
            f"(salt ratio {salt_ratio!r}, water ratio {water_ratio!r}); the "
            "saturated quasi-steady Prat closure is invalid under undersupply "
            "and is stopped rather than clipped."
        )
    return SupplyDiagnostic(
        salt_demand=salt_demand,
        salt_available=available_salt_flux,
        water_demand=water_demand,
        water_available=available_water_flux,
        salt_ratio=salt_ratio,
        water_ratio=water_ratio,
        salt_limited=salt_limited,
        water_limited=water_limited,
    )
