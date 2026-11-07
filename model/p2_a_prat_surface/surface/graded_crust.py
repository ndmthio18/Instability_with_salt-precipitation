"""P2-C: graded (top-refined) depth-resolved crust with evolving porosity.

Why this module exists
----------------------
`evolving_crust.py` (P2-B) resolves Prat equation (2) on a **uniform** grid of
fixed cell size. That is adequate at the P2-A declared parameters, where the
reaction-diffusion length

    delta_r = sqrt( eps D_s,eff / ((1-C_eq) a_v k_r) ) = h / sqrt(Da)

is 227 micrometres, i.e. *larger* than the 100-micrometre seed, so the whole
crust reacts almost uniformly. Transferring the `salt_grow` parameters
(``eps_0 = 0.43``, ``k_r = 2.3e-3 m/s``, ``a_v = 36(1-eps_0)/d_sub``) moves
``delta_r`` to about 0.54 micrometres. A uniform grid that resolves it over a
0.1-1 mm crust needs 10^4-10^5 cells and is not affordable in pure Python, and
an unresolved run would be meaningless. This module therefore generalises the
P2-B geometry from a single ``cell_size`` to an explicit list of cell
thicknesses, so the grid can be fine at the moving top and coarse in the inert
wake.

Nothing in `prat_crust.py` or `evolving_crust.py` is modified. P2-B stays the
uniform-grid reference that this module must reproduce cell for cell, which
`verification/test_p2c_graded.py` asserts.

What is the same as P2-B
------------------------
Everything physical. ``z`` is still a material coordinate (the crust base is
fixed and new material appears only at the moving top), so equation (2) still
applies pointwise with no Lagrangian mapping; the quasi-steady, ``Pe << Da``,
linearised-sink reduced transport problem is unchanged; the top closure is
still exact in ``C_top`` and still reduces to equation (28); the tracked state
per cell is still the solid volume per unit projected area, so no porosity
update can leak mass.

What is new, and why
--------------------
1. **Non-uniform cells.** ``thicknesses[i]`` replaces a single ``cell_size``.
   New material is created at the declared ``top_cell_size``.

2. **Wake coarsening.** Cells deeper than ``freeze_depth`` below the moving
   front may be merged pairwise, conserving thickness and solid volume
   exactly. Beyond a few ``delta_r`` the supersaturation is smaller than
   ``exp(-freeze_depth/delta_r)`` and the porosity profile is flat, so this is
   a declared numerical device whose only cost is resolution of the wake
   porosity. It is off by default and its influence is measured, never
   assumed.

3. **Explicit `salt_grow` closures.** ``a_v = 36(1-eps)/d_sub`` and
   ``D_eff = sqrt(eps) D_molecular`` are the B2 closures of the reference
   project. They are *not* expressible as the P2-B power laws
   ``(eps/eps_0)**m`` -- in particular ``a_v`` **increases** as porosity falls
   -- so they are provided as named closures and selected explicitly.

4. **Two corrections to the P2-B base-flux diagnostics.** These change no
   governing equation and no growth velocity; they change only the
   independently reconstructed soil-side fluxes:

   a. *Pore-volume release.* Precipitation removes salt mass ``R`` from the
      solution but frees pore volume ``R/rho_cr`` which refills with solution,
      so the quasi-steady solution-mass balance is

          d(rho U)/dz = -R (1 - rho S / rho_cr) ,

      not ``-R``. P2-B used the latter, which overstates ``rho U(0)`` by
      ``R_int rho S / rho_cr`` and is the main part of the 7.5-10.4 % route
      discrepancy reported in the P2-B verification.

   b. *Birth porosity at the moving interface.* The new pore volume created at
      the top is filled at the **deposition** porosity ``eps_0``, not at the
      current top-cell porosity, so

          rho U(h) = (rho eps_0 S + rho_cr (1-eps_0) - rho_v) v_s + j .

   Both legacy forms are still evaluated and reported, so the correction is
   visible rather than silently applied.

5. **Connectivity is the primary stop.** P2-B stopped at the geometric
   clogging porosity. `INTERFACE_CONSERVATION.md` and the audit treat loss of
   hydraulic connectivity as the point beyond which the saturated capillary
   supply closure is void, and it is reached first. Both are recorded.

Nothing here caps, clips or normalises the growth velocity.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from .prat_crust import CrustModelDomainExit
from .thermo import BrineDensityLaw

AREA_CLOSURES = ("power", "salt_grow_b2")
DIFFUSIVITY_CLOSURES = ("power", "sqrt_eps")


@dataclass(frozen=True)
class GradedCrustParameters:
    """Crust parameters, a superset of `evolving_crust.EvolvingCrustParameters`.

    With the default closures (``power`` with exponent zero) every method here
    returns exactly what the P2-B parameter object returns.
    """

    equilibrium_mass_fraction: float
    deposition_porosity: float
    crystal_density: float
    specific_area: float
    reaction_coefficient: float
    effective_salt_diffusivity: float
    beta: float = 2.0
    solution_density: float | None = None
    vapour_density: float = 0.0
    crust_saturation: float = 1.0
    clogging_porosity: float = 0.01
    connectivity_porosity: float = 0.02
    diffusivity_exponent: float = 0.0
    specific_area_exponent: float = 0.0
    linearised_top_closure: bool = False
    #: Named closures. ``power`` reproduces P2-B exactly.
    area_closure: str = "power"
    diffusivity_closure: str = "power"
    #: Required by ``salt_grow_b2``: the substrate/crystal size d_sub [m].
    substrate_size: float | None = None
    #: Required by ``sqrt_eps``: the molecular salt diffusivity [m^2/s].
    molecular_diffusivity: float | None = None
    #: Legacy P2-B base-flux closures, for measuring the correction.
    legacy_liquid_closure: bool = False

    def __post_init__(self) -> None:
        if not (0.0 < self.equilibrium_mass_fraction < 1.0):
            raise ValueError("C_eq must lie strictly in (0, 1)")
        if not (0.0 < self.deposition_porosity < 1.0):
            raise ValueError("deposition porosity must lie strictly in (0, 1)")
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
        if not (0.0 <= self.clogging_porosity < self.deposition_porosity):
            raise ValueError("clogging porosity must lie in [0, eps_0)")
        if self.area_closure not in AREA_CLOSURES:
            raise ValueError(f"unknown area closure {self.area_closure!r}")
        if self.diffusivity_closure not in DIFFUSIVITY_CLOSURES:
            raise ValueError(
                f"unknown diffusivity closure {self.diffusivity_closure!r}"
            )
        if self.area_closure == "salt_grow_b2" and not self.substrate_size:
            raise ValueError("the salt_grow_b2 area closure needs substrate_size")
        if self.diffusivity_closure == "sqrt_eps" and not self.molecular_diffusivity:
            raise ValueError("the sqrt_eps closure needs molecular_diffusivity")

    @property
    def rho(self) -> float:
        if self.solution_density is not None:
            return self.solution_density
        return BrineDensityLaw().density(self.equilibrium_mass_fraction)

    # -- declared porosity-dependent closures ------------------------------

    def diffusivity(self, porosity: float) -> float:
        eps = max(porosity, 0.0)
        if self.diffusivity_closure == "sqrt_eps":
            return math.sqrt(eps) * self.molecular_diffusivity
        if self.diffusivity_exponent == 0.0:
            return self.effective_salt_diffusivity
        ratio = eps / self.deposition_porosity
        return self.effective_salt_diffusivity * ratio**self.diffusivity_exponent

    def area(self, porosity: float) -> float:
        eps = max(porosity, 0.0)
        if self.area_closure == "salt_grow_b2":
            return 36.0 * (1.0 - eps) / self.substrate_size
        if self.specific_area_exponent == 0.0:
            return self.specific_area
        ratio = eps / self.deposition_porosity
        return self.specific_area * ratio**self.specific_area_exponent

    # -- declared a-priori scales ------------------------------------------

    def reaction_length(self, porosity: float | None = None) -> float:
        """``delta_r = sqrt(eps D / ((1-C_eq) a_v k_r))`` at a given porosity."""
        eps = self.deposition_porosity if porosity is None else porosity
        one = 1.0 - self.equilibrium_mass_fraction
        return math.sqrt(
            eps
            * self.diffusivity(eps)
            / (one * self.area(eps) * self.reaction_coefficient)
        )

    def damkohler_star(self, porosity: float | None = None) -> float:
        """``Da* = beta^2 (1-eps)^2 (1-C_eq) k_r / (eps D a_v)``."""
        eps = self.deposition_porosity if porosity is None else porosity
        one = 1.0 - self.equilibrium_mass_fraction
        return (
            self.beta**2
            * (1.0 - eps) ** 2
            * one
            * self.reaction_coefficient
            / (eps * self.diffusivity(eps) * self.area(eps))
        )

    def predicted_wake_porosity_drop(self) -> float:
        """``(1-eps_0)/sqrt(Da*)``: the thick-crust (``Da >> 1``) prediction.

        Derivation (registered before measurement, see
        ``config/acceptance_criteria_p2c.json``): a material point deposited at
        the top at ``eps_0`` sees ``C*(d) = C*_top exp(-d/delta_r)`` as it is
        buried at speed ``v_s``, so

            d eps = -(a_v rho k_r C*_top / rho_cr) exp(-v_s t/delta_r) dt

        integrates to ``-a_v delta_r / beta == -(1-eps_0)/sqrt(Da*)``, which is
        independent of the evaporation flux, of ``C*_top`` and of ``rho``.
        """
        return (1.0 - self.deposition_porosity) / math.sqrt(self.damkohler_star())

    def predicted_sealing_growth(self, stop_porosity: float | None = None) -> float:
        """``int beta/a_v(eps) d eps``: the thin-crust (``Da << 1``) law."""
        stop = (
            self.connectivity_porosity if stop_porosity is None else stop_porosity
        )
        eps_0 = self.deposition_porosity
        if self.area_closure == "salt_grow_b2":
            return (
                self.beta
                * self.substrate_size
                / 36.0
                * math.log((1.0 - stop) / (1.0 - eps_0))
            )
        m = self.specific_area_exponent
        pre = self.beta * eps_0**m / self.specific_area
        if m == 1.0:
            return pre * math.log(eps_0 / stop)
        return pre * (eps_0 ** (1.0 - m) - stop ** (1.0 - m)) / (1.0 - m)


def _thomas(lower: list, diag: list, upper: list, rhs: list) -> list:
    """Solve a tridiagonal system (Thomas algorithm, no scipy)."""
    n = len(diag)
    c = [0.0] * n
    d = [0.0] * n
    beta = diag[0]
    if beta == 0.0:
        raise CrustModelDomainExit("singular crust transport system")
    c[0] = upper[0] / beta
    d[0] = rhs[0] / beta
    for i in range(1, n):
        beta = diag[i] - lower[i] * c[i - 1]
        if beta == 0.0:
            raise CrustModelDomainExit("singular crust transport system")
        c[i] = upper[i] / beta if i < n - 1 else 0.0
        d[i] = (rhs[i] - lower[i] * d[i - 1]) / beta
    x = [0.0] * n
    x[-1] = d[-1]
    for i in range(n - 2, -1, -1):
        x[i] = d[i] - c[i] * x[i + 1]
    return x


# --- grid construction --------------------------------------------------------


def graded_seed_grid(
    height: float,
    top_cell: float,
    window: float = 0.0,
    ratio: float = 1.0,
    maximum_cell: float | None = None,
) -> list:
    """Cell thicknesses, bottom to top, tiling ``[0, height]`` exactly.

    The top ``window`` metres are uniform at ``top_cell``; below that the cells
    grow geometrically by ``ratio`` up to ``maximum_cell``. ``ratio = 1`` and
    ``window >= height`` give a uniform grid, which is the P2-B geometry.
    """
    if height <= 0.0:
        raise ValueError("seed height must be positive")
    if top_cell <= 0.0:
        raise ValueError("top cell size must be positive")
    if ratio < 1.0:
        raise ValueError("grading ratio must be at least 1")
    top_cell = min(top_cell, height)
    cap = maximum_cell if maximum_cell else height
    from_top = []
    covered = 0.0
    size = top_cell
    while covered < height - 1.0e-15 * height:
        if covered >= window:
            size = min(size * ratio, cap)
        take = min(size, height - covered)
        from_top.append(take)
        covered += take
        if len(from_top) > 10_000_000:
            raise CrustModelDomainExit("graded grid exceeded the cell cap")
    # Absorb a sliver bottom cell into its neighbour so no cell is degenerate.
    # The pop must happen first: ``from_top[-2] += from_top.pop()`` resolves
    # the target index after the pop and would merge into the wrong cell.
    if len(from_top) > 1 and from_top[-1] < 0.5 * from_top[-2]:
        sliver = from_top.pop()
        from_top[-1] += sliver
    return list(reversed(from_top))


# --- state --------------------------------------------------------------------


@dataclass
class GradedCrustState:
    """Depth-resolved crust on an explicit, possibly non-uniform, grid.

    ``solid[i]`` is the solid volume per unit projected area in cell ``i``
    (metres of solid) and ``thicknesses[i]`` its occupied thickness, so
    ``height == sum(thicknesses)`` by construction.
    """

    top_cell_size: float
    thicknesses: list = field(default_factory=list)
    solid: list = field(default_factory=list)
    height: float = 0.0
    time: float = 0.0
    active: bool = False
    activation_time: float | None = None
    cumulative_salt_from_soil: float = 0.0
    cumulative_water_from_soil: float = 0.0
    cumulative_evaporation: float = 0.0
    cumulative_vapour_sweep: float = 0.0
    cumulative_internal_reaction: float = 0.0
    cumulative_top_skeleton: float = 0.0
    seed_salt: float = 0.0
    seed_water: float = 0.0
    sealed: bool = False
    sealed_time: float | None = None
    disconnected: bool = False
    disconnected_time: float | None = None
    stopped: bool = False
    stop_reason: str | None = None
    merges: int = 0
    solid_salt_mass: float = 0.0
    dissolved_salt: float = 0.0
    water: float = 0.0
    history: list = field(default_factory=list)

    # -- geometry -----------------------------------------------------------

    def centres(self) -> list:
        out = []
        lower = 0.0
        for t in self.thicknesses:
            out.append(lower + 0.5 * t)
            lower += t
        return out

    def porosity(self, index: int) -> float:
        t = self.thicknesses[index]
        if t <= 0.0:
            return 0.0
        return 1.0 - self.solid[index] / t

    def porosity_profile(self) -> list:
        return [self.porosity(i) for i in range(len(self.thicknesses))]

    @property
    def minimum_porosity(self) -> float:
        profile = self.porosity_profile()
        return min(profile) if profile else 1.0

    @property
    def solid_salt(self) -> float:
        """Exact sum of the tracked solid state [m of solid per m^2]."""
        return sum(self.solid)

    def copy(self) -> "GradedCrustState":
        clone = GradedCrustState(top_cell_size=self.top_cell_size)
        for key, value in self.__dict__.items():
            setattr(clone, key, list(value) if isinstance(value, list) else value)
        return clone


@dataclass(frozen=True)
class GradedProfileSolution:
    """One quasi-steady solve of the variable-porosity, variable-grid crust."""

    height: float
    evaporation_flux: float
    top_supersaturation: float
    top_mass_fraction: float
    growth_velocity: float
    top_gradient_coefficient: float
    profile: list
    centres: list
    porosity: list
    thicknesses: list
    internal_reaction: float
    top_skeleton_rate: float
    dissolved_salt: float
    water: float
    solution_mass_flux_top: float
    solution_mass_flux_base: float
    solution_mass_flux_base_legacy: float
    base_gradient: float
    salt_flux_base_route: float
    salt_flux_base_route_legacy: float
    minimum_porosity: float
    top_porosity: float
    reaction_length: float
    cells_across_reaction_length: float
    reaction_depth_95: float
    #: Closed form of (transport route - inventory route); see `solve_profile`.
    #: Derived for a halite-equilibrated base and reported only there.
    route_difference_identity: float
    solution_mass_flux_top_relative: float
    #: The NaCl mass fraction actually imposed at z = 0.
    base_mass_fraction: float
    #: True when the top closure has no positive supersaturation root.
    growth_arrested: bool


def solve_profile(
    state: GradedCrustState,
    params: GradedCrustParameters,
    evaporation_flux: float,
    base_mass_fraction: float | None = None,
) -> GradedProfileSolution:
    """Quasi-steady crust transport and the exact top closure.

    ``base_mass_fraction`` is the NaCl mass fraction imposed at the crust base
    ``z = 0``. ``None`` means ``C_eq``, i.e. Prat's own halite-equilibrated
    base, and then every line below is the committed P2-C computation
    unchanged. A different value is needed only by the P2-D soil/crust
    coupling, where the interface concentration is an unknown to be solved
    for rather than a prescribed equilibrium.

    With ``u = C - C_eq`` the reduced problem is linear with ``u(0) = u_b``
    and ``u(h) = s``, so it separates as ``u = u_b psi + s phi`` where
    ``phi(0)=0, phi(h)=1`` and ``psi(0)=1, psi(h)=0``. ``psi`` is solved only
    when ``u_b != 0``, so the default path is bit-identical to P2-C.
    """
    if not state.active or state.height <= 0.0:
        raise CrustModelDomainExit(
            "the crust transport problem requires an active crust with h > 0"
        )
    n = len(state.thicknesses)
    if n == 0:
        raise CrustModelDomainExit("no occupied crust cells")

    rho = params.rho
    sat = params.crust_saturation
    c_eq = params.equilibrium_mass_fraction
    one = 1.0 - c_eq

    t_cells = list(state.thicknesses)
    z_cells = state.centres()
    eps_cells = [state.porosity(i) for i in range(n)]
    cond = [eps_cells[k] * params.diffusivity(eps_cells[k]) for k in range(n)]
    lam_cells = [
        one * params.area(eps_cells[k]) * params.reaction_coefficient
        for k in range(n)
    ]

    lower = [0.0] * n
    diag = [0.0] * n
    upper = [0.0] * n
    rhs = [0.0] * n
    for k in range(n):
        if k == 0:
            distance = z_cells[0]
            coefficient = cond[0] / distance if distance > 0.0 else 0.0
            diag[k] -= coefficient              # phi(0) = 0
        else:
            distance = z_cells[k] - z_cells[k - 1]
            face = 0.5 * (cond[k] + cond[k - 1])
            coefficient = face / distance
            diag[k] -= coefficient
            lower[k] += coefficient
        if k == n - 1:
            distance = state.height - z_cells[k]
            coefficient = cond[k] / distance if distance > 0.0 else 0.0
            diag[k] -= coefficient
            rhs[k] -= coefficient               # phi(h) = 1
        else:
            distance = z_cells[k + 1] - z_cells[k]
            face = 0.5 * (cond[k] + cond[k + 1])
            coefficient = face / distance
            diag[k] -= coefficient
            upper[k] += coefficient
        diag[k] -= lam_cells[k] * t_cells[k]

    phi = _thomas(lower, diag, upper, rhs)

    c_eq_base = c_eq if base_mass_fraction is None else base_mass_fraction
    u_base = c_eq_base - c_eq
    if u_base != 0.0:
        rhs_psi = [0.0] * n
        # phi(0)=0 contributes nothing; psi(0)=1 contributes at the bottom
        # face only. The matrix is the same, so only the right-hand side
        # changes: move the Dirichlet value to the right-hand side.
        distance0 = z_cells[0]
        coefficient0 = cond[0] / distance0 if distance0 > 0.0 else 0.0
        rhs_psi[0] -= coefficient0
        psi = _thomas(lower, diag, upper, rhs_psi)
    else:
        psi = [0.0] * n

    top_distance = state.height - z_cells[-1]
    if top_distance <= 0.0:
        raise CrustModelDomainExit("degenerate top cell in the crust grid")
    top_gradient = (1.0 - phi[-1]) / top_distance
    eps_top = eps_cells[-1]
    a_coefficient = rho * cond[-1] * top_gradient
    b_coefficient = rho * cond[-1] * (0.0 - psi[-1]) / top_distance

    k_kinetic = params.beta * params.reaction_coefficient * rho / params.crystal_density
    skeleton = params.crystal_density * (1.0 - params.deposition_porosity)

    def residual(s: float) -> float:
        v_s = k_kinetic * s
        c_top = c_eq if params.linearised_top_closure else c_eq + s
        return (
            a_coefficient * s
            + b_coefficient * u_base
            - c_top * evaporation_flux
            + skeleton * (1.0 - c_top) * v_s
            + c_top * params.vapour_density * v_s
        )

    # residual(0) >= 0 means the top is not supersaturated: the Prat growth
    # closure has no positive root and growth is arrested. Reported, not
    # silently returned as zero velocity.
    growth_arrested = residual(0.0) >= 0.0
    s_top = _solve_top_supersaturation(residual, c_eq, evaporation_flux)
    v_s = k_kinetic * s_top
    c_top = c_eq + s_top

    c_star = [s_top * phi[k] + u_base * psi[k] for k in range(n)]
    internal = 0.0
    dissolved = 0.0
    water = 0.0
    cell_reaction = [0.0] * n
    for k in range(n):
        cell_reaction[k] = (
            params.area(eps_cells[k])
            * rho
            * params.reaction_coefficient
            * c_star[k]
            * t_cells[k]
        )
        internal += cell_reaction[k]
        pore = eps_cells[k] * sat * rho * t_cells[k]
        dissolved += pore * (c_eq + c_star[k])
        water += pore * (1.0 - c_eq - c_star[k])

    # Moving-interface solution flux. The new pore volume is created at the
    # DEPOSITION porosity, so eps_0 appears here, not the top-cell porosity.
    eps_interface = eps_top if params.legacy_liquid_closure else params.deposition_porosity
    rho_u_top = (
        rho * eps_interface * sat + skeleton - params.vapour_density
    ) * v_s + evaporation_flux
    # Pore-volume release: d(rho U)/dz = -R (1 - rho S / rho_cr).
    release = 1.0 - rho * sat / params.crystal_density
    rho_u_base = rho_u_top + internal * release
    rho_u_base_legacy = rho_u_top + internal

    base_gradient = (c_star[0] - u_base) / z_cells[0] if z_cells[0] > 0.0 else 0.0
    diffusive_base = (
        -rho * eps_cells[0] * sat * params.diffusivity(eps_cells[0]) * base_gradient
    )
    salt_base = c_eq_base * rho_u_base + diffusive_base
    salt_base_legacy = c_eq_base * rho_u_base_legacy + diffusive_base

    # Closed form of the remaining route difference.
    #
    # Summing the discrete cell balances of the reduced problem gives
    # ``A s + D_base = (1-C_eq) R_int`` exactly, and substituting that, the top
    # closure and rho U(0) into ``C_eq rho U(0) + D_base`` and subtracting the
    # exact inventory identity ``d(M_solid + M_d)/dt`` leaves
    #
    #     F_s,transport - F_s,inventory = - s L_top ,
    #     L_top = j + rho_cr (1-eps_0) v_s - rho_v v_s ,
    #
    # i.e. the solution flux at the top measured relative to the moving
    # interface, times the surface supersaturation. This is the residual of
    # the LINEARISED SINK ``(1-C_eq)`` that Prat's own reduction uses in place
    # of ``(1-C)`` (audit section 3), not a bookkeeping error: it vanishes as
    # ``s -> 0`` and is independent of the grid and of the time step. It is
    # reported so the quasi-steady reduction can be audited rather than hidden.
    route_identity = -s_top * (
        evaporation_flux + skeleton * v_s - params.vapour_density * v_s
    )

    delta_r = params.reaction_length(eps_top)
    across = 0.0
    depth = 0.0
    for k in range(n - 1, -1, -1):
        if depth >= delta_r:
            break
        across += min(t_cells[k], delta_r - depth) / t_cells[k]
        depth += t_cells[k]
    # Depth below the front containing 95 % of the instantaneous reaction.
    target = 0.95 * internal
    running = 0.0
    depth_95 = 0.0
    for k in range(n - 1, -1, -1):
        if running + cell_reaction[k] >= target:
            short = (target - running) / cell_reaction[k] if cell_reaction[k] > 0 else 0.0
            depth_95 += short * t_cells[k]
            break
        running += cell_reaction[k]
        depth_95 += t_cells[k]

    return GradedProfileSolution(
        height=state.height,
        evaporation_flux=evaporation_flux,
        top_supersaturation=s_top,
        top_mass_fraction=c_top,
        growth_velocity=v_s,
        top_gradient_coefficient=a_coefficient,
        profile=c_star,
        centres=z_cells,
        porosity=eps_cells,
        thicknesses=t_cells,
        internal_reaction=internal,
        top_skeleton_rate=skeleton * v_s,
        dissolved_salt=dissolved,
        water=water,
        solution_mass_flux_top=rho_u_top,
        solution_mass_flux_base=rho_u_base,
        solution_mass_flux_base_legacy=rho_u_base_legacy,
        base_gradient=base_gradient,
        salt_flux_base_route=salt_base,
        salt_flux_base_route_legacy=salt_base_legacy,
        minimum_porosity=min(eps_cells),
        top_porosity=eps_top,
        reaction_length=delta_r,
        cells_across_reaction_length=across,
        reaction_depth_95=depth_95,
        route_difference_identity=route_identity,
        base_mass_fraction=c_eq_base,
        growth_arrested=growth_arrested,
        solution_mass_flux_top_relative=(
            evaporation_flux + skeleton * v_s - params.vapour_density * v_s
        ),
    )


def _solve_top_supersaturation(residual, c_eq: float, flux: float) -> float:
    """Bracketed bisection for the scalar top closure."""
    if flux <= 0.0:
        return 0.0
    low = 0.0
    if residual(low) >= 0.0:
        return low
    high = max(1.0e-12, 1.0e-3 * c_eq)
    for _ in range(200):
        if residual(high) >= 0.0:
            break
        high *= 2.0
        if c_eq + high >= 1.0:
            raise CrustModelDomainExit(
                "the crust-top closure has no root below C_top = 1"
            )
    else:
        raise CrustModelDomainExit("failed to bracket the crust-top closure")
    for _ in range(200):
        mid = 0.5 * (low + high)
        if residual(mid) < 0.0:
            low = mid
        else:
            high = mid
        if high - low <= 1.0e-16 * max(high, 1.0e-30):
            break
    return 0.5 * (low + high)


# --- activation ---------------------------------------------------------------


def activate(
    state: GradedCrustState,
    params: GradedCrustParameters,
    seed_height: float,
    evaporation_flux: float,
    time: float,
    window: float = 0.0,
    ratio: float = 1.0,
    maximum_cell: float | None = None,
    base_mass_fraction: float | None = None,
) -> GradedCrustState:
    """Create the seed crust at the declared deposition porosity.

    As in P2-A and P2-B the seed inventory is charged to the soil so that the
    component budgets close across the activation event.
    """
    if state.active:
        return state
    if seed_height <= 0.0:
        raise ValueError("seed height must be positive; h = 0 is not a crust")
    state.thicknesses = graded_seed_grid(
        seed_height, state.top_cell_size, window, ratio, maximum_cell
    )
    state.solid = [(1.0 - params.deposition_porosity) * t for t in state.thicknesses]
    state.height = sum(state.thicknesses)
    state.active = True
    state.activation_time = time
    solution = solve_profile(state, params, evaporation_flux, base_mass_fraction)
    state.solid_salt_mass = params.crystal_density * state.solid_salt
    state.dissolved_salt = solution.dissolved_salt
    state.water = solution.water
    state.seed_salt = state.solid_salt_mass + state.dissolved_salt
    state.seed_water = state.water
    state.cumulative_salt_from_soil += state.seed_salt
    state.cumulative_water_from_soil += state.seed_water
    return state


# --- wake coarsening ----------------------------------------------------------


def coarsen(
    state: GradedCrustState,
    freeze_depth: float,
    maximum_cell: float,
) -> int:
    """Merge adjacent cell pairs deeper than ``freeze_depth`` below the front.

    Thickness and solid volume are both additive, so a merge is exact in mass
    and exact in porosity wherever the merged cells already share a porosity.
    The number of merges is recorded on the state.
    """
    if freeze_depth <= 0.0 or maximum_cell <= 0.0:
        return 0
    limit = state.height - freeze_depth
    if limit <= 0.0:
        return 0
    thicknesses = []
    solid = []
    merged = 0
    lower = 0.0
    index = 0
    n = len(state.thicknesses)
    while index < n:
        t = state.thicknesses[index]
        s = state.solid[index]
        upper = lower + t
        # Only merge a pair entirely inside the frozen wake.
        while (
            index + 1 < n
            and upper + state.thicknesses[index + 1] <= limit
            and t + state.thicknesses[index + 1] <= maximum_cell
        ):
            index += 1
            t += state.thicknesses[index]
            s += state.solid[index]
            upper = lower + t
            merged += 1
        thicknesses.append(t)
        solid.append(s)
        lower = upper
        index += 1
    if merged:
        state.thicknesses = thicknesses
        state.solid = solid
        state.merges += merged
    return merged


# --- time integration ---------------------------------------------------------


def step(
    state: GradedCrustState,
    params: GradedCrustParameters,
    evaporation_flux: float,
    dt: float,
    freeze_depth: float = 0.0,
    maximum_cell: float = 0.0,
    maximum_cells: int = 2_000_000,
    base_mass_fraction: float | None = None,
) -> dict:
    """Advance porosity, the moving top and the inventories by ``dt``.

    The soil-side fluxes are the inventory (balance) route, evaluated from the
    actual change of the tracked states, so

        F_s,soil = d(M_s,solid + M_s,d)/dt
        F_w,soil = j - rho_v v_s + dM_w,c/dt

    hold exactly. The independent transport-route reconstruction is evaluated
    alongside and reported, in both the corrected and the legacy P2-B form.
    """
    if dt <= 0.0 or not math.isfinite(dt):
        raise ValueError("time step must be positive and finite")
    if not state.active:
        raise CrustModelDomainExit("cannot step an inactive crust")
    if state.stopped:
        raise CrustModelDomainExit(
            f"crust stopped at t = {state.stop_reason!r}; the saturated "
            "capillary-supply closure is disabled rather than extrapolated"
        )

    before = solve_profile(state, params, evaporation_flux, base_mass_fraction)
    solid_before = params.crystal_density * state.solid_salt
    dissolved_before = before.dissolved_salt
    water_before = before.water
    rho = params.rho

    # 1. internal reaction, cell by cell (equation 2, pointwise in z)
    internal_volume = 0.0
    clamped = False
    for k in range(len(state.thicknesses)):
        rate = (
            params.area(before.porosity[k])
            * rho
            * params.reaction_coefficient
            * before.profile[k]
        )
        added = rate * before.thicknesses[k] * dt / params.crystal_density
        available = before.thicknesses[k] - state.solid[k]
        if added > available:
            added = max(0.0, available)
            clamped = True
        state.solid[k] += added
        internal_volume += added

    # 2. moving top: new material at the declared deposition porosity
    growth = before.growth_velocity * dt
    remaining = growth
    skeleton_fraction = 1.0 - params.deposition_porosity
    # Fill the current top cell up to top_cell_size, then open new cells.
    room = state.top_cell_size - state.thicknesses[-1]
    if room > 0.0 and remaining > 0.0:
        take = min(room, remaining)
        state.thicknesses[-1] += take
        state.solid[-1] += skeleton_fraction * take
        remaining -= take
    while remaining > 0.0:
        if len(state.thicknesses) >= maximum_cells:
            raise CrustModelDomainExit("crust grid exceeded the declared cell cap")
        take = min(state.top_cell_size, remaining)
        state.thicknesses.append(take)
        state.solid.append(skeleton_fraction * take)
        remaining -= take
    state.height += growth
    state.time += dt

    if freeze_depth > 0.0 and maximum_cell > 0.0:
        coarsen(state, freeze_depth, maximum_cell)

    after = solve_profile(state, params, evaporation_flux, base_mass_fraction)
    solid_after = params.crystal_density * state.solid_salt
    state.solid_salt_mass = solid_after
    state.dissolved_salt = after.dissolved_salt
    state.water = after.water

    salt_flux = (
        (solid_after - solid_before) + (after.dissolved_salt - dissolved_before)
    ) / dt
    water_storage_rate = (after.water - water_before) / dt
    water_flux = (
        evaporation_flux
        - params.vapour_density * before.growth_velocity
        + water_storage_rate
    )

    state.cumulative_salt_from_soil += salt_flux * dt
    state.cumulative_water_from_soil += water_flux * dt
    state.cumulative_evaporation += evaporation_flux * dt
    state.cumulative_vapour_sweep += params.vapour_density * before.growth_velocity * dt
    state.cumulative_internal_reaction += internal_volume * params.crystal_density
    state.cumulative_top_skeleton += (
        skeleton_fraction * params.crystal_density * growth
    )

    minimum = after.minimum_porosity
    if minimum <= params.connectivity_porosity and not state.disconnected:
        state.disconnected = True
        state.disconnected_time = state.time
        state.stopped = True
        state.stop_reason = "hydraulic_connectivity_loss"
    if minimum <= params.clogging_porosity and not state.sealed:
        state.sealed = True
        state.sealed_time = state.time
        state.stopped = True
        state.stop_reason = "geometric_clogging"

    scale = max(abs(salt_flux), 1.0e-300)
    record = {
        "time_s": state.time,
        "dt_s": dt,
        "crust_height_m": state.height,
        "cells": len(state.thicknesses),
        "merges": state.merges,
        "evaporation_flux_kg_m2_s": evaporation_flux,
        "growth_velocity_m_s": before.growth_velocity,
        "top_supersaturation": before.top_supersaturation,
        "top_mass_fraction": before.top_mass_fraction,
        "top_porosity": after.top_porosity,
        "minimum_porosity": minimum,
        "base_porosity": after.porosity[0],
        "deposition_porosity": params.deposition_porosity,
        "mean_porosity": sum(
            p * t for p, t in zip(after.porosity, after.thicknesses)
        )
        / state.height,
        "reaction_length_m": before.reaction_length,
        "cells_across_reaction_length": before.cells_across_reaction_length,
        "reaction_depth_95_m": before.reaction_depth_95,
        "internal_reaction_kg_m2_s": before.internal_reaction,
        "internal_reaction_applied_kg_m2_s": (
            internal_volume * params.crystal_density / dt
        ),
        "internal_reaction_clamped": clamped,
        "top_skeleton_kg_m2_s": before.top_skeleton_rate,
        "internal_over_top": (
            before.internal_reaction / before.top_skeleton_rate
            if before.top_skeleton_rate > 0.0
            else math.inf
        ),
        "salt_flux_soil_kg_m2_s": salt_flux,
        "water_flux_soil_kg_m2_s": water_flux,
        "salt_flux_base_route_kg_m2_s": before.salt_flux_base_route,
        "salt_flux_base_route_legacy_kg_m2_s": before.salt_flux_base_route_legacy,
        "route_difference_kg_m2_s": before.salt_flux_base_route - salt_flux,
        "route_difference_relative": abs(before.salt_flux_base_route - salt_flux) / scale,
        "route_difference_legacy_relative": (
            abs(before.salt_flux_base_route_legacy - salt_flux) / scale
        ),
        "route_difference_identity_kg_m2_s": before.route_difference_identity,
        "route_identity_error_relative": (
            abs(
                (before.salt_flux_base_route - salt_flux)
                - before.route_difference_identity
            )
            / max(abs(before.route_difference_identity), 1.0e-300)
        ),
        "solution_mass_flux_top_relative_kg_m2_s": (
            before.solution_mass_flux_top_relative
        ),
        "solution_mass_flux_base_kg_m2_s": before.solution_mass_flux_base,
        "solid_salt_kg_m2": solid_after,
        "dissolved_salt_kg_m2": after.dissolved_salt,
        "water_kg_m2": after.water,
        "cumulative_salt_from_soil_kg_m2": state.cumulative_salt_from_soil,
        "cumulative_water_from_soil_kg_m2": state.cumulative_water_from_soil,
        "cumulative_evaporation_kg_m2": state.cumulative_evaporation,
        "cumulative_internal_reaction_kg_m2": state.cumulative_internal_reaction,
        "cumulative_top_skeleton_kg_m2": state.cumulative_top_skeleton,
        "base_mass_fraction": before.base_mass_fraction,
        "growth_arrested": before.growth_arrested,
        "salt_flux_base_route_corrected_kg_m2_s": before.salt_flux_base_route,
        "sealed": state.sealed,
        "disconnected": state.disconnected,
        "stopped": state.stopped,
        "stop_reason": state.stop_reason,
    }
    state.history.append(record)
    return record


def budget(state: GradedCrustState, initial: GradedCrustState) -> dict:
    """Component budget residuals between two states of the same patch."""
    d_solid = state.solid_salt_mass - initial.solid_salt_mass
    d_dissolved = state.dissolved_salt - initial.dissolved_salt
    d_water = state.water - initial.water
    salt_transfer = state.cumulative_salt_from_soil - initial.cumulative_salt_from_soil
    water_transfer = (
        state.cumulative_water_from_soil - initial.cumulative_water_from_soil
    )
    evaporated = state.cumulative_evaporation - initial.cumulative_evaporation
    swept = state.cumulative_vapour_sweep - initial.cumulative_vapour_sweep
    internal = (
        state.cumulative_internal_reaction - initial.cumulative_internal_reaction
    )
    top = state.cumulative_top_skeleton - initial.cumulative_top_skeleton
    salt_residual = salt_transfer - (d_solid + d_dissolved)
    water_residual = water_transfer - (evaporated - swept + d_water)
    salt_scale = max(abs(salt_transfer), abs(d_solid + d_dissolved), 1.0e-300)
    water_scale = max(abs(water_transfer), abs(evaporated), 1.0e-300)
    solid_scale = max(abs(d_solid), 1.0e-300)
    return {
        "delta_solid_salt_kg_m2": d_solid,
        "delta_dissolved_salt_kg_m2": d_dissolved,
        "delta_water_kg_m2": d_water,
        "integrated_salt_from_soil_kg_m2": salt_transfer,
        "integrated_water_from_soil_kg_m2": water_transfer,
        "integrated_evaporation_kg_m2": evaporated,
        "vapour_swept_kg_m2": swept,
        "internal_reaction_kg_m2": internal,
        "top_skeleton_kg_m2": top,
        "internal_fraction_of_new_solid": internal / max(internal + top, 1.0e-300),
        "solid_split_residual_relative": abs(d_solid - internal - top) / solid_scale,
        "salt_residual_kg_m2": salt_residual,
        "water_residual_kg_m2": water_residual,
        "salt_residual_relative": abs(salt_residual) / salt_scale,
        "water_residual_relative": abs(water_residual) / water_scale,
    }
