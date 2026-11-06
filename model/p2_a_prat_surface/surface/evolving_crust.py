"""P2-B: depth-resolved crust with evolving porosity.

Why this module exists
----------------------
P2-A implements the Prat (2026) closed growth law exactly as audited, with a
uniform **frozen** crust porosity (``docs/PRAT_PHYSICS_AUDIT.md`` section 3,
reduction 3). Its verification showed that assumption fails at the declared
parameters: Prat equation (2) forces porosity to change wherever
``C != C_eq``, and the depth-integrated internal reaction exceeds the top
skeleton rate by ``1/sqrt(Da*)``. This module removes the frozen-porosity
assumption and keeps everything else.

Nothing in `prat_crust.py` is modified. That module stays the audited
closed-form reference and the limiting case this one must reproduce.

The material coordinate is already ``z``
----------------------------------------
In Prat's geometry the crust base is fixed at ``z = 0`` and new material is
created only at the moving top ``z = h(t)``. Existing material neither moves
nor deforms. Therefore ``z`` is itself a material coordinate and equation (2)
applies pointwise with no Lagrangian mapping:

    d eps(z,t) / dt = - a_v rho k_r C*(z,t) / rho_cr ,     0 < z < h(t)

New material enters at ``z = h(t)`` with the declared deposition porosity
``eps_0``. A fluid particle's porosity history is thus the time integral of the
local supersaturation it has been exposed to since it was deposited. That is
the deposition-age bookkeeping which `INTERFACE_CONSERVATION.md` section 4 says
is required before a top/internal solid split may be reported, so this module
can report that split legitimately, whereas P2-A could not.

Generalised reduced transport problem
-------------------------------------
With ``eps`` depending on depth, the quasi-steady, ``Pe << Da``, linearised
reaction-diffusion problem of audit section 3 becomes a variable-coefficient
two-point problem rather than a constant-coefficient one:

    d/dz [ eps(z) D_s,eff(z) dC*/dz ] = (1 - C_eq) a_v k_r C* ,
    C*(0) = 0 ,     C*(h) = C*_top .

It is linear and homogeneous in ``C*``, so the solution separates as
``C*(z) = C*_top phi(z)`` where ``phi`` solves the same equation with
``phi(0) = 0`` and ``phi(1) = 1``. One tridiagonal solve per step gives
``phi`` and hence the top gradient coefficient

    A = rho (eps D)_top dphi/dz|_h .

Closing relations, kept EXACT this time
---------------------------------------
P2-A inherited equation (8), which was derived with ``C ~ C_eq`` at the top.
Here equations (6) and (7) are combined without that linearisation. Writing
``C_top = C_eq + s`` and eliminating ``rho U(h)``:

    rho eps_top D_top dC/dz|_h
        = C_top j - rho_cr (1-eps_top)(1-C_top) v_s - C_top rho_v v_s

and the outer-surface kinetic relation of equation (20), in the form verified
in P2-A to be the one compatible with equation (28), is

    v_s = beta k_r rho C*_top / rho_cr  ==  K s .

Substituting gives a single scalar equation for ``s``,

    g(s) = A s - (C_eq + s) j
           + rho_cr (1-eps_top)(1 - C_eq - s) K s
           + (C_eq + s) rho_v K s  =  0 ,

which is quadratic and concave in ``s`` with ``g(0) = -C_eq j < 0`` for
``j > 0``. It is solved by bracketed bisection, so no Newton tuning is needed.
Dropping the ``s^2`` terms and setting ``eps`` uniform and frozen recovers
Prat equation (28) exactly; that reduction is proved algebraically in
`docs/P2_B_EVOLVING_POROSITY.md` and checked numerically in
`verification/test_p2b_evolving.py`.

Declared closures, default off
------------------------------
Reducing porosity physically changes the effective diffusivity and the
specific crystal area. Prat (2026) supplies no closure for either, and the
audit lists both among the unresolved items. This module therefore exposes

    D_s,eff(eps) = D_ref (eps/eps_0) ** diffusivity_exponent
    a_v(eps)     = a_v0  (eps/eps_0) ** specific_area_exponent

with **both exponents defaulting to zero**, i.e. constant, so that the bare
effect of equation (2) is isolated. Non-zero exponents are declared
sensitivities, never fits.

Discretisation
--------------
A fixed grid of spacing ``dz`` on ``z`` with a partial top cell. The state per
cell is the **solid thickness** ``s_i`` (solid volume per unit projected area),
so total solid is ``rho_cr sum(s_i)`` by construction and porosity is a derived
quantity ``eps_i = 1 - s_i / t_i``. Growth adds ``(1-eps_0) v_s dt`` of solid
together with ``v_s dt`` of occupied thickness; the internal reaction adds
``a_v rho k_r C*_i t_i dt / rho_cr``. Both are exact increments of the same
tracked variable, so no porosity update can leak solid mass.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

from .prat_crust import CrustModelDomainExit, CrustParameters
from .thermo import BrineDensityLaw


@dataclass(frozen=True)
class EvolvingCrustParameters:
    """Parameters of the depth-resolved crust.

    The first six entries have the same meaning and units as in
    `prat_crust.CrustParameters`. ``deposition_porosity`` is the porosity of
    newly created top layers, i.e. the ``eps_s`` that Prat's self-replicating
    construction assumes; it is a declared input.
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
    #: Porosity at which the pore network is declared sealed.
    clogging_porosity: float = 0.01
    #: Porosity at which capillary supply is declared lost (audit 5.4).
    connectivity_porosity: float = 0.02
    #: Declared sensitivities; zero means "constant", which is the default.
    diffusivity_exponent: float = 0.0
    specific_area_exponent: float = 0.0
    #: Reinstate the ``C ~ C_eq`` top linearisation of Prat equation (8).
    #: Only for the reduction test: with this True and a uniform frozen
    #: porosity the solver must reproduce equation (28) to discretisation
    #: error. Production runs keep it False, i.e. exact in ``C_top``.
    linearised_top_closure: bool = False

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

    @property
    def rho(self) -> float:
        if self.solution_density is not None:
            return self.solution_density
        return BrineDensityLaw().density(self.equilibrium_mass_fraction)

    # -- declared porosity-dependent closures ------------------------------

    def diffusivity(self, porosity: float) -> float:
        if self.diffusivity_exponent == 0.0:
            return self.effective_salt_diffusivity
        ratio = max(porosity, 0.0) / self.deposition_porosity
        return self.effective_salt_diffusivity * ratio**self.diffusivity_exponent

    def area(self, porosity: float) -> float:
        if self.specific_area_exponent == 0.0:
            return self.specific_area
        ratio = max(porosity, 0.0) / self.deposition_porosity
        return self.specific_area * ratio**self.specific_area_exponent

    # -- the frozen-porosity counterpart -----------------------------------

    def frozen(self) -> CrustParameters:
        """The audited P2-A parameter set this module must reduce to."""
        return CrustParameters(
            equilibrium_mass_fraction=self.equilibrium_mass_fraction,
            crust_porosity=self.deposition_porosity,
            crystal_density=self.crystal_density,
            specific_area=self.specific_area,
            reaction_coefficient=self.reaction_coefficient,
            effective_salt_diffusivity=self.effective_salt_diffusivity,
            beta=self.beta,
            solution_density=self.solution_density,
            vapour_density=self.vapour_density,
            crust_saturation=self.crust_saturation,
        )


def _thomas(lower: list, diag: list, upper: list, rhs: list) -> list:
    """Solve a tridiagonal system in place (Thomas algorithm, no scipy)."""
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


@dataclass
class EvolvingCrustState:
    """Depth-resolved crust state on a fixed grid with a moving top.

    ``solid[i]`` is the solid volume per unit projected area in cell ``i``
    (metres of solid). ``height`` is the position of the moving top.
    """

    cell_size: float
    height: float = 0.0
    solid: list = field(default_factory=list)
    time: float = 0.0
    active: bool = False
    activation_time: float | None = None
    eligibility_time: float | None = None
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
    history: list = field(default_factory=list)

    # -- geometry -----------------------------------------------------------

    def thickness(self, index: int) -> float:
        low = index * self.cell_size
        return min(self.cell_size, max(0.0, self.height - low))

    def centre(self, index: int) -> float:
        return index * self.cell_size + 0.5 * self.thickness(index)

    def occupied(self, minimum_fraction: float = 1.0e-9) -> list:
        """Indices of cells with a non-degenerate occupied thickness."""
        return [
            i
            for i in range(len(self.solid))
            if self.thickness(i) > minimum_fraction * self.cell_size
        ]

    def porosity(self, index: int) -> float:
        t = self.thickness(index)
        if t <= 0.0:
            return 0.0
        return 1.0 - self.solid[index] / t

    def porosity_profile(self) -> list:
        return [self.porosity(i) for i in self.occupied()]

    @property
    def minimum_porosity(self) -> float:
        profile = self.porosity_profile()
        return min(profile) if profile else 1.0

    @property
    def solid_salt(self) -> float:
        """``M_s,solid`` [kg NaCl/m^2]; exact sum of the tracked state."""
        return sum(self.solid)

    @property
    def total_salt(self) -> float:
        return self.solid_salt_mass + self.dissolved_salt

    # These two are filled by `step` so budgets can be taken from the state.
    solid_salt_mass: float = 0.0
    dissolved_salt: float = 0.0
    water: float = 0.0


@dataclass(frozen=True)
class ProfileSolution:
    """One quasi-steady solve of the variable-porosity crust."""

    height: float
    evaporation_flux: float
    top_supersaturation: float
    top_mass_fraction: float
    growth_velocity: float
    top_gradient_coefficient: float      # A = rho (eps D)_top dphi/dz|_h
    profile: list                        # C*(z) at the occupied cell centres
    centres: list
    porosity: list
    thicknesses: list
    internal_reaction: float             # R_int [kg NaCl/m^2/s]
    top_skeleton_rate: float             # (1-eps_0) rho_cr v_s
    dissolved_salt: float                # M_s,d [kg NaCl/m^2]
    water: float                         # M_w,c [kg H2O/m^2]
    solution_mass_flux_top: float        # rho U(h)
    solution_mass_flux_base: float       # rho U(0)
    base_gradient: float                 # dC/dz|_0+
    salt_flux_base_route: float          # independent reconstruction
    minimum_porosity: float
    top_porosity: float


def solve_profile(
    state: EvolvingCrustState,
    params: EvolvingCrustParameters,
    evaporation_flux: float,
) -> ProfileSolution:
    """Quasi-steady crust transport and top closure at the current state."""
    if not state.active or state.height <= 0.0:
        raise CrustModelDomainExit(
            "the crust transport problem requires an active crust with h > 0"
        )
    cells = state.occupied()
    if not cells:
        raise CrustModelDomainExit("no occupied crust cells")

    rho = params.rho
    sat = params.crust_saturation
    c_eq = params.equilibrium_mass_fraction
    lam_cells = []
    eps_cells = []
    t_cells = []
    z_cells = []
    for i in cells:
        eps = state.porosity(i)
        eps_cells.append(eps)
        t_cells.append(state.thickness(i))
        z_cells.append(state.centre(i))
        lam_cells.append((1.0 - c_eq) * params.area(eps) * params.reaction_coefficient)

    n = len(cells)
    cond = [eps_cells[k] * params.diffusivity(eps_cells[k]) for k in range(n)]

    # Finite volume with Dirichlet phi(0) = 0 and phi(h) = 1.
    lower = [0.0] * n
    diag = [0.0] * n
    upper = [0.0] * n
    rhs = [0.0] * n
    for k in range(n):
        # lower face
        if k == 0:
            distance = z_cells[0]
            face = cond[0]
            coefficient = face / distance if distance > 0.0 else 0.0
            diag[k] -= coefficient            # phi_bottom = 0
        else:
            distance = z_cells[k] - z_cells[k - 1]
            face = 0.5 * (cond[k] + cond[k - 1])
            coefficient = face / distance
            diag[k] -= coefficient
            lower[k] += coefficient
        # upper face
        if k == n - 1:
            distance = state.height - z_cells[k]
            face = cond[k]
            coefficient = face / distance if distance > 0.0 else 0.0
            diag[k] -= coefficient
            rhs[k] -= coefficient             # phi_top = 1
        else:
            distance = z_cells[k + 1] - z_cells[k]
            face = 0.5 * (cond[k] + cond[k + 1])
            coefficient = face / distance
            diag[k] -= coefficient
            upper[k] += coefficient
        diag[k] -= lam_cells[k] * t_cells[k]

    phi = _thomas(lower, diag, upper, rhs)

    top_distance = state.height - z_cells[-1]
    if top_distance <= 0.0:
        raise CrustModelDomainExit("degenerate top cell in the crust grid")
    top_gradient = (1.0 - phi[-1]) / top_distance
    eps_top = eps_cells[-1]
    a_coefficient = rho * cond[-1] * top_gradient

    # Scalar top closure, exact in C_top (no C ~ C_eq linearisation).
    k_kinetic = params.beta * params.reaction_coefficient * rho / params.crystal_density
    skeleton = params.crystal_density * (1.0 - params.deposition_porosity)

    def residual(s: float) -> float:
        v_s = k_kinetic * s
        # The linearised branch replaces C_top by C_eq exactly where Prat
        # equation (8) does, and is used only to prove the reduction.
        c_top = c_eq if params.linearised_top_closure else c_eq + s
        return (
            a_coefficient * s
            - c_top * evaporation_flux
            + skeleton * (1.0 - c_top) * v_s
            + c_top * params.vapour_density * v_s
        )

    s_top = _solve_top_supersaturation(residual, c_eq, evaporation_flux)
    v_s = k_kinetic * s_top
    c_top = c_eq + s_top

    c_star = [s_top * value for value in phi]
    internal = 0.0
    dissolved = 0.0
    water = 0.0
    for k in range(n):
        internal += (
            params.area(eps_cells[k])
            * rho
            * params.reaction_coefficient
            * c_star[k]
            * t_cells[k]
        )
        pore = eps_cells[k] * sat * rho * t_cells[k]
        dissolved += pore * (c_eq + c_star[k])
        water += pore * (1.0 - c_eq - c_star[k])

    rho_u_top = (
        rho * eps_top * sat + skeleton - params.vapour_density
    ) * v_s + evaporation_flux
    rho_u_base = rho_u_top + internal
    base_gradient = s_top * phi[0] / z_cells[0] if z_cells[0] > 0.0 else 0.0
    salt_base = (
        c_eq * rho_u_base
        - rho * eps_cells[0] * sat * params.diffusivity(eps_cells[0]) * base_gradient
    )

    return ProfileSolution(
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
        base_gradient=base_gradient,
        salt_flux_base_route=salt_base,
        minimum_porosity=min(eps_cells),
        top_porosity=eps_top,
    )


def _solve_top_supersaturation(residual, c_eq: float, flux: float) -> float:
    """Bracketed bisection for the scalar top closure."""
    if flux <= 0.0:
        return 0.0
    low = 0.0
    g_low = residual(low)
    if g_low >= 0.0:
        return low
    high = max(1.0e-12, 1.0e-3 * c_eq)
    for _ in range(200):
        if residual(high) >= 0.0:
            break
        high *= 2.0
        if c_eq + high >= 1.0:
            raise CrustModelDomainExit(
                "the crust-top closure has no root below C_top = 1; the "
                "linearised reaction closure is invalid here"
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
    state: EvolvingCrustState,
    params: EvolvingCrustParameters,
    seed_height: float,
    evaporation_flux: float,
    time: float,
) -> EvolvingCrustState:
    """Create the seed crust at the declared deposition porosity.

    As in P2-A, the seed inventory is charged to the soil so that component
    budgets close across the activation event.
    """
    if state.active:
        return state
    if seed_height <= 0.0:
        raise ValueError("seed height must be positive; h = 0 is not a crust")
    state.height = seed_height
    state.active = True
    state.activation_time = time
    count = max(1, math.ceil(seed_height / state.cell_size - 1.0e-12))
    state.solid = [0.0] * count
    for i in range(count):
        state.solid[i] = (1.0 - params.deposition_porosity) * state.thickness(i)
    solution = solve_profile(state, params, evaporation_flux)
    state.solid_salt_mass = params.crystal_density * state.solid_salt
    state.dissolved_salt = solution.dissolved_salt
    state.water = solution.water
    state.seed_salt = state.solid_salt_mass + state.dissolved_salt
    state.seed_water = state.water
    state.cumulative_salt_from_soil += state.seed_salt
    state.cumulative_water_from_soil += state.seed_water
    return state


# --- time integration ---------------------------------------------------------


def step(
    state: EvolvingCrustState,
    params: EvolvingCrustParameters,
    evaporation_flux: float,
    dt: float,
    maximum_cells: int = 100000,
) -> dict:
    """Advance porosity, the moving top and the inventories by ``dt``.

    The reported base fluxes are the **inventory (balance) route**, evaluated
    from the actual change of the tracked states, so

        F_s,soil = d(M_s,solid + M_s,d)/dt
        F_w,soil = j - rho_v v_s + dM_w,c/dt

    hold exactly. The independent base-flux reconstruction of
    `INTERFACE_CONSERVATION.md` section 3 is evaluated alongside and reported
    as ``salt_flux_base_route`` so the two can be compared.
    """
    if dt <= 0.0 or not math.isfinite(dt):
        raise ValueError("time step must be positive and finite")
    if not state.active:
        raise CrustModelDomainExit("cannot step an inactive crust")
    if state.sealed:
        raise CrustModelDomainExit(
            f"crust sealed at t = {state.sealed_time!r}; the saturated "
            "capillary-supply closure is disabled rather than extrapolated"
        )

    before = solve_profile(state, params, evaporation_flux)
    solid_before = params.crystal_density * state.solid_salt
    dissolved_before = before.dissolved_salt
    water_before = before.water

    cells = state.occupied()
    rho = params.rho

    # 1. internal reaction, cell by cell (equation 2, pointwise in z)
    internal_volume = 0.0
    clamped = False
    for position, index in enumerate(cells):
        rate = (
            params.area(before.porosity[position])
            * rho
            * params.reaction_coefficient
            * before.profile[position]
        )
        added = rate * before.thicknesses[position] * dt / params.crystal_density
        # Porosity cannot go below zero. The clamp can only bind on the step
        # that seals a cell; it is reported rather than hidden, and the
        # budgets use the clamped increment so conservation is unaffected.
        available = before.thicknesses[position] - state.solid[index]
        if added > available:
            added = max(0.0, available)
            clamped = True
        state.solid[index] += added
        internal_volume += added

    # 2. moving top: new material at the declared deposition porosity
    growth = before.growth_velocity * dt
    new_height = state.height + growth
    while math.ceil(new_height / state.cell_size - 1.0e-12) > len(state.solid):
        if len(state.solid) >= maximum_cells:
            raise CrustModelDomainExit("crust grid exceeded the declared cell cap")
        state.solid.append(0.0)
    remaining = growth
    cursor = state.height
    while remaining > 1.0e-18:
        cell = min(int(cursor / state.cell_size + 1.0e-12), len(state.solid) - 1)
        cell_top = (cell + 1) * state.cell_size
        take = min(remaining, max(0.0, cell_top - cursor))
        if take <= 0.0:
            cursor = cell_top
            continue
        state.solid[cell] += (1.0 - params.deposition_porosity) * take
        cursor += take
        remaining -= take
    state.height = new_height
    state.time += dt

    after = solve_profile(state, params, evaporation_flux)
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
        (1.0 - params.deposition_porosity) * params.crystal_density * growth
    )

    minimum = after.minimum_porosity
    if minimum <= params.clogging_porosity and not state.sealed:
        state.sealed = True
        state.sealed_time = state.time

    record = {
        "time_s": state.time,
        "dt_s": dt,
        "crust_height_m": state.height,
        "cells": len(state.occupied()),
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
        "route_difference_kg_m2_s": before.salt_flux_base_route - salt_flux,
        "solid_salt_kg_m2": solid_after,
        "dissolved_salt_kg_m2": after.dissolved_salt,
        "water_kg_m2": after.water,
        "cumulative_salt_from_soil_kg_m2": state.cumulative_salt_from_soil,
        "cumulative_water_from_soil_kg_m2": state.cumulative_water_from_soil,
        "cumulative_evaporation_kg_m2": state.cumulative_evaporation,
        "cumulative_internal_reaction_kg_m2": state.cumulative_internal_reaction,
        "cumulative_top_skeleton_kg_m2": state.cumulative_top_skeleton,
        "sealed": state.sealed,
        "connected": minimum > params.connectivity_porosity,
    }
    state.history.append(record)
    return record


def budget(state: EvolvingCrustState, initial: EvolvingCrustState) -> dict:
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
        "seed_salt_kg_m2": state.seed_salt - initial.seed_salt,
        "seed_water_kg_m2": state.seed_water - initial.seed_water,
        "internal_reaction_kg_m2": (
            state.cumulative_internal_reaction - initial.cumulative_internal_reaction
        ),
        "top_skeleton_kg_m2": (
            state.cumulative_top_skeleton - initial.cumulative_top_skeleton
        ),
        "solid_split_residual_relative": abs(
            d_solid
            - (state.cumulative_internal_reaction - initial.cumulative_internal_reaction)
            - (state.cumulative_top_skeleton - initial.cumulative_top_skeleton)
        )
        / solid_scale,
        "salt_residual_kg_m2": salt_residual,
        "water_residual_kg_m2": water_residual,
        "salt_residual_relative": abs(salt_residual) / salt_scale,
        "water_residual_relative": abs(water_residual) / water_scale,
    }
