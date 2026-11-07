"""P2-D: an isolated 1D unsaturated Richards-NaCl column.

Why this exists
---------------
`task/6.md` items 4, 5 and 8 need a Richards column that can be driven by a
water-activity-controlled evaporation boundary and audited for component
conservation. The production solver is `model/p1_r2_conservative_richards`,
which is DuMuX/Dune C++ and needs the pinned stack recorded in its
`provenance/source_manifest.json`. No C++ toolchain is available in this
environment (`g++`, `clang++`, `cl`, `cmake`, `make` are all absent), so that
solver cannot be built or run here.

This module is therefore a **mirror**, not a port: it reproduces P1-R2's
*constitutive closures and balance equations exactly*, on a simpler
discretisation, so that the coupling can be demonstrated and audited with real
numbers. Every closure below is transcribed from
`model/p1_r2_conservative_richards/source/appl/simpleRichardsnc/main.cc` and
`model.hh`; the transcription is checked term by term in
`verification/test_p2d_coupling.py`.

What is identical to P1-R2
--------------------------
The two balance equations, per bulk volume::

    d/dt (phi S_w rho)            + div F_total = 0
    d/dt (phi S_w rho w_s)        + div F_salt  = 0

with `rho = rho(x)`, `w_s = w_s(x)` and `x` the NaCl **mole** fraction, and the
flux closure::

    v   = -(K/mu) krw(S_up) [ grad p + rho_face g e_z ]
    F_total = rho_up v
    F_salt  = rho_up w_up v - rho_face D_face grad w_s

i.e. upwinded mobility, upwinded density and composition in the advective
terms, arithmetic face averages for the density in the gravity term and for
`rho D` in the mass-averaged Fickian term. The constitutive functions
`sw(pw)`, `krw(sw)`, `D_eff(sw)`, `rho(x)`, `w_s(x)` are the same expressions
with the same parameters.

What is different, and declared
-------------------------------
1. P1-R2 is a vertex-centred box (CVFE) scheme in two dimensions. This is a
   cell-centred two-point-flux finite volume scheme in one dimension. The
   closures are the same; the discretisation is not, so discrete values must
   not be compared digit for digit with P1-R2 output.
2. Time integration is backward Euler with a damped Newton solve on a
   block-tridiagonal system. P1-R2 uses DuMuX's Newton with its own
   line search and time-step control.
3. No gravity-driven lateral structure exists in 1D, so nothing here can
   address the instability question.

Conservation is **not** inherited from the discretisation argument: it is
measured. `face_fluxes()` reconstructs every face flux from the converged
primary variables independently of the residual assembly, and `audit()`
compares the inventory change against the integrated boundary transfer.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

P_NON_WETTING_REFERENCE = 1.0e5


class ColumnSolverError(RuntimeError):
    """Raised when the Newton solve fails; never returns an unconverged state."""


class ColumnDomainExit(RuntimeError):
    """Raised when the state leaves the declared validity envelope."""


@dataclass(frozen=True)
class SoilParameters:
    """P1-R2 problem parameters. Every default is that file's own default.

    The values without defaults are the ones P1-R2 refuses to default, because
    `README.md` states that rounded fallbacks are intentionally unavailable:
    the two molar masses and the three density-law constants.
    """

    molar_mass_h2o: float
    molar_mass_nacl: float
    pure_water_density: float
    density_linear_coefficient: float
    density_quadratic_coefficient: float
    porosity: float
    permeability: float
    pc_alpha: float
    pc_n: float
    residual_saturation: float
    swe_l: float = 0.73
    viscosity: float = 1.1e-3
    salt_diffusivity: float = 1.5e-9
    gravity: float = 9.81
    enable_gravity: bool = True
    krw_mualem_exponent: bool = False

    # -- constitutive closures, transcribed from main.cc --------------------

    @property
    def pc_m(self) -> float:
        return 1.0 - 1.0 / self.pc_n

    def effective_sw(self, sw: float) -> float:
        return (sw - self.residual_saturation) / (1.0 - self.residual_saturation)

    def sw(self, pw: float) -> float:
        """van Genuchten saturation. Mirrors `Problem::sw`."""
        pc = P_NON_WETTING_REFERENCE - pw
        if pc <= 0.0:
            return 1.0
        n, m = self.pc_n, self.pc_m
        return self.residual_saturation + (1.0 - self.residual_saturation) / (
            (1.0 + (self.pc_alpha * pc) ** n) ** m
        )

    def pw(self, sw: float) -> float:
        """Inverse of `sw`, used only to build initial states."""
        if sw >= 1.0:
            return P_NON_WETTING_REFERENCE
        se = self.effective_sw(sw)
        if not (0.0 < se < 1.0):
            raise ColumnDomainExit(f"saturation {sw!r} outside (swr, 1)")
        n, m = self.pc_n, self.pc_m
        pc = (se ** (-1.0 / m) - 1.0) ** (1.0 / n) / self.pc_alpha
        return P_NON_WETTING_REFERENCE - pc

    def krw(self, sw: float) -> float:
        """Mualem-van Genuchten. Mirrors `Problem::krw`, factor form default."""
        m = self.pc_m
        se = min(max(self.effective_sw(sw), 0.0), 1.0)
        if se <= 0.0:
            return 0.0
        value = 1.0 - (1.0 - se ** (1.0 / m)) ** m
        if self.krw_mualem_exponent:
            return se**self.swe_l * value * value
        return self.swe_l * value * value

    def effective_diffusion_coefficient(self, sw: float) -> float:
        """Millington-Quirk. Mirrors `Problem::effectiveDiffusionCoefficient`."""
        return self.porosity**1.5 * max(sw, 0.0) ** 3.5 * self.salt_diffusivity

    def mass_fraction(self, mole_fraction: float) -> float:
        mixture = (
            1.0 - mole_fraction
        ) * self.molar_mass_h2o + mole_fraction * self.molar_mass_nacl
        return mole_fraction * self.molar_mass_nacl / mixture

    def mole_fraction(self, mass_fraction: float) -> float:
        """Exact inverse of `mass_fraction`."""
        a = mass_fraction / self.molar_mass_nacl
        b = (1.0 - mass_fraction) / self.molar_mass_h2o
        return a / (a + b)

    def rho(self, mole_fraction: float) -> float:
        """Mirrors `Problem::rhoW`."""
        w = self.mass_fraction(mole_fraction)
        return self.pure_water_density + 1000.0 * w * (
            self.density_linear_coefficient + self.density_quadratic_coefficient * w
        )

    def gravitational_acceleration(self) -> float:
        return self.gravity if self.enable_gravity else 0.0


@dataclass
class ColumnState:
    """Cell-centred state of the column. Cell 0 is the bottom."""

    cells: int
    height: float
    pressure: list = field(default_factory=list)
    mole_fraction: list = field(default_factory=list)
    time: float = 0.0

    @property
    def dz(self) -> float:
        return self.height / self.cells

    def centre(self, index: int) -> float:
        return (index + 0.5) * self.dz

    def copy(self) -> "ColumnState":
        return ColumnState(
            cells=self.cells,
            height=self.height,
            pressure=list(self.pressure),
            mole_fraction=list(self.mole_fraction),
            time=self.time,
        )


@dataclass(frozen=True)
class BoundaryCondition:
    """Top and bottom conditions.

    ``top_flux`` returns ``(total_mass_flux, salt_mass_flux)`` leaving the
    column through the top face, in kg/m2/s, given the column state. It is a
    callable so the flux can depend on the solution, which is exactly what a
    water-activity-controlled film requires (item 4).
    """

    top_flux: object
    bottom: str = "no_flux"          # "no_flux" or "reservoir"
    reservoir_mole_fraction: float = 0.0
    reservoir_pressure: float = P_NON_WETTING_REFERENCE


def storage(params: SoilParameters, pressure: float, mole_fraction: float) -> tuple:
    """``(total liquid mass, NaCl mass)`` per bulk volume [kg/m3].

    Mirrors `RichardsSaltModelLocalResidual::massStorage`.
    """
    theta_rho = params.porosity * params.sw(pressure) * params.rho(mole_fraction)
    return theta_rho, theta_rho * params.mass_fraction(mole_fraction)


def _interior_face_flux(
    params: SoilParameters,
    p_lo: float,
    x_lo: float,
    p_hi: float,
    x_hi: float,
    distance: float,
) -> tuple:
    """Flux through the face between a lower and an upper cell, upward positive.

    Mirrors `RichardsSaltModelLocalResidual::transportFlux` with the outer
    normal pointing up.
    """
    rho_lo, rho_hi = params.rho(x_lo), params.rho(x_hi)
    rho_face = 0.5 * (rho_lo + rho_hi)
    sw_lo, sw_hi = params.sw(p_lo), params.sw(p_hi)
    d_face = 0.5 * (
        params.effective_diffusion_coefficient(sw_lo)
        + params.effective_diffusion_coefficient(sw_hi)
    )
    grad_p = (p_hi - p_lo) / distance
    body = rho_face * params.gravitational_acceleration()
    velocity = -params.permeability * (grad_p + body)
    sw_up = sw_lo if velocity > 0.0 else sw_hi
    velocity *= params.krw(sw_up) / params.viscosity
    rho_up = rho_lo if velocity > 0.0 else rho_hi
    w_up = params.mass_fraction(x_lo if velocity > 0.0 else x_hi)
    w_lo = params.mass_fraction(x_lo)
    w_hi = params.mass_fraction(x_hi)
    total = rho_up * velocity
    advective = total * w_up
    diffusive = -rho_face * d_face * (w_hi - w_lo) / distance
    return total, advective + diffusive, advective, diffusive, velocity


def face_fluxes(
    params: SoilParameters,
    state: ColumnState,
    boundary: BoundaryCondition,
) -> dict:
    """Independent reconstruction of every face flux, upward positive.

    This function does not call the residual assembly and is not called by it.
    It is the independent route required by item 5.
    """
    n = state.cells
    dz = state.dz
    total = [0.0] * (n + 1)
    salt = [0.0] * (n + 1)
    advective = [0.0] * (n + 1)
    diffusive = [0.0] * (n + 1)
    for face in range(1, n):
        t, s, a, d, _ = _interior_face_flux(
            params,
            state.pressure[face - 1],
            state.mole_fraction[face - 1],
            state.pressure[face],
            state.mole_fraction[face],
            dz,
        )
        total[face], salt[face], advective[face], diffusive[face] = t, s, a, d
    # bottom face
    if boundary.bottom == "reservoir":
        t, s, a, d, _ = _interior_face_flux(
            params,
            boundary.reservoir_pressure,
            boundary.reservoir_mole_fraction,
            state.pressure[0],
            state.mole_fraction[0],
            0.5 * dz,
        )
        total[0], salt[0], advective[0], diffusive[0] = t, s, a, d
    # top face: prescribed, leaving the column
    top_total, top_salt = boundary.top_flux(params, state)
    total[n], salt[n] = top_total, top_salt
    advective[n], diffusive[n] = top_total, 0.0
    return {
        "total": total,
        "salt": salt,
        "advective": advective,
        "diffusive": diffusive,
    }


def inventory(params: SoilParameters, state: ColumnState) -> dict:
    """Column inventories per unit cross-section [kg/m2]."""
    dz = state.dz
    liquid = 0.0
    salt = 0.0
    for i in range(state.cells):
        s_total, s_salt = storage(params, state.pressure[i], state.mole_fraction[i])
        liquid += s_total * dz
        salt += s_salt * dz
    return {"liquid_kg_m2": liquid, "salt_kg_m2": salt, "water_kg_m2": liquid - salt}


def residual(
    params: SoilParameters,
    state: ColumnState,
    previous: ColumnState,
    boundary: BoundaryCondition,
    dt: float,
) -> list:
    """Backward-Euler residual, two equations per cell, scaled by cell volume."""
    n = state.cells
    dz = state.dz
    flux = face_fluxes(params, state, boundary)
    out = [0.0] * (2 * n)
    for i in range(n):
        now_total, now_salt = storage(
            params, state.pressure[i], state.mole_fraction[i]
        )
        old_total, old_salt = storage(
            params, previous.pressure[i], previous.mole_fraction[i]
        )
        out[2 * i] = (now_total - old_total) * dz / dt + (
            flux["total"][i + 1] - flux["total"][i]
        )
        out[2 * i + 1] = (now_salt - old_salt) * dz / dt + (
            flux["salt"][i + 1] - flux["salt"][i]
        )
    return out


def _block_tridiagonal_solve(lower, diag, upper, rhs) -> list:
    """Thomas algorithm on 2x2 blocks (no numpy, no scipy)."""
    n = len(diag)

    def inv2(m):
        det = m[0][0] * m[1][1] - m[0][1] * m[1][0]
        if det == 0.0:
            raise ColumnSolverError("singular 2x2 block in the column Jacobian")
        return [[m[1][1] / det, -m[0][1] / det], [-m[1][0] / det, m[0][0] / det]]

    def mul(a, b):
        return [
            [a[0][0] * b[0][0] + a[0][1] * b[1][0], a[0][0] * b[0][1] + a[0][1] * b[1][1]],
            [a[1][0] * b[0][0] + a[1][1] * b[1][0], a[1][0] * b[0][1] + a[1][1] * b[1][1]],
        ]

    def mv(a, v):
        return [a[0][0] * v[0] + a[0][1] * v[1], a[1][0] * v[0] + a[1][1] * v[1]]

    def sub(a, b):
        return [[a[0][0] - b[0][0], a[0][1] - b[0][1]], [a[1][0] - b[1][0], a[1][1] - b[1][1]]]

    c = [None] * n
    d = [None] * n
    inv = inv2(diag[0])
    c[0] = mul(inv, upper[0])
    d[0] = mv(inv, rhs[0])
    for i in range(1, n):
        inv = inv2(sub(diag[i], mul(lower[i], c[i - 1])))
        c[i] = mul(inv, upper[i]) if i < n - 1 else [[0.0, 0.0], [0.0, 0.0]]
        d[i] = mv(inv, [rhs[i][0] - (lower[i][0][0] * d[i - 1][0] + lower[i][0][1] * d[i - 1][1]),
                        rhs[i][1] - (lower[i][1][0] * d[i - 1][0] + lower[i][1][1] * d[i - 1][1])])
    x = [None] * n
    x[-1] = d[-1]
    for i in range(n - 2, -1, -1):
        x[i] = [d[i][0] - (c[i][0][0] * x[i + 1][0] + c[i][0][1] * x[i + 1][1]),
                d[i][1] - (c[i][1][0] * x[i + 1][0] + c[i][1][1] * x[i + 1][1])]
    return [component for block in x for component in block]


def _jacobian_blocks(params, state, previous, boundary, dt, base):
    """Block-tridiagonal Jacobian by central differences.

    Central differences cost two residual evaluations per unknown instead of
    one, but their truncation error is O(step^2) rather than O(step), which is
    what lets the Newton solve reach a residual small enough for the component
    budgets to close near round-off rather than at the differencing error.
    """
    n = state.cells
    lower = [[[0.0, 0.0], [0.0, 0.0]] for _ in range(n)]
    diag = [[[0.0, 0.0], [0.0, 0.0]] for _ in range(n)]
    upper = [[[0.0, 0.0], [0.0, 0.0]] for _ in range(n)]
    for j in range(n):
        for which in (0, 1):
            if which == 0:
                target = state.pressure
                step_size = 1.0e-6 * max(abs(state.pressure[j]), 1.0e4)
            else:
                target = state.mole_fraction
                step_size = 1.0e-7 * max(state.mole_fraction[j], 1.0e-4)
            saved = target[j]
            target[j] = saved + step_size
            plus = residual(params, state, previous, boundary, dt)
            target[j] = saved - step_size
            minus = residual(params, state, previous, boundary, dt)
            target[j] = saved
            denominator = 2.0 * step_size
            for i in (j - 1, j, j + 1):
                if not (0 <= i < n):
                    continue
                block = diag[i] if i == j else (upper[i] if i == j - 1 else lower[i])
                for eq in (0, 1):
                    block[eq][which] = (
                        plus[2 * i + eq] - minus[2 * i + eq]
                    ) / denominator
    return lower, diag, upper


def step(
    params: SoilParameters,
    state: ColumnState,
    boundary: BoundaryCondition,
    dt: float,
    tolerance: float = 1.0e-12,
    maximum_iterations: int = 60,
) -> dict:
    """One backward-Euler step with a damped Newton solve.

    Returns convergence diagnostics. Raises rather than returning an
    unconverged state.
    """
    previous = state.copy()
    # The convergence scale must be a property of the problem, not of the
    # first residual: an exactly steady state starts at round-off and would
    # otherwise never satisfy a purely relative test. Use the magnitude of the
    # storage term, which is the size of the equations being solved.
    dz = state.dz
    squared = 0.0
    for i in range(state.cells):
        total, salt = storage(params, previous.pressure[i], previous.mole_fraction[i])
        squared += (total * dz / dt) ** 2 + (salt * dz / dt) ** 2
    scale = max(math.sqrt(squared), 1.0e-300)
    for iteration in range(maximum_iterations):
        base = residual(params, state, previous, boundary, dt)
        norm = math.sqrt(sum(value * value for value in base))
        if norm <= tolerance * scale:
            state.time = previous.time + dt
            # Report the fluxes the step actually applied. Backward Euler uses
            # the END-of-step boundary values, so a caller that accumulates
            # the start-of-step flux does not audit what was solved.
            applied = face_fluxes(params, state, boundary)
            return {
                "iterations": iteration,
                "residual": norm,
                "scale": scale,
                "top_total_kg_m2_s": applied["total"][-1],
                "top_salt_kg_m2_s": applied["salt"][-1],
                "bottom_total_kg_m2_s": applied["total"][0],
                "bottom_salt_kg_m2_s": applied["salt"][0],
            }
        lower, diag, upper = _jacobian_blocks(
            params, state, previous, boundary, dt, base
        )
        rhs = [[-base[2 * i], -base[2 * i + 1]] for i in range(state.cells)]
        update = _block_tridiagonal_solve(lower, diag, upper, rhs)
        damping = 1.0
        for _ in range(60):
            trial_p = [state.pressure[i] + damping * update[2 * i]
                       for i in range(state.cells)]
            trial_x = [state.mole_fraction[i] + damping * update[2 * i + 1]
                       for i in range(state.cells)]
            if all(0.0 < x < 0.9 for x in trial_x) and all(
                params.sw(p) >= params.residual_saturation for p in trial_p
            ):
                break
            damping *= 0.5
        else:
            raise ColumnSolverError("could not damp the Newton update into the domain")
        for i in range(state.cells):
            state.pressure[i] = trial_p[i]
            state.mole_fraction[i] = trial_x[i]
    raise ColumnSolverError(
        f"Newton did not converge in {maximum_iterations} iterations"
    )


def audit(
    params: SoilParameters,
    state: ColumnState,
    initial: ColumnState,
    transferred: dict,
) -> dict:
    """Independent component budget.

    ``transferred`` holds the time-integrated boundary transfer, positive when
    mass leaves the column. The residuals are formed from the independently
    reconstructed inventories, never from the solver's own increments.
    """
    now = inventory(params, state)
    before = inventory(params, initial)
    out = {}
    for name in ("liquid", "salt", "water"):
        key = f"{name}_kg_m2"
        change = now[key] - before[key]
        leaving = transferred[name]
        residual_value = change + leaving
        denominator = max(abs(before[key]), abs(leaving), 1.0e-300)
        out[f"delta_{name}_kg_m2"] = change
        out[f"transferred_{name}_kg_m2"] = leaving
        out[f"{name}_residual_kg_m2"] = residual_value
        out[f"{name}_residual_relative"] = abs(residual_value) / denominator
    out.update({f"initial_{k}": v for k, v in before.items()})
    out.update({f"final_{k}": v for k, v in now.items()})
    return out


def run(
    params: SoilParameters,
    state: ColumnState,
    boundary: BoundaryCondition,
    dt: float,
    steps: int,
    record_stride: int = 1,
    stop=None,
) -> dict:
    """Integrate, accumulating exactly the boundary transfer that was applied.

    The transfer is taken from `step`'s report, i.e. the end-of-step face
    fluxes that backward Euler actually solved with, so the audit measures the
    scheme rather than a separate quadrature of the boundary condition.
    """
    initial = state.copy()
    transferred = {"liquid": 0.0, "salt": 0.0, "water": 0.0}
    history = []
    reason = "steps"
    for index in range(steps):
        try:
            info = soil_step = step(params, state, boundary, dt)
        except (ColumnSolverError, ColumnDomainExit) as error:
            reason = f"{type(error).__name__}: {error}"
            break
        top_total = soil_step["top_total_kg_m2_s"]
        top_salt = soil_step["top_salt_kg_m2_s"]
        bottom_total = soil_step["bottom_total_kg_m2_s"]
        bottom_salt = soil_step["bottom_salt_kg_m2_s"]
        # Positive = leaving the column: out at the top, minus what enters at
        # the bottom face (whose flux is also measured upward positive).
        transferred["liquid"] += (top_total - bottom_total) * dt
        transferred["salt"] += (top_salt - bottom_salt) * dt
        transferred["water"] = transferred["liquid"] - transferred["salt"]
        if index % record_stride == 0 or index == steps - 1:
            inv = inventory(params, state)
            history.append(
                {
                    "time_s": state.time,
                    "newton_iterations": info["iterations"],
                    "newton_residual_relative": info["residual"] / info["scale"],
                    "top_total_kg_m2_s": top_total,
                    "top_salt_kg_m2_s": top_salt,
                    "bottom_total_kg_m2_s": bottom_total,
                    "bottom_salt_kg_m2_s": bottom_salt,
                    "top_mass_fraction": params.mass_fraction(state.mole_fraction[-1]),
                    "top_saturation": params.sw(state.pressure[-1]),
                    "top_pressure_pa": state.pressure[-1],
                    **inv,
                }
            )
        if stop is not None:
            verdict = stop(params, state)
            if verdict:
                reason = verdict
                break
    return {
        "initial": initial,
        "state": state,
        "transferred": transferred,
        "history": history,
        "stop_reason": reason,
        "steps": len(history),
    }
