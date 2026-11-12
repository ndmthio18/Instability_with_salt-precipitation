"""Conservative prescribed-rate dual-moving-front ALE prototype.

This is an intentionally diagnostic model.  It supplies no predictive law for
bottom dissolution.  Instead it verifies the geometric conservation law and
the salt/water/solid ledgers for a declared pair of front velocities.  The
lower interfacial reaction is independent of, and non-overlapping with, the
signed bulk reaction.
"""

from __future__ import annotations

from dataclasses import dataclass, replace
from typing import Callable

import numpy as np


Array = np.ndarray
ReactionLaw = Callable[[Array, Array, Array, float], Array]


@dataclass(frozen=True)
class BoundaryKinematics:
    top_velocity_m_s: float
    bottom_velocity_m_s: float
    base_dissolved_fraction_to_crust: float = 0.0

    def validate(self) -> None:
        if self.top_velocity_m_s < 0.0 or self.bottom_velocity_m_s < 0.0:
            raise ValueError("the diagnostic tests require nonnegative front velocities")
        if not 0.0 <= self.base_dissolved_fraction_to_crust <= 1.0:
            raise ValueError("base partition must lie in [0,1]")


@dataclass
class DualFrontState:
    z_bottom_m: float
    z_top_m: float
    pore_volume_m: Array
    dissolved_salt_kg_m2: Array
    external_salt_kg_m2: float = 0.0
    external_water_kg_m2: float = 0.0
    cumulative_top_solid_kg_m2: float = 0.0
    cumulative_base_dissolution_kg_m2: float = 0.0
    cumulative_bulk_reaction_kg_m2: float = 0.0
    time_s: float = 0.0

    def copy(self) -> "DualFrontState":
        return replace(
            self,
            pore_volume_m=self.pore_volume_m.copy(),
            dissolved_salt_kg_m2=self.dissolved_salt_kg_m2.copy(),
        )


@dataclass
class _Rates:
    z_bottom: float
    z_top: float
    pore: Array
    salt: Array
    external_salt: float
    external_water: float
    top_solid: float
    base_dissolution: float
    bulk_reaction: float


class DualFrontModel:
    """One-dimensional finite-volume ALE with affine mesh velocity."""

    def __init__(
        self,
        cells: int,
        *,
        liquid_density_kg_m3: float = 1202.168577441632,
        crystal_density_kg_m3: float = 2165.0,
        top_porosity: float = 0.43,
        bottom_porosity: float = 0.43,
        reaction_law: ReactionLaw | None = None,
        bulk_exclusion_fraction: float = 0.0,
        minimum_height_m: float = 1.0e-9,
    ) -> None:
        if cells < 4:
            raise ValueError("cells must be at least four")
        if not 0.0 < top_porosity < 1.0 or not 0.0 < bottom_porosity < 1.0:
            raise ValueError("interface porosities must lie in (0,1)")
        if not 0.0 <= bulk_exclusion_fraction < 1.0:
            raise ValueError("bulk exclusion fraction must lie in [0,1)")
        self.cells = cells
        self.rho = liquid_density_kg_m3
        self.rho_cr = crystal_density_kg_m3
        self.top_porosity = top_porosity
        self.bottom_porosity = bottom_porosity
        self.reaction_law = reaction_law
        self.bulk_exclusion_fraction = bulk_exclusion_fraction
        self.minimum_height_m = minimum_height_m
        self.faces = np.linspace(0.0, 1.0, cells + 1)
        self.widths = np.diff(self.faces)
        self.centers = 0.5 * (self.faces[:-1] + self.faces[1:])

    def initial_state(
        self,
        height_m: float,
        concentration: float | Array,
        porosity: float | Array = 0.43,
        z_bottom_m: float = 0.0,
    ) -> DualFrontState:
        eps = np.broadcast_to(np.asarray(porosity, dtype=float), (self.cells,)).copy()
        c = np.broadcast_to(np.asarray(concentration, dtype=float), (self.cells,)).copy()
        if height_m <= self.minimum_height_m or np.any((eps <= 0.0) | (eps >= 1.0)):
            raise ValueError("invalid initial geometry or porosity")
        if np.any((c < 0.0) | (c >= 1.0)):
            raise ValueError("invalid initial mass fraction")
        pore = eps * height_m * self.widths
        return DualFrontState(
            z_bottom_m=z_bottom_m,
            z_top_m=z_bottom_m + height_m,
            pore_volume_m=pore,
            dissolved_salt_kg_m2=self.rho * c * pore,
        )

    def primitive(self, state: DualFrontState) -> tuple[float, Array, Array]:
        height = state.z_top_m - state.z_bottom_m
        volumes = height * self.widths
        eps = state.pore_volume_m / volumes
        c = state.dissolved_salt_kg_m2 / (self.rho * state.pore_volume_m)
        return height, eps, c

    def inventory(self, state: DualFrontState) -> dict[str, float]:
        height, _, _ = self.primitive(state)
        pore = float(np.sum(state.pore_volume_m))
        dissolved = float(np.sum(state.dissolved_salt_kg_m2))
        water = self.rho * pore - dissolved
        solid = self.rho_cr * (height - pore)
        return {
            "height_m": height,
            "pore_volume_m": pore,
            "dissolved_salt_kg_m2": dissolved,
            "water_kg_m2": water,
            "solid_kg_m2": solid,
            "total_salt_with_external_kg_m2": solid
            + dissolved
            + state.external_salt_kg_m2,
            "total_water_with_external_kg_m2": water + state.external_water_kg_m2,
            "average_porosity": pore / height,
        }

    def _upwind_faces(self, values: Array, transport: Array, left: float, right: float) -> Array:
        result = np.empty(self.cells + 1)
        result[0] = left
        result[-1] = right
        for face in range(1, self.cells):
            result[face] = values[face - 1] if transport[face] >= 0.0 else values[face]
        return result

    def rates(self, state: DualFrontState, motion: BoundaryKinematics) -> _Rates:
        motion.validate()
        height, eps, c = self.primitive(state)
        if height <= self.minimum_height_m:
            raise RuntimeError("crust thickness reached the safety limit")

        vt = motion.top_velocity_m_s
        vb = motion.bottom_velocity_m_s
        w = vb + self.faces * (vt - vb)
        mesh_transport = -w
        eps_face = self._upwind_faces(
            eps, mesh_transport, self.bottom_porosity, self.top_porosity
        )
        h_flux = mesh_transport * eps_face
        c_face = self._upwind_faces(c, h_flux, float(c[0]), float(c[-1]))
        salt_flux = self.rho * c_face * h_flux

        g_top = self.rho_cr * (1.0 - self.top_porosity) * vt
        d_base = self.rho_cr * (1.0 - self.bottom_porosity) * vb
        alpha = motion.base_dissolved_fraction_to_crust
        # The prescribed-rate diagnostic routes the salt needed for top
        # crystallisation from the external/soil-side reservoir as a constant
        # through-crust flux and consumes it at the upper interface.  Thus
        # G_top changes solid inventory but not dissolved inventory by itself.
        # New-pore brine is carried by the geometric flux.
        salt_flux += g_top
        salt_flux[0] += alpha * d_base

        volumes = height * self.widths
        if self.reaction_law is None:
            reaction = np.zeros(self.cells)
        else:
            reaction = np.asarray(
                self.reaction_law(self.centers.copy(), eps.copy(), c.copy(), state.time_s),
                dtype=float,
            )
            if reaction.shape != (self.cells,):
                raise ValueError("reaction law returned an incorrect shape")
        excluded = self.centers < self.bulk_exclusion_fraction
        reaction[excluded] = 0.0
        reaction_mass = reaction * volumes

        pore_rate = h_flux[:-1] - h_flux[1:] - reaction_mass / self.rho_cr
        salt_rate = salt_flux[:-1] - salt_flux[1:] - reaction_mass

        pore_salt_top = self.rho * c_face[-1] * self.top_porosity * vt
        pore_water_top = self.rho * (1.0 - c_face[-1]) * self.top_porosity * vt
        pore_salt_base = self.rho * c_face[0] * self.bottom_porosity * vb
        pore_water_base = self.rho * (1.0 - c_face[0]) * self.bottom_porosity * vb
        external_salt = (
            -g_top - pore_salt_top + pore_salt_base + (1.0 - alpha) * d_base
        )
        # Adding alpha*D_base salt to a saturated, constant-density pore volume
        # displaces the same mass of water to the external reservoir.  This is
        # a declared diagnostic closure, not an activity/density model.
        external_water = (
            -pore_water_top
            + pore_water_base
            + alpha * d_base
            - (1.0 - self.rho / self.rho_cr) * float(np.sum(reaction_mass))
        )

        return _Rates(
            z_bottom=vb,
            z_top=vt,
            pore=pore_rate,
            salt=salt_rate,
            external_salt=external_salt,
            external_water=external_water,
            top_solid=g_top,
            base_dissolution=d_base,
            bulk_reaction=float(np.sum(reaction_mass)),
        )

    @staticmethod
    def _increment(state: DualFrontState, rates: _Rates, dt: float) -> DualFrontState:
        result = state.copy()
        result.z_bottom_m += dt * rates.z_bottom
        result.z_top_m += dt * rates.z_top
        result.pore_volume_m += dt * rates.pore
        result.dissolved_salt_kg_m2 += dt * rates.salt
        result.external_salt_kg_m2 += dt * rates.external_salt
        result.external_water_kg_m2 += dt * rates.external_water
        result.cumulative_top_solid_kg_m2 += dt * rates.top_solid
        result.cumulative_base_dissolution_kg_m2 += dt * rates.base_dissolution
        result.cumulative_bulk_reaction_kg_m2 += dt * rates.bulk_reaction
        result.time_s += dt
        return result

    def step(self, state: DualFrontState, motion: BoundaryKinematics, dt: float) -> DualFrontState:
        """Advance with conservative explicit midpoint/Heun (SSPRK2)."""

        if dt <= 0.0:
            raise ValueError("dt must be positive")
        first = self.rates(state, motion)
        predictor = self._increment(state, first, dt)
        second = self.rates(predictor, motion)
        result = state.copy()
        for name in (
            "z_bottom_m",
            "z_top_m",
            "external_salt_kg_m2",
            "external_water_kg_m2",
            "cumulative_top_solid_kg_m2",
            "cumulative_base_dissolution_kg_m2",
            "cumulative_bulk_reaction_kg_m2",
        ):
            setattr(
                result,
                name,
                getattr(state, name)
                + 0.5 * dt * (getattr(first, self._rate_name(name)) + getattr(second, self._rate_name(name))),
            )
        result.pore_volume_m = state.pore_volume_m + 0.5 * dt * (first.pore + second.pore)
        result.dissolved_salt_kg_m2 = (
            state.dissolved_salt_kg_m2 + 0.5 * dt * (first.salt + second.salt)
        )
        result.time_s = state.time_s + dt
        self.validate(result)
        return result

    @staticmethod
    def _rate_name(state_name: str) -> str:
        return {
            "z_bottom_m": "z_bottom",
            "z_top_m": "z_top",
            "external_salt_kg_m2": "external_salt",
            "external_water_kg_m2": "external_water",
            "cumulative_top_solid_kg_m2": "top_solid",
            "cumulative_base_dissolution_kg_m2": "base_dissolution",
            "cumulative_bulk_reaction_kg_m2": "bulk_reaction",
        }[state_name]

    def validate(self, state: DualFrontState) -> None:
        height, eps, c = self.primitive(state)
        if height <= self.minimum_height_m:
            raise RuntimeError("nonphysical crust thickness")
        if np.any(state.pore_volume_m <= 0.0) or np.any((eps <= 0.0) | (eps >= 1.0)):
            raise RuntimeError("nonphysical pore volume or porosity")
        if np.any((c < -1.0e-12) | (c >= 1.0)):
            raise RuntimeError("nonphysical dissolved-salt mass fraction")

    def run(
        self,
        state: DualFrontState,
        motion: BoundaryKinematics,
        duration_s: float,
        dt_s: float,
    ) -> tuple[DualFrontState, list[dict[str, float]]]:
        if duration_s <= 0.0:
            raise ValueError("duration must be positive")
        current = state.copy()
        history = [self._record(current, motion)]
        while current.time_s < state.time_s + duration_s - 1.0e-14:
            dt = min(dt_s, state.time_s + duration_s - current.time_s)
            current = self.step(current, motion, dt)
            history.append(self._record(current, motion))
        return current, history

    def _record(self, state: DualFrontState, motion: BoundaryKinematics) -> dict[str, float]:
        inv = self.inventory(state)
        return {
            "time_s": state.time_s,
            "z_bottom_m": state.z_bottom_m,
            "z_top_m": state.z_top_m,
            "top_velocity_m_s": motion.top_velocity_m_s,
            "bottom_velocity_m_s": motion.bottom_velocity_m_s,
            "thickness_rate_m_s": motion.top_velocity_m_s - motion.bottom_velocity_m_s,
            **inv,
            "external_salt_kg_m2": state.external_salt_kg_m2,
            "external_water_kg_m2": state.external_water_kg_m2,
            "cumulative_top_solid_kg_m2": state.cumulative_top_solid_kg_m2,
            "cumulative_base_dissolution_kg_m2": state.cumulative_base_dissolution_kg_m2,
            "cumulative_bulk_reaction_kg_m2": state.cumulative_bulk_reaction_kg_m2,
        }
