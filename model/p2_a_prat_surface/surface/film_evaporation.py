"""Physical gas-film evaporation flux (reused from the heterogeneity project).

Formulation source: `problem.hh::filmWaterFlux_` of the heterogeneity cylinder
revision v12. In that code the evaporative water mass flux through a stagnant
gas film of thickness `delta` is

    j = rho_ref * D_v / delta * ( X_w(y_s) - X_w(y_ref) )            [kg/m2/s]

with

    y_s   = a_w * p_sat(T) / p_atm        interfacial vapour mole fraction
    y_ref = RH * p_sat(T) / p_atm         ambient vapour mole fraction
    X_w(y) = M_H2O*y / (M_H2O*y + M_air*(1-y))      gas water mass fraction
    rho_ref = gas density at the ambient reference composition y_ref
    D_v    = binary H2O-air gas diffusion coefficient at (T, p_atm)

This is a mass-fraction-gradient (mass-averaged) Fickian film closure, the same
reference frame used by the P1-R2 liquid-phase closure
``j_s = -rho_l D_eff grad(w_s)``. The two therefore compose without a frame
change.

Differences from the upstream heterogeneity call, all required by P2-A:

* The upstream `Revision.WaterFluxAreaCorrection` factor and the partial-step
  activation ramp are experiment-specific and are NOT reproduced.
* No calibration against an initially prescribed flux is performed
  (P2-A item 2). `delta` is a declared physical input.
* No normalisation of local fluxes to a fixed total evaporation rate is applied
  (P2-A item 4). The total evaporation rate is free to change with salinity.

The flux is signed. A negative value is condensation and is returned as such
rather than silently clipped, following the upstream comment
"Signed flux permits condensation; no normalization to observed total loss."
"""

from __future__ import annotations

import math
from dataclasses import dataclass

from . import thermo
from .water_activity import (
    water_activity_from_mass_fraction,
    water_activity_from_mole_fraction,
)


@dataclass(frozen=True)
class GasFilm:
    """Declared physical inputs of the gas-film evaporation model.

    Attributes
    ----------
    temperature:
        Isothermal film/surface temperature [K].
    relative_humidity:
        Ambient relative humidity, `RH` in [0, 1).
    thickness:
        Gas-film (boundary-layer) thickness `delta` [m].
    atmospheric_pressure:
        Total gas pressure [Pa].
    vapour_diffusion_coefficient:
        `D_v` [m^2/s]. `None` means "evaluate the documented DuMuX H2O-air
        correlation at (T, p_atm)"; an explicit number overrides it.
    """

    temperature: float = 293.15
    relative_humidity: float = 0.4
    thickness: float = 0.0083
    atmospheric_pressure: float = 1.0e5
    vapour_diffusion_coefficient: float | None = None

    def __post_init__(self) -> None:
        if not math.isfinite(self.temperature) or self.temperature <= 0.0:
            raise ValueError("gas film requires a positive finite temperature")
        if not (0.0 <= self.relative_humidity < 1.0):
            raise ValueError("gas film requires relative humidity in [0, 1)")
        if not math.isfinite(self.thickness) or self.thickness <= 0.0:
            raise ValueError("gas film requires a positive finite thickness")
        if (
            not math.isfinite(self.atmospheric_pressure)
            or self.atmospheric_pressure <= 0.0
        ):
            raise ValueError("gas film requires a positive finite pressure")
        if self.vapour_diffusion_coefficient is not None and (
            not math.isfinite(self.vapour_diffusion_coefficient)
            or self.vapour_diffusion_coefficient <= 0.0
        ):
            raise ValueError("explicit vapour diffusion coefficient must be positive")

    # -- derived reference state ------------------------------------------

    @property
    def saturation_pressure(self) -> float:
        return thermo.saturation_pressure(self.temperature)

    @property
    def diffusivity(self) -> float:
        if self.vapour_diffusion_coefficient is not None:
            return self.vapour_diffusion_coefficient
        return thermo.vapour_diffusion_coefficient(
            self.temperature, self.atmospheric_pressure
        )

    @property
    def ambient_vapour_mole_fraction(self) -> float:
        return (
            self.relative_humidity
            * self.saturation_pressure
            / self.atmospheric_pressure
        )

    @property
    def ambient_vapour_mass_fraction(self) -> float:
        return thermo.gas_water_mass_fraction(self.ambient_vapour_mole_fraction)

    @property
    def reference_gas_density(self) -> float:
        """rho_ref, evaluated at the ambient composition, as upstream."""
        return thermo.ideal_gas_density(
            self.temperature,
            self.atmospheric_pressure,
            self.ambient_vapour_mole_fraction,
        )

    @property
    def transfer_coefficient(self) -> float:
        """rho_ref * D_v / delta  [kg m^-2 s^-1 per unit mass-fraction difference]."""
        return self.reference_gas_density * self.diffusivity / self.thickness

    # -- flux --------------------------------------------------------------

    def interfacial_vapour_mole_fraction(self, activity: float) -> float:
        y = activity * self.saturation_pressure / self.atmospheric_pressure
        if not (0.0 <= y < 1.0):
            raise ValueError(
                f"invalid film interfacial vapour state: y_s = {y!r} outside [0, 1)"
            )
        return y

    def driving_mass_fraction_difference(self, activity: float) -> float:
        """X_w(y_s) - X_w(y_ref); signed, positive means evaporation."""
        y_s = self.interfacial_vapour_mole_fraction(activity)
        return thermo.gas_water_mass_fraction(y_s) - self.ambient_vapour_mass_fraction

    def flux_from_activity(self, activity: float) -> float:
        """Water-component mass flux [kg H2O m^-2 s^-1]; signed."""
        if not math.isfinite(activity) or activity < 0.0 or activity > 1.0:
            raise ValueError("water activity must lie in [0, 1]")
        return self.transfer_coefficient * self.driving_mass_fraction_difference(
            activity
        )

    def flux_from_mass_fraction(self, w_salt: float) -> float:
        """Evaporation flux from the liquid NaCl *mass* fraction (Prat's C)."""
        return self.flux_from_activity(water_activity_from_mass_fraction(w_salt))

    def flux_from_mole_fraction(self, x_salt: float) -> float:
        """Evaporation flux from the liquid NaCl *mole* fraction (Richards x)."""
        return self.flux_from_activity(water_activity_from_mole_fraction(x_salt))

    # -- diagnostics -------------------------------------------------------

    def equivalent_thickness_for_flux(
        self, activity: float, target_flux: float
    ) -> float:
        """Report the `delta` that would reproduce `target_flux` at `activity`.

        DIAGNOSTIC ONLY. This is *not* used to set `delta`; P2-A item 2 forbids
        calibrating the evaporation model against a prescribed flux. It exists
        so that the distance between the declared physical film and a legacy
        prescribed rate can be quantified and reported.
        """
        if not math.isfinite(target_flux) or target_flux == 0.0:
            raise ValueError("equivalent thickness requires a non-zero target flux")
        delta_x = self.driving_mass_fraction_difference(activity)
        return self.reference_gas_density * self.diffusivity * delta_x / target_flux
