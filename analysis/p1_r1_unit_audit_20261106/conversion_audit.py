#!/usr/bin/env python3
"""Deterministic unit-conversion checks for the P1-R1 Richards audit."""

from __future__ import annotations

import json
from pathlib import Path


M_H2O = 0.018
M_NACL = 0.0584428
RHO_WATER = 1000.0
CASE_36877_EVAPORATION = 3.065544359e-9


def properties(x_nacl: float) -> dict[str, float]:
    mbar = (1.0 - x_nacl) * M_H2O + x_nacl * M_NACL
    w_nacl = x_nacl * M_NACL / mbar
    rho = RHO_WATER + RHO_WATER * w_nacl * (0.6708 + 0.3739 * w_nacl)
    molar_density = rho / mbar
    salt_mass_from_mass_basis = rho * w_nacl
    salt_mass_from_molar_basis = molar_density * x_nacl * M_NACL
    water_mass_from_mass_basis = rho * (1.0 - w_nacl)
    water_mass_from_molar_basis = molar_density * (1.0 - x_nacl) * M_H2O

    assert abs(salt_mass_from_mass_basis - salt_mass_from_molar_basis) < 1e-10
    assert abs(water_mass_from_mass_basis - water_mass_from_molar_basis) < 1e-10
    assert abs(rho - (salt_mass_from_mass_basis + water_mass_from_mass_basis)) < 1e-10

    return {
        "x_NaCl": x_nacl,
        "Mbar_kg_per_mol": mbar,
        "w_NaCl": w_nacl,
        "rho_kg_per_m3_liquid": rho,
        "c_mol_per_m3_liquid": molar_density,
        "salt_mass_concentration_kg_per_m3_liquid": salt_mass_from_mass_basis,
        "water_mass_concentration_kg_per_m3_liquid": water_mass_from_mass_basis,
    }


def central_derivative(func, x: float, h: float = 1e-7) -> float:
    return (func(x + h) - func(x - h)) / (2.0 * h)


def non_equivalence_metrics(x_nacl: float) -> dict[str, float]:
    p = properties(x_nacl)
    dc_dx = central_derivative(lambda value: properties(value)["c_mol_per_m3_liquid"], x_nacl)
    drho_dx = central_derivative(lambda value: properties(value)["rho_kg_per_m3_liquid"], x_nacl)
    return {
        "x_NaCl": x_nacl,
        "dc_dx_mol_per_m3": dc_dx,
        "extra_molar_salt_storage_fraction_x_dc_dx_over_c": x_nacl * dc_dx / p["c_mol_per_m3_liquid"],
        "drho_dx_kg_per_m3": drho_dx,
        "density_change_fraction_per_unit_x_drho_dx_over_rho": x_nacl * drho_dx / p["rho_kg_per_m3_liquid"],
    }


def main() -> None:
    states = {
        "case36877_initial": 0.0036,
        "representative_p1_top_state": 0.008018127628,
        "high_concentration_audit_state": 0.1,
    }
    result = {
        "purpose": "P1-R1 deterministic conversion audit; not a solver validation",
        "constants": {
            "M_H2O_kg_per_mol": M_H2O,
            "M_NaCl_kg_per_mol": M_NACL,
            "rho_reference_kg_per_m3": RHO_WATER,
        },
        "states": {name: properties(value) for name, value in states.items()},
        "non_equivalence": {
            name: non_equivalence_metrics(value) for name, value in states.items()
        },
        "case36877_evaporation_if_interpreted_as_pure_water_volume_demand": {
            "e_m_per_s": CASE_36877_EVAPORATION,
            "H2O_mass_flux_kg_per_m2_s": RHO_WATER * CASE_36877_EVAPORATION,
            "H2O_molar_flux_mol_per_m2_s": RHO_WATER * CASE_36877_EVAPORATION / M_H2O,
            "assumption": "e is an equivalent pure-water volume demand referenced to 1000 kg/m3",
        },
        "checks": {
            "mass_molar_component_identities": "PASS",
            "total_mass_closure": "PASS",
        },
    }
    output = Path(__file__).with_name("conversion_audit.json")
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()
