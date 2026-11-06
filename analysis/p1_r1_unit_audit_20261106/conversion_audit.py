#!/usr/bin/env python3
"""Deterministic Case 36877 unit-conversion checks for the P1-R1/R2 audit."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[2]
CASE_INPUT = (
    REPO_ROOT
    / "analysis"
    / "matched_pair_36877_20261102"
    / "legacy_nx28"
    / "sandy--boreal--s1354--f1080--global_default_L0p5--native.input"
)
REQUIRED_KEYS = {
    "MolarMassH2O",
    "MolarMassNaCl",
    "PureWaterDensity",
    "DensityLinearCoefficient",
    "DensityQuadraticCoefficient",
    "EvaporationRate",
    "SaltMoleFraction",
}
EXPECTED_CASE_PARAMETERS = {
    "MolarMassH2O": "0.01801528",
    "MolarMassNaCl": "0.05844",
    "PureWaterDensity": "998.203",
    "DensityLinearCoefficient": "0.670804",
    "DensityQuadraticCoefficient": "0.373854",
    "EvaporationRate": "3.065544359023483e-09",
    "SaltMoleFraction": "0.0036",
}


def read_case_parameters(path: Path = CASE_INPUT) -> tuple[dict[str, str], str]:
    """Read required values from the frozen input and reject rounded substitutes."""
    raw = path.read_bytes()
    values: dict[str, str] = {}
    for line in raw.decode("utf-8").splitlines():
        payload = line.split("#", 1)[0].strip()
        if "=" not in payload:
            continue
        key, value = (part.strip() for part in payload.split("=", 1))
        if key in REQUIRED_KEYS:
            values[key] = value
    missing = REQUIRED_KEYS - values.keys()
    if missing:
        raise ValueError(f"missing frozen Case 36877 parameters: {sorted(missing)}")
    inconsistent = {
        key: {"expected": expected, "actual": values[key]}
        for key, expected in EXPECTED_CASE_PARAMETERS.items()
        if values[key] != expected
    }
    if inconsistent:
        raise ValueError(f"inconsistent Case 36877 constants: {inconsistent}")
    return values, hashlib.sha256(raw).hexdigest()


CASE_PARAMETERS, CASE_INPUT_SHA256 = read_case_parameters()
M_H2O = float(CASE_PARAMETERS["MolarMassH2O"])
M_NACL = float(CASE_PARAMETERS["MolarMassNaCl"])
RHO_WATER = float(CASE_PARAMETERS["PureWaterDensity"])
DENSITY_LINEAR = float(CASE_PARAMETERS["DensityLinearCoefficient"])
DENSITY_QUADRATIC = float(CASE_PARAMETERS["DensityQuadraticCoefficient"])
CASE_36877_EVAPORATION = float(CASE_PARAMETERS["EvaporationRate"])
CASE_36877_INITIAL_X = float(CASE_PARAMETERS["SaltMoleFraction"])


def properties(x_nacl: float) -> dict[str, float]:
    mbar = (1.0 - x_nacl) * M_H2O + x_nacl * M_NACL
    w_nacl = x_nacl * M_NACL / mbar
    # Match the frozen C++ law exactly: the salinity increment is referenced
    # to 1000 kg/m3, not multiplied by PureWaterDensity.
    rho = RHO_WATER + 1000.0 * w_nacl * (
        DENSITY_LINEAR + DENSITY_QUADRATIC * w_nacl
    )
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
    dc_dx = central_derivative(
        lambda value: properties(value)["c_mol_per_m3_liquid"], x_nacl
    )
    drho_dx = central_derivative(
        lambda value: properties(value)["rho_kg_per_m3_liquid"], x_nacl
    )
    return {
        "x_NaCl": x_nacl,
        "dc_dx_mol_per_m3": dc_dx,
        "extra_molar_salt_storage_fraction_x_dc_dx_over_c": (
            x_nacl * dc_dx / p["c_mol_per_m3_liquid"]
        ),
        "drho_dx_kg_per_m3": drho_dx,
        "density_change_fraction_per_unit_x_drho_dx_over_rho": (
            x_nacl * drho_dx / p["rho_kg_per_m3_liquid"]
        ),
    }


def build_result() -> dict[str, object]:
    states = {
        "case36877_initial": CASE_36877_INITIAL_X,
        "representative_p1_top_state": 0.008018127628,
        "high_concentration_audit_state": 0.1,
    }
    return {
        "purpose": "P1-R1/R2 deterministic conversion audit; not a solver validation",
        "source": {
            "frozen_case36877_input": CASE_INPUT.relative_to(REPO_ROOT).as_posix(),
            "sha256": CASE_INPUT_SHA256,
            "parameters_as_read": CASE_PARAMETERS,
        },
        "constants": {
            "M_H2O_kg_per_mol": M_H2O,
            "M_NaCl_kg_per_mol": M_NACL,
            "rho_reference_kg_per_m3": RHO_WATER,
            "density_linear_coefficient": DENSITY_LINEAR,
            "density_quadratic_coefficient": DENSITY_QUADRATIC,
        },
        "states": {name: properties(value) for name, value in states.items()},
        "non_equivalence": {
            name: non_equivalence_metrics(value) for name, value in states.items()
        },
        "case36877_evaporation_if_interpreted_as_pure_water_volume_demand": {
            "e_m_per_s": CASE_36877_EVAPORATION,
            "H2O_mass_flux_kg_per_m2_s": RHO_WATER * CASE_36877_EVAPORATION,
            "H2O_molar_flux_mol_per_m2_s": (
                RHO_WATER * CASE_36877_EVAPORATION / M_H2O
            ),
            "assumption": (
                "E is the equivalent pure-water volume demand from the frozen "
                "Case 36877 input, referenced to its PureWaterDensity"
            ),
        },
        "checks": {
            "frozen_parameter_identity": "PASS",
            "mass_molar_component_identities": "PASS",
            "total_mass_closure": "PASS",
        },
    }


def main() -> None:
    output = Path(__file__).with_name("conversion_audit.json")
    output.write_text(
        json.dumps(build_result(), indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(output)


if __name__ == "__main__":
    main()
