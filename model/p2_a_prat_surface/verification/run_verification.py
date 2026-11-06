#!/usr/bin/env python3
"""Generate the P2-A numerical evidence under the declared reference inputs.

Writes CSV tables and a JSON summary into ``results/<stamp>/``. Nothing here
calibrates or tunes: every parameter comes from
``config/acceptance_criteria.json``.

Usage::

    python run_verification.py [output_directory]
"""

from __future__ import annotations

import csv
import json
import math
import platform
import sys
from datetime import date
from pathlib import Path

from reference_config import (
    CRITERIA,
    CRITERIA_PATH,
    PARAMS,
    reference_crust,
    reference_film,
    reference_patch,
    reference_seed,
)
from surface import coupling as cp
from surface import prat_crust as pc
from surface import thermo
from surface import water_activity as wa

C_EQ = PARAMS["equilibrium_mass_fraction"]
MODULE_ROOT = Path(__file__).resolve().parent.parent


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"refusing to write an empty table: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)


# --- V-A: water activity versus NaCl concentration ---------------------------


def table_water_activity() -> list[dict]:
    film = reference_film()
    rows = []
    n = 121
    for i in range(n):
        x = 0.1 * i / (n - 1)
        w = thermo.mass_fraction_from_mole_fraction(x)
        m = thermo.molality_from_mole_fraction(x)
        if m > wa.MOLALITY_LIMIT:
            continue
        a_w = wa.water_activity(m)
        rows.append(
            {
                "x_nacl": x,
                "w_nacl": w,
                "molality_mol_per_kg": m,
                "osmotic_coefficient": wa.osmotic_coefficient(m),
                "water_activity": a_w,
                "liquid_density_kg_m3": thermo.BrineDensityLaw().density(w),
                "vapour_mole_fraction_interface": (
                    film.interfacial_vapour_mole_fraction(a_w)
                ),
                "driving_mass_fraction_difference": (
                    film.driving_mass_fraction_difference(a_w)
                ),
                "evaporation_flux_kg_m2_s": film.flux_from_activity(a_w),
                "flux_relative_to_fresh_water": (
                    film.flux_from_activity(a_w) / film.flux_from_activity(1.0)
                ),
            }
        )
    return rows


# --- V-C: growth-law regimes --------------------------------------------------


def table_growth_regimes() -> list[dict]:
    crust = reference_crust()
    film = reference_film()
    j = film.flux_from_mass_fraction(C_EQ)
    da_star = crust.damkohler_star
    rows = []
    for i in range(61):
        h = 10.0 ** (-6.0 + 5.0 * i / 60.0)
        da = crust.damkohler(h)
        f_finite = pc.efficiency(da, da_star)
        f_large = pc.efficiency_large_damkohler(da_star)
        flags = pc.validity_flags(crust, h, j)
        fluxes = pc.interface_fluxes(crust, h, j)
        rows.append(
            {
                "height_m": h,
                "damkohler": da,
                "damkohler_top": crust.damkohler_top(h),
                "damkohler_star": da_star,
                "peclet": flags["peclet"],
                "pe_over_da": flags["pe_over_da"],
                "efficiency_finite_da": f_finite,
                "efficiency_large_da": f_large,
                "efficiency_relative_difference": abs(f_finite - f_large) / f_large,
                "growth_velocity_finite_da_m_s": pc.growth_velocity(crust, h, j),
                "growth_velocity_large_da_m_s": pc.growth_velocity(
                    crust, h, j, pc.GrowthLaw.LARGE_DA
                ),
                "top_supersaturation": fluxes.top_supersaturation,
                "top_mass_fraction": fluxes.top_mass_fraction,
                "internal_reaction_kg_m2_s": fluxes.internal_reaction,
                "top_skeleton_kg_m2_s": fluxes.top_skeleton_rate,
                "internal_over_top": flags["internal_over_top_precipitation"],
                "inverse_sqrt_damkohler_star": 1.0 / math.sqrt(da_star),
                "salt_flux_soil_kg_m2_s": fluxes.salt_flux,
                "water_flux_soil_kg_m2_s": fluxes.water_flux,
                "route_difference_kg_m2_s": fluxes.route_difference,
                "route_difference_identity_kg_m2_s": fluxes.route_difference_identity,
                "porosity_consumption_timescale_s": flags[
                    "porosity_consumption_timescale_s"
                ],
                "growth_timescale_s": flags["growth_timescale_s"],
                "timescale_separation": flags["timescale_separation"],
                "large_da_valid": flags["large_da_valid"],
                "pe_much_smaller_than_da": flags["pe_much_smaller_than_da"],
                "constant_porosity_defensible": flags["constant_porosity_defensible"],
            }
        )
    return rows


# --- V-D/V-E: height evolution and budgets -----------------------------------


def run_trajectory(interface, law, total=2.0e5, dt=200.0, scheme="heun") -> dict:
    patch = reference_patch(interface=interface, law=law)
    patch.consider_activation(C_EQ, 0.0)
    initial = patch.state.copy()
    elapsed = 0.0
    while elapsed < total - 1.0e-9:
        step = min(dt, total - elapsed)
        patch.step(C_EQ, step, scheme=scheme)
        elapsed += step
    return {"patch": patch, "budget": patch.budget(initial)}


#: The trajectories are integrated at full resolution; only the committed CSV
#: is thinned, so budgets and convergence orders are unaffected.
OUTPUT_STRIDE = 20


def table_trajectories() -> tuple[list[dict], list[dict]]:
    rows: list[dict] = []
    budgets: list[dict] = []
    for interface in cp.EvaporationInterface:
        for law in pc.GrowthLaw:
            result = run_trajectory(interface, law)
            patch = result["patch"]
            tag = f"{interface.value}|{law.value}"
            history = patch.state.history
            kept = set(range(0, len(history), OUTPUT_STRIDE)) | {len(history) - 1}
            for index, record in enumerate(history):
                if index not in kept:
                    continue
                row = {"case": tag, "output_stride": OUTPUT_STRIDE}
                row.update(record)
                rows.append(row)
            budget = {"case": tag}
            budget.update(result["budget"])
            last = patch.state.history[-1]
            budget["final_height_m"] = last["crust_height_m"]
            budget["final_evaporation_flux_kg_m2_s"] = last[
                "evaporation_flux_kg_m2_s"
            ]
            budget["final_growth_velocity_m_s"] = last["growth_velocity_m_s"]
            budget["dissolved_storage_gap_relative"] = abs(
                last["dissolved_storage_gap_kg_m2"]
            ) / abs(last["dissolved_salt_kg_m2"])
            budgets.append(budget)
    return rows, budgets


def table_time_step_convergence() -> list[dict]:
    rows = []
    for scheme in ("euler", "heun"):
        reference = run_trajectory(
            cp.EvaporationInterface.SOIL_TOP,
            pc.GrowthLaw.FINITE_DA,
            dt=12.5,
            scheme="heun",
        )["patch"].state.height
        previous_error = None
        for dt in (3200.0, 1600.0, 800.0, 400.0, 200.0):
            height = run_trajectory(
                cp.EvaporationInterface.SOIL_TOP,
                pc.GrowthLaw.FINITE_DA,
                dt=dt,
                scheme=scheme,
            )["patch"].state.height
            error = abs(height - reference)
            order = (
                math.log2(previous_error / error)
                if previous_error is not None and error > 0.0
                else float("nan")
            )
            rows.append(
                {
                    "scheme": scheme,
                    "dt_s": dt,
                    "final_height_m": height,
                    "reference_height_m": reference,
                    "absolute_error_m": error,
                    "observed_order": order,
                }
            )
            previous_error = error
    return rows


# --- item 9: interface comparison --------------------------------------------


def table_interface_comparison() -> list[dict]:
    film = reference_film()
    rows = []
    for exponent in range(-6, -1):
        for sub in (0.0, 0.5):
            h = 10.0 ** (exponent + sub)
            for a_v in (1.0e3, 1.0e4, 1.0e5, 1.0e6, 1.0e7):
                for k_r in (1.0e-9, 1.0e-8, 1.0e-7, 1.0e-6):
                    crust = reference_crust(
                        specific_area=a_v, reaction_coefficient=k_r
                    )
                    flux_a = film.flux_from_mass_fraction(C_EQ)
                    flags = pc.validity_flags(crust, h, flux_a)
                    row = {
                        "height_m": h,
                        "specific_area_1_m": a_v,
                        "reaction_coefficient_m_s": k_r,
                        "damkohler": flags["damkohler"],
                        "damkohler_star": flags["damkohler_star"],
                        "peclet": flags["peclet"],
                        "flux_option_a_kg_m2_s": flux_a,
                        "large_da_valid": flags["large_da_valid"],
                        "constant_porosity_defensible": flags[
                            "constant_porosity_defensible"
                        ],
                    }
                    try:
                        comparison = cp.compare_interfaces(film, crust, h, C_EQ)
                    except pc.CrustModelDomainExit:
                        # Option B left the validity range of the inherited
                        # activity correlation. Recorded, never extrapolated.
                        row.update(
                            {
                                "flux_option_b_kg_m2_s": float("nan"),
                                "relative_flux_difference": float("nan"),
                                "relative_growth_difference": float("nan"),
                                "top_supersaturation": float("nan"),
                                "top_mass_fraction": float("nan"),
                                "acceptable_at_1_percent": False,
                                "model_domain_exit": True,
                            }
                        )
                    else:
                        row.update(
                            {
                                "flux_option_b_kg_m2_s": comparison.flux_b,
                                "relative_flux_difference": (
                                    comparison.relative_flux_difference
                                ),
                                "relative_growth_difference": (
                                    comparison.relative_growth_difference
                                ),
                                "top_supersaturation": (
                                    comparison.top_supersaturation_b
                                ),
                                "top_mass_fraction": comparison.top_mass_fraction_b,
                                "acceptable_at_1_percent": comparison.acceptable,
                                "model_domain_exit": False,
                            }
                        )
                    rows.append(row)
    return rows


def flux_sensitivity_to_concentration(mass_fraction: float = C_EQ) -> float:
    """``S = -d ln j / dC`` [-], the local salinity sensitivity of the film.

    Option A replaces the crust-top mass fraction by the base value, i.e. it
    omits ``C*_top``. To first order the resulting relative evaporation-flux
    error is therefore ``S * C*_top``, which is the practical acceptance
    criterion for item 9.
    """
    film = reference_film()
    delta = 1.0e-7
    hi = film.flux_from_mass_fraction(mass_fraction + delta)
    lo = film.flux_from_mass_fraction(mass_fraction - delta)
    centre = film.flux_from_mass_fraction(mass_fraction)
    return -(hi - lo) / (2.0 * delta) / centre


# --- supply limitation --------------------------------------------------------


def table_supply_response() -> list[dict]:
    crust = reference_crust()
    film = reference_film()
    j = film.flux_from_mass_fraction(C_EQ)
    fluxes = pc.interface_fluxes(crust, 1.0e-3, j)
    rows = []
    for fraction in (2.0, 1.5, 1.0, 0.9, 0.5, 0.1):
        diag = pc.supply_limitation_diagnostic(
            fluxes,
            fraction * fluxes.salt_flux,
            fraction * fluxes.water_flux,
        )
        raised = False
        try:
            pc.supply_limitation_diagnostic(
                fluxes,
                fraction * fluxes.salt_flux,
                fraction * fluxes.water_flux,
                raise_on_shortfall=True,
            )
        except pc.CrustModelDomainExit:
            raised = True
        rows.append(
            {
                "available_fraction_of_demand": fraction,
                "salt_demand_kg_m2_s": diag.salt_demand,
                "salt_available_kg_m2_s": diag.salt_available,
                "salt_ratio": diag.salt_ratio,
                "salt_limited": diag.salt_limited,
                "water_limited": diag.water_limited,
                "model_domain_exit_raised": raised,
                "growth_velocity_unchanged_m_s": fluxes.growth_velocity,
            }
        )
    return rows


def table_connectivity_response() -> list[dict]:
    film = reference_film()
    j = film.flux_from_mass_fraction(C_EQ)
    rows = []
    for saturation in (1.0, 0.995, 0.99, 0.95, 0.8, 0.4):
        crust = reference_crust(crust_saturation=saturation)
        try:
            fluxes = pc.interface_fluxes(crust, 1.0e-3, j)
            rows.append(
                {
                    "crust_saturation": saturation,
                    "closure_active": True,
                    "salt_flux_soil_kg_m2_s": fluxes.salt_flux,
                    "water_flux_soil_kg_m2_s": fluxes.water_flux,
                    "growth_velocity_m_s": fluxes.growth_velocity,
                }
            )
        except pc.CrustModelDomainExit:
            rows.append(
                {
                    "crust_saturation": saturation,
                    "closure_active": False,
                    "salt_flux_soil_kg_m2_s": float("nan"),
                    "water_flux_soil_kg_m2_s": float("nan"),
                    "growth_velocity_m_s": float("nan"),
                }
            )
    return rows


# --- summary ------------------------------------------------------------------


def build_summary(tables: dict[str, list[dict]]) -> dict:
    film = reference_film()
    crust = reference_crust()
    j_eq = film.flux_from_mass_fraction(C_EQ)
    j_fresh = film.flux_from_mass_fraction(0.0)
    legacy_rate = 3.065544359023483e-09
    legacy_flux = 998.203 * legacy_rate

    interface_rows = tables["interface_comparison"]
    evaluated = [r for r in interface_rows if not r["model_domain_exit"]]
    exits = [r for r in interface_rows if r["model_domain_exit"]]
    acceptable = [r for r in evaluated if r["acceptable_at_1_percent"]]
    worst = max(evaluated, key=lambda r: r["relative_flux_difference"])

    regimes = tables["growth_regimes"]
    large_da_rows = [r for r in regimes if r["large_da_valid"]]
    porosity_ok = [r for r in regimes if r["constant_porosity_defensible"]]

    sensitivity = flux_sensitivity_to_concentration()
    convergence = tables["time_step_convergence"]
    orders = {
        scheme: [
            r["observed_order"]
            for r in convergence
            if r["scheme"] == scheme and not math.isnan(r["observed_order"])
        ]
        for scheme in ("euler", "heun")
    }

    return {
        "generated": date.today().isoformat(),
        "python": platform.python_version(),
        "criteria_file": str(CRITERIA_PATH.relative_to(MODULE_ROOT)),
        "declared_reference": PARAMS,
        "evaporation": {
            "flux_fresh_water_kg_m2_s": j_fresh,
            "flux_at_equilibrium_kg_m2_s": j_eq,
            "salinity_reduction_fraction": (j_fresh - j_eq) / j_fresh,
            "water_activity_at_equilibrium": wa.water_activity_from_mass_fraction(
                C_EQ
            ),
            "transfer_coefficient_kg_m2_s": film.transfer_coefficient,
            "vapour_diffusion_coefficient_m2_s": film.diffusivity,
            "saturation_pressure_pa": film.saturation_pressure,
            "legacy_p1_prescribed_flux_kg_m2_s": legacy_flux,
            "ratio_declared_film_to_legacy": j_eq / legacy_flux,
            "equivalent_film_thickness_for_legacy_flux_m": (
                film.equivalent_thickness_for_flux(
                    wa.water_activity_from_mass_fraction(C_EQ), legacy_flux
                )
            ),
            "note": (
                "The equivalent thickness is reported as a diagnostic only. "
                "P2-A item 2 forbids calibrating the film against a prescribed "
                "flux, and no calibration was performed."
            ),
        },
        "growth_law": {
            "damkohler_star": crust.damkohler_star,
            "efficiency_large_da": pc.efficiency_large_damkohler(
                crust.damkohler_star
            ),
            "internal_over_top_large_da_limit": 1.0
            / math.sqrt(crust.damkohler_star),
            "identity": (
                "In the large-Da limit the depth-integrated internal reaction "
                "divided by the top skeleton rate equals 1/sqrt(Da*), so "
                "f_inf = 1/(1 + internal/top). Derived in this work from the "
                "analytic profile; Prat (2026) does not state it."
            ),
            "smallest_height_with_da_above_10_m": (
                min(r["height_m"] for r in large_da_rows) if large_da_rows else None
            ),
            "heights_with_defensible_constant_porosity": len(porosity_ok),
            "heights_examined": len(regimes),
        },
        "interface_comparison": {
            "rows": len(interface_rows),
            "evaluated": len(evaluated),
            "model_domain_exits": len(exits),
            "acceptable_at_1_percent": len(acceptable),
            "max_relative_flux_difference": worst["relative_flux_difference"],
            "max_difference_case": {
                "height_m": worst["height_m"],
                "specific_area_1_m": worst["specific_area_1_m"],
                "reaction_coefficient_m_s": worst["reaction_coefficient_m_s"],
                "damkohler": worst["damkohler"],
            },
            "largest_acceptable_peclet": (
                max(abs(r["peclet"]) for r in acceptable) if acceptable else None
            ),
            "smallest_unacceptable_peclet": (
                min(
                    abs(r["peclet"])
                    for r in evaluated
                    if not r["acceptable_at_1_percent"]
                )
                if len(acceptable) < len(evaluated)
                else None
            ),
            "largest_acceptable_top_supersaturation": (
                max(r["top_supersaturation"] for r in acceptable)
                if acceptable
                else None
            ),
            "flux_sensitivity_to_concentration": sensitivity,
            "first_order_supersaturation_threshold": 0.01 / sensitivity,
            "first_order_prediction_max_relative_error": max(
                abs(
                    r["relative_flux_difference"]
                    - sensitivity * r["top_supersaturation"]
                )
                / max(r["relative_flux_difference"], 1.0e-30)
                for r in evaluated
                if r["relative_flux_difference"] > 1.0e-9
            ),
            "criterion": (
                "Option A omits the crust-top supersaturation C*_top, so to "
                "first order its relative evaporation-flux error is "
                "S * C*_top with S = -d ln j / dC evaluated at C_eq. "
                "At the declared film, S is reported above and the 1 per cent "
                "acceptance level corresponds to the quoted C*_top threshold. "
                "Peclet alone is NOT the criterion: the sweep contains "
                "acceptable rows at larger Pe than some unacceptable rows, "
                "because C*_top scales as Pe/sqrt(Da) for Da >> 1."
            ),
        },
        "budgets": tables["budgets"],
        "time_step_convergence": {
            scheme: {
                "orders": values,
                "min_order": min(values) if values else None,
            }
            for scheme, values in orders.items()
        },
        "supply_response": tables["supply_response"],
        "connectivity_response": tables["connectivity_response"],
        "limitations": CRITERIA["notes"],
    }


def main(argv: list[str]) -> int:
    out = Path(argv[1]) if len(argv) > 1 else MODULE_ROOT / "results" / "p2a_verification"
    out.mkdir(parents=True, exist_ok=True)

    trajectories, budgets = table_trajectories()
    tables = {
        "water_activity": table_water_activity(),
        "growth_regimes": table_growth_regimes(),
        "trajectories": trajectories,
        "budgets": budgets,
        "time_step_convergence": table_time_step_convergence(),
        "interface_comparison": table_interface_comparison(),
        "supply_response": table_supply_response(),
        "connectivity_response": table_connectivity_response(),
    }
    for name, rows in tables.items():
        write_csv(out / f"{name}.csv", rows)

    summary = build_summary(tables)
    with (out / "summary.json").open("w", encoding="utf-8") as handle:
        json.dump(summary, handle, indent=2, sort_keys=False)
        handle.write("\n")

    print(f"wrote {len(tables)} tables and summary.json to {out}")
    for name, rows in tables.items():
        print(f"  {name}.csv: {len(rows)} rows")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
