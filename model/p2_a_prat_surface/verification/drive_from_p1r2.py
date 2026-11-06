#!/usr/bin/env python3
"""Drive the P2-A one-way chain from the frozen P1-R2 surface salinity record.

P2-A item 7 requires the chain

    soil surface salinity -> water activity -> evaporation flux
                          -> Prat growth velocity -> crust height

to be implemented and verified as a **one-way diagnostic**. This script does
exactly that, reading the committed P1-R2 `surface_events.csv` records without
touching them.

Declared inconsistency (stated, not hidden)
-------------------------------------------
The P1-R2 trajectories were produced with a *fixed prescribed* evaporation
flux (Mode A: ``rho_w,ref * E`` with ``E = 3.065544359023483e-9 m/s``). This
module computes a *salinity-dependent* flux from the gas film. Feeding the
P1-R2 salinity history into the film model therefore mixes two different
atmospheric boundary conditions. That is acceptable for a one-way diagnostic
and is the reason the chain is labelled one-way; it is NOT a coupled result,
and no growth or onset claim may be drawn from it.

Second declared limitation
--------------------------
Within the committed P1-R2 window the soil surface NaCl mass fraction never
reaches halite equilibrium, so the declared seed event never fires and the
crust stays absent. To exercise the remaining links of the chain, a clearly
labelled **synthetic** concentration continuation is run as a separate case.

Usage::

    python drive_from_p1r2.py [output_directory]
"""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path

from reference_config import (
    PARAMS,
    reference_crust,
    reference_film,
    reference_patch,
)
from surface import coupling as cp
from surface import prat_crust as pc
from surface import thermo
from surface import water_activity as wa

C_EQ = PARAMS["equilibrium_mass_fraction"]
MODULE_ROOT = Path(__file__).resolve().parent.parent
REPO_ROOT = MODULE_ROOT.parent.parent

P1R2_RUNS = {
    "V5_mass_mode_A_1d": (
        "model/p1_r2_conservative_richards/results/stage_c_20261106/"
        "runs/V5_mass_mode_A_1d/surface_events.csv"
    ),
    "V5_mass_mode_B_1d": (
        "model/p1_r2_conservative_richards/results/stage_c_20261106/"
        "runs/V5_mass_mode_B_1d/surface_events.csv"
    ),
    "V8_case36877_R24_2d": (
        "model/p1_r2_conservative_richards/results/stage_d_20261106/"
        "run/surface_events.csv"
    ),
}

LEGACY_EVAPORATION_RATE = 3.065544359023483e-09  # m/s, frozen Case 36877 input
LEGACY_REFERENCE_DENSITY = 998.203               # kg/m3, Mode A reference
LEGACY_FLUX = LEGACY_REFERENCE_DENSITY * LEGACY_EVAPORATION_RATE

#: Thinning of the committed trajectory CSVs only; integration is unthinned.
OUTPUT_STRIDE = 20


def read_surface_events(path: Path) -> list[dict]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def drive_one_way(records: list[dict], case: str) -> tuple[list[dict], dict]:
    """Evaluate the chain on a frozen soil-salinity trajectory."""
    film = reference_film()
    crust = reference_crust()
    rows: list[dict] = []
    eligible_time = None
    peak_mass_fraction = 0.0
    for index, record in enumerate(records):
        t = float(record["time_s"])
        x_mean = float(record["x_top_area_mean"])
        x_max = float(record["x_top_max"])
        w_mean = thermo.mass_fraction_from_mole_fraction(x_mean)
        w_max = thermo.mass_fraction_from_mole_fraction(x_max)
        m = thermo.molality_from_mass_fraction(w_mean)
        a_w = wa.water_activity(m)
        j = film.flux_from_activity(a_w)
        peak_mass_fraction = max(peak_mass_fraction, w_max)
        if w_max >= C_EQ and eligible_time is None:
            eligible_time = t
        if index % OUTPUT_STRIDE and index != len(records) - 1:
            continue
        rows.append(
            {
                "case": case,
                "time_s": t,
                "x_top_area_mean": x_mean,
                "x_top_max": x_max,
                "w_top_area_mean": w_mean,
                "molality_mol_per_kg": m,
                "water_activity": a_w,
                "interfacial_vapour_mole_fraction": (
                    film.interfacial_vapour_mole_fraction(a_w)
                ),
                "driving_mass_fraction_difference": (
                    film.driving_mass_fraction_difference(a_w)
                ),
                "film_evaporation_flux_kg_m2_s": j,
                "p1r2_prescribed_flux_kg_m2_s": LEGACY_FLUX,
                "film_over_prescribed": j / LEGACY_FLUX,
                "solubility_eligible": bool(w_max >= C_EQ),
                "crust_active": False,
                "crust_height_m": 0.0,
                "growth_velocity_m_s": 0.0,
            }
        )
    first, last = rows[0], rows[-1]
    summary = {
        "case": case,
        "records_read": len(records),
        "records_written": len(rows),
        "output_stride": OUTPUT_STRIDE,
        "time_window_s": [first["time_s"], last["time_s"]],
        "soil_mass_fraction_range": [
            first["w_top_area_mean"],
            last["w_top_area_mean"],
        ],
        "water_activity_range": [first["water_activity"], last["water_activity"]],
        "film_flux_range_kg_m2_s": [
            first["film_evaporation_flux_kg_m2_s"],
            last["film_evaporation_flux_kg_m2_s"],
        ],
        "film_flux_relative_decline": (
            first["film_evaporation_flux_kg_m2_s"]
            - last["film_evaporation_flux_kg_m2_s"]
        )
        / first["film_evaporation_flux_kg_m2_s"],
        "equilibrium_mass_fraction": C_EQ,
        "max_soil_mass_fraction": peak_mass_fraction,
        "fraction_of_equilibrium_reached": peak_mass_fraction / C_EQ,
        "seed_event_fired": eligible_time is not None,
        "eligibility_time_s": eligible_time,
        "crust_height_m": 0.0,
        "note": (
            "One-way diagnostic. The salinity trajectory was produced under a "
            "fixed prescribed evaporation flux; the film flux reported here is "
            "salinity dependent and was not fed back. No onset or growth claim."
        ),
    }
    return rows, summary


def synthetic_continuation(
    interface: cp.EvaporationInterface,
    law: pc.GrowthLaw,
    ramp_seconds: float = 1.0e5,
    hold_seconds: float = 4.0e5,
    dt: float = 250.0,
) -> tuple[list[dict], dict]:
    """Declared synthetic soil-salinity forcing that reaches equilibrium.

    This is NOT a P1-R2 result. The soil-top mass fraction is prescribed to
    rise linearly from the frozen initial Case 36877 value to halite
    equilibrium over ``ramp_seconds`` and is then held, purely so that the
    seed event, the growth law, the inventories and the base fluxes are all
    exercised on a single trajectory.
    """
    w_start = thermo.mass_fraction_from_mole_fraction(0.0036)
    patch = reference_patch(interface=interface, law=law)
    rows: list[dict] = []
    initial = patch.state.copy()
    t = 0.0
    index = 0
    total = ramp_seconds + hold_seconds
    while t < total - 1.0e-9:
        step = min(dt, total - t)
        if t < ramp_seconds:
            w = w_start + (C_EQ - w_start) * (t / ramp_seconds)
        else:
            w = C_EQ
        activated = patch.consider_activation(w, t)
        record = patch.step(w, step)
        # Integrated at full resolution; the committed CSV is thinned, and the
        # activation step is always kept.
        if index % OUTPUT_STRIDE == 0 or activated:
            row = {
                "case": f"synthetic|{interface.value}|{law.value}",
                "output_stride": OUTPUT_STRIDE,
            }
            row.update(record)
            rows.append(row)
        index += 1
        t += step
    rows.append(
        {
            "case": f"synthetic|{interface.value}|{law.value}",
            "output_stride": OUTPUT_STRIDE,
            **patch.state.history[-1],
        }
    )
    budget = patch.budget(initial)
    last = patch.state.history[-1]
    summary = {
        "case": f"synthetic|{interface.value}|{law.value}",
        "declared_synthetic": True,
        "ramp_seconds": ramp_seconds,
        "hold_seconds": hold_seconds,
        "dt_s": dt,
        "activation_time_s": patch.state.activation_time,
        "eligibility_time_s": patch.state.eligibility_time,
        "seed_height_m": PARAMS["seed_height"],
        "final_height_m": patch.state.height,
        "height_growth_factor": patch.state.height / PARAMS["seed_height"],
        "final_evaporation_flux_kg_m2_s": last["evaporation_flux_kg_m2_s"],
        "final_growth_velocity_m_s": last["growth_velocity_m_s"],
        "final_salt_flux_soil_kg_m2_s": last["salt_flux_soil_kg_m2_s"],
        "final_water_flux_soil_kg_m2_s": last["water_flux_soil_kg_m2_s"],
        "final_top_supersaturation": last["top_supersaturation"],
        "final_damkohler": last["damkohler"],
        "salt_residual_relative": budget["salt_residual_relative"],
        "water_residual_relative": budget["water_residual_relative"],
        "dissolved_storage_gap_relative": abs(
            last["dissolved_storage_gap_kg_m2"]
        )
        / abs(last["dissolved_salt_kg_m2"]),
        "constant_porosity_defensible": last["validity_constant_porosity_defensible"],
        "large_da_valid": last["validity_large_da_valid"],
        "pe_much_smaller_than_da": last["validity_pe_much_smaller_than_da"],
        "internal_over_top_precipitation": last[
            "validity_internal_over_top_precipitation"
        ],
    }
    return rows, summary


def write_csv(path: Path, rows: list[dict]) -> None:
    if not rows:
        raise ValueError(f"refusing to write an empty table: {path}")
    fieldnames: list[str] = []
    for row in rows:
        for key in row:
            if key not in fieldnames:
                fieldnames.append(key)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, restval="")
        writer.writeheader()
        writer.writerows(rows)


def main(argv: list[str]) -> int:
    out = (
        Path(argv[1])
        if len(argv) > 1
        else MODULE_ROOT / "results" / "p1r2_driven_chain"
    )
    out.mkdir(parents=True, exist_ok=True)

    driven_rows: list[dict] = []
    summaries: list[dict] = []
    missing: list[str] = []
    for case, relative in P1R2_RUNS.items():
        path = REPO_ROOT / relative
        if not path.is_file():
            missing.append(relative)
            continue
        rows, summary = drive_one_way(read_surface_events(path), case)
        summary["source"] = relative
        driven_rows.extend(rows)
        summaries.append(summary)
    if not driven_rows:
        raise SystemExit(f"no P1-R2 surface records found; looked for {missing}")
    write_csv(out / "p1r2_driven_chain.csv", driven_rows)

    synthetic_rows: list[dict] = []
    synthetic_summaries: list[dict] = []
    for interface in cp.EvaporationInterface:
        for law in pc.GrowthLaw:
            rows, summary = synthetic_continuation(interface, law)
            synthetic_rows.extend(rows)
            synthetic_summaries.append(summary)
    write_csv(out / "synthetic_continuation.csv", synthetic_rows)

    payload = {
        "p1r2_driven": summaries,
        "missing_sources": missing,
        "synthetic_continuation": synthetic_summaries,
        "declared_limitations": [
            "The P1-R2 salinity histories were produced with a fixed prescribed "
            "evaporation flux; the film flux computed here was not fed back.",
            "No crust forms on any committed P1-R2 trajectory because the soil "
            "surface never reaches halite equilibrium within the frozen window.",
            "The synthetic continuation prescribes the soil concentration and is "
            "not a solution of the Richards-salt problem.",
        ],
    }
    with (out / "summary.json").open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2)
        handle.write("\n")

    for summary in summaries:
        print(
            f"{summary['case']}: w {summary['soil_mass_fraction_range'][0]:.6f} -> "
            f"{summary['soil_mass_fraction_range'][1]:.6f} "
            f"({100 * summary['fraction_of_equilibrium_reached']:.2f}% of C_eq), "
            f"film flux decline {100 * summary['film_flux_relative_decline']:.2f}%, "
            f"seed fired: {summary['seed_event_fired']}"
        )
    for summary in synthetic_summaries:
        print(
            f"{summary['case']}: h {PARAMS['seed_height']:.3e} -> "
            f"{summary['final_height_m']:.6e} m, "
            f"salt budget {summary['salt_residual_relative']:.2e}, "
            f"water budget {summary['water_residual_relative']:.2e}"
        )
    print(f"wrote results to {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
