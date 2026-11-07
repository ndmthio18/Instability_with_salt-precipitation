#!/usr/bin/env python3
"""Emit reference vectors for the C++ port from the verified Python model.

No C++ toolchain exists in this environment, so `cpp/pratcrust.hh` cannot be
compiled or run here. These vectors are the gate that replaces that: they are
produced by the Python model that P2-A, P2-B and P2-C verified, and
`cpp/verify_vectors.cc` recomputes every one of them from the header and
reports the largest relative difference. Until that check passes on a machine
with a compiler, the C++ module is unverified.

The declared configuration below is compiled into the verifier as well, so the
two sides cannot drift apart silently.

Usage::

    python generate_cpp_vectors.py [output_path]
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import reference_config_p2d as rc
from surface import graded_crust as gc
from surface import thermo
from surface.water_activity import osmotic_coefficient, water_activity_from_mass_fraction

HERE = Path(__file__).resolve().parent
DEFAULT_OUTPUT = HERE.parent / "cpp" / "reference_vectors.txt"

#: Declared vector configuration. Mirrored verbatim in verify_vectors.cc.
SEED_HEIGHT = 1.0e-5
CELLS_PER_DELTA = 8.0
WINDOW_LENGTHS = 40.0
RATIO = 1.25
MAXIMUM_CELL_LENGTHS = 10.0
DT = 5.0
STEPS = 10


def build():
    params = rc.crust_parameters()
    film = rc.film()
    delta = params.reaction_length()
    state = gc.GradedCrustState(
        top_cell_size=min(delta / CELLS_PER_DELTA, SEED_HEIGHT / 8.0)
    )
    gc.activate(
        state,
        params,
        SEED_HEIGHT,
        film.flux_from_mass_fraction(rc.C_EQ),
        0.0,
        window=WINDOW_LENGTHS * delta,
        ratio=RATIO,
        maximum_cell=MAXIMUM_CELL_LENGTHS * delta,
    )
    return params, film, state


def vectors() -> list:
    params, film, state = build()
    c_eq = params.equilibrium_mass_fraction
    out = [
        ("constant_C_eq", c_eq),
        ("constant_rho", params.rho),
        ("constant_a_v0", params.area(params.deposition_porosity)),
        ("thermo_psat_293_15", thermo.saturation_pressure(293.15)),
        ("thermo_dv_293_15_1e5", thermo.vapour_diffusion_coefficient(293.15, 1.0e5)),
        ("thermo_gaswmf_0_02", thermo.gas_water_mass_fraction(0.02)),
        ("thermo_w_from_x_0_1", thermo.mass_fraction_from_mole_fraction(0.1)),
        ("thermo_x_from_w_Ceq", thermo.mole_fraction_from_mass_fraction(c_eq)),
        ("thermo_molality_Ceq", thermo.molality_from_mass_fraction(c_eq)),
        ("thermo_brine_density_Ceq", thermo.BrineDensityLaw().density(c_eq)),
        ("activity_osmotic_1", osmotic_coefficient(1.0)),
        ("activity_osmotic_6", osmotic_coefficient(6.0)),
        ("activity_from_w_Ceq", water_activity_from_mass_fraction(c_eq)),
        ("film_transfer_coefficient", film.transfer_coefficient),
        ("film_flux_fresh", film.flux_from_mass_fraction(0.0)),
        ("film_flux_Ceq", film.flux_from_mass_fraction(c_eq)),
        ("closure_area_0_43", params.area(0.43)),
        ("closure_area_0_25", params.area(0.25)),
        ("closure_diffusivity_0_43", params.diffusivity(0.43)),
        ("closure_diffusivity_0_25", params.diffusivity(0.25)),
        ("closure_reaction_length_0_43", params.reaction_length(0.43)),
        ("closure_damkohler_star_0_43", params.damkohler_star(0.43)),
        ("seed_cells", float(len(state.thicknesses))),
        ("seed_height", state.height),
        ("seed_top_cell", state.top_cell_size),
        ("seed_solid_salt", state.solid_salt_mass),
        ("seed_dissolved_salt", state.dissolved_salt),
        ("seed_water", state.water),
    ]

    flux = film.flux_from_mass_fraction(c_eq)
    for label, base in (("Ceq", None), ("0_99Ceq", 0.99 * c_eq),
                        ("1_01Ceq", 1.01 * c_eq)):
        solution = gc.solve_profile(state, params, flux, base)
        prefix = f"profile_{label}_"
        out += [
            (prefix + "s_top", solution.top_supersaturation),
            (prefix + "v_s", solution.growth_velocity),
            (prefix + "internal_reaction", solution.internal_reaction),
            (prefix + "dissolved_salt", solution.dissolved_salt),
            (prefix + "water", solution.water),
            (prefix + "rho_u_base", solution.solution_mass_flux_base),
            (prefix + "salt_flux_transport", solution.salt_flux_base_route),
            (prefix + "base_gradient", solution.base_gradient),
        ]

    # Ten steps at a halite-equilibrated base.
    evolving = state
    for _ in range(STEPS):
        record = gc.step(evolving, params, flux, DT)
    out += [
        ("step10_time", evolving.time),
        ("step10_height", evolving.height),
        ("step10_cells", float(len(evolving.thicknesses))),
        ("step10_minimum_porosity", record["minimum_porosity"]),
        ("step10_solid_salt", evolving.solid_salt_mass),
        ("step10_dissolved_salt", evolving.dissolved_salt),
        ("step10_water", evolving.water),
        ("step10_salt_flux_soil", record["salt_flux_soil_kg_m2_s"]),
        ("step10_water_flux_soil", record["water_flux_soil_kg_m2_s"]),
        ("step10_salt_flux_transport", record["salt_flux_base_route_kg_m2_s"]),
        ("step10_growth_velocity", record["growth_velocity_m_s"]),
        ("step10_internal_reaction", record["internal_reaction_kg_m2_s"]),
    ]
    return out


def main(argv: list) -> int:
    path = Path(argv[1]) if len(argv) > 1 else DEFAULT_OUTPUT
    path.parent.mkdir(parents=True, exist_ok=True)
    rows = vectors()
    lines = [
        "# P2-D reference vectors for the C++ port of the Prat crust module.",
        "# Produced by verification/generate_cpp_vectors.py from the Python",
        "# model verified in P2-A, P2-B and P2-C. cpp/verify_vectors.cc must",
        "# reproduce every value from cpp/pratcrust.hh.",
        "#",
        "# Declared configuration (mirrored in verify_vectors.cc):",
        f"#   seed height        {SEED_HEIGHT!r} m",
        f"#   cells per delta_r  {CELLS_PER_DELTA!r}",
        f"#   window             {WINDOW_LENGTHS!r} delta_r",
        f"#   grading ratio      {RATIO!r}",
        f"#   maximum cell       {MAXIMUM_CELL_LENGTHS!r} delta_r",
        f"#   time step          {DT!r} s, {STEPS!r} steps",
        "#",
        "# name value",
    ]
    for name, value in rows:
        lines.append(f"{name} {value!r}")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")
    print(f"wrote {len(rows)} reference vectors to {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
