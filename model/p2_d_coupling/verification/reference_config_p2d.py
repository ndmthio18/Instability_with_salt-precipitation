"""The declared P2-D configuration.

Soil parameters are transcribed from the P1-R2 Case 36877 verification inputs
(`model/p1_r2_conservative_richards/verification/inputs/V8_case36877_...input`
and `V5_mass_mode_A_1d.input`, which share them). Crust parameters follow
Prat (2026) equations (16)/(17) and the paper's own values, as recorded in
`docs/PRAT_PHYSICS_AUDIT.md` section 3.3. The gas film is the P2-A declared
film. Nothing here is fitted.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
P2D_ROOT = HERE.parent
REPO_ROOT = P2D_ROOT.parent.parent
P2A_ROOT = REPO_ROOT / "model" / "p2_a_prat_surface"
for path in (str(P2D_ROOT), str(P2A_ROOT), str(P2A_ROOT / "verification")):
    if path not in sys.path:
        sys.path.insert(0, path)

from column import soil  # noqa: E402
from surface import graded_crust as gc  # noqa: E402
from surface.film_evaporation import GasFilm  # noqa: E402
from surface.water_activity import MOLALITY_LIMIT  # noqa: E402

CRITERIA_PATH = HERE / "config" / "acceptance_criteria_p2d.json"
with CRITERIA_PATH.open("r", encoding="utf-8") as _handle:
    CRITERIA_P2D = json.load(_handle)
TOL = CRITERIA_P2D["tolerances"]
DECLARED = {k: v["value"] for k, v in CRITERIA_P2D["declared_parameters"].items()}

#: P1-R2 Case 36877 soil, byte-for-byte from its verification inputs.
SOIL = soil.SoilParameters(
    molar_mass_h2o=0.01801528,
    molar_mass_nacl=0.05844,
    pure_water_density=998.203,
    density_linear_coefficient=0.670804,
    density_quadratic_coefficient=0.373854,
    porosity=0.4057687912195745,
    permeability=1.261386079996293e-11,
    pc_alpha=0.0003121789225428109,
    pc_n=2.98602455003942,
    residual_saturation=0.1282679353072967,
    swe_l=-0.8395681012108449,
    viscosity=0.0011,
    salt_diffusivity=1.5e-09,
    gravity=9.800000000000001,
    enable_gravity=True,
    krw_mualem_exponent=True,
)

COLUMN_HEIGHT = 0.2
WATER_LEVEL = 0.6307756332876138
INITIAL_MOLE_FRACTION = 0.0036
SOLUBILITY_MOLE_FRACTION = 0.1          # frozen P1 Transport.SolubilityMoleFraction
C_EQ = SOIL.mass_fraction(SOLUBILITY_MOLE_FRACTION)

#: Prat (2026) values, audit section 3.3.
PRAT = {
    "k_r": 2.3e-3,
    "eps_0": 0.43,
    "d_sub": 25e-6,
    "d_c_over_d_sub": 0.14,
    "D_s": 1.3e-9,
    "beta": 2.0,
}
A_V0 = 36.0 * (1.0 - PRAT["eps_0"]) / PRAT["d_sub"]


def film(**overrides) -> GasFilm:
    """The P2-A declared gas film."""
    kwargs = dict(
        temperature=293.15,
        relative_humidity=0.4,
        thickness=0.0083,
        atmospheric_pressure=1.0e5,
    )
    kwargs.update(overrides)
    return GasFilm(**kwargs)


def crust_parameters(**overrides) -> gc.GradedCrustParameters:
    """Prat equation (17) area closure and the eps^(1/2) tortuosity closure.

    ``C_eq`` is the soil's own halite mass fraction so the soil and the crust
    agree on what saturation means; ``rho`` is the soil's brine density law at
    that mass fraction, for the same reason.
    """
    kwargs = dict(
        equilibrium_mass_fraction=C_EQ,
        deposition_porosity=PRAT["eps_0"],
        crystal_density=2165.0,
        specific_area=A_V0,
        reaction_coefficient=PRAT["k_r"],
        effective_salt_diffusivity=math.sqrt(PRAT["eps_0"]) * PRAT["D_s"],
        beta=PRAT["beta"],
        solution_density=SOIL.rho(SOLUBILITY_MOLE_FRACTION),
        vapour_density=0.0,
        clogging_porosity=0.01,
        connectivity_porosity=0.02,
        area_closure="salt_grow_b2",
        substrate_size=PRAT["d_sub"],
        diffusivity_closure="sqrt_eps",
        molecular_diffusivity=PRAT["D_s"],
    )
    kwargs.update(overrides)
    return gc.GradedCrustParameters(**kwargs)


def hydrostatic_state(
    cells: int | None = None,
    mole_fraction: float | None = None,
    height: float | None = None,
) -> soil.ColumnState:
    """P1-R2's own initial profile: ``p = p_n - (waterLevel + y) rho g``."""
    n = int(DECLARED["column_cells"]) if cells is None else cells
    h = COLUMN_HEIGHT if height is None else height
    x = INITIAL_MOLE_FRACTION if mole_fraction is None else mole_fraction
    state = soil.ColumnState(cells=n, height=h)
    rho = SOIL.rho(x)
    g = SOIL.gravitational_acceleration()
    for index in range(n):
        y = state.centre(index)
        state.pressure.append(
            soil.P_NON_WETTING_REFERENCE - (WATER_LEVEL + y) * rho * g
        )
        state.mole_fraction.append(x)
    return state


def solubility_stop(params, state):
    """Declared seeding trigger: the soil surface reaches halite saturation."""
    if params.mass_fraction(state.mole_fraction[-1]) >= C_EQ:
        return "soil_top_reached_solubility"
    return None


def no_flux(_params, _state):
    return 0.0, 0.0


def prescribed_flux(total: float, salt: float):
    def closure(_params, _state):
        return total, salt
    return closure


def film_flux(gas_film: GasFilm | None = None):
    """Water-activity-controlled evaporation (task item 4).

    The flux is evaluated from the **top cell's** NaCl mass fraction through
    the Pitzer water activity, so it falls as the surface concentrates. No
    rate is prescribed and no weight is normalised.
    """
    chosen = film() if gas_film is None else gas_film

    def closure(params: soil.SoilParameters, state: soil.ColumnState):
        w = params.mass_fraction(state.mole_fraction[-1])
        molality = w / ((1.0 - w) * params.molar_mass_nacl)
        if molality > MOLALITY_LIMIT:
            raise soil.ColumnDomainExit(
                f"molality {molality!r} exceeds the declared Pitzer limit "
                f"{MOLALITY_LIMIT!r}; the activity correlation is not "
                "extrapolated"
            )
        return chosen.flux_from_mass_fraction(w), 0.0

    return closure
