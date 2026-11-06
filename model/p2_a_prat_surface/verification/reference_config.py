"""The single declared P2-A reference configuration.

Every number here is taken from `config/acceptance_criteria.json`, which was
registered before the verification runs. Nothing in this file is fitted and
nothing is calibrated against an observed evaporation rate.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
MODULE_ROOT = HERE.parent
if str(MODULE_ROOT) not in sys.path:
    sys.path.insert(0, str(MODULE_ROOT))

from surface import prat_crust as pc  # noqa: E402
from surface.coupling import EvaporationInterface, SurfacePatch  # noqa: E402
from surface.film_evaporation import GasFilm  # noqa: E402

CRITERIA_PATH = HERE / "config" / "acceptance_criteria.json"


def load_criteria() -> dict:
    with CRITERIA_PATH.open("r", encoding="utf-8") as handle:
        return json.load(handle)


CRITERIA = load_criteria()
PARAMS = {k: v["value"] for k, v in CRITERIA["declared_reference_parameters"].items()}
TOL = CRITERIA["tolerances"]


def reference_crust(**overrides) -> pc.CrustParameters:
    kwargs = dict(
        equilibrium_mass_fraction=PARAMS["equilibrium_mass_fraction"],
        crust_porosity=PARAMS["crust_porosity"],
        crystal_density=PARAMS["crystal_density"],
        specific_area=PARAMS["specific_area"],
        reaction_coefficient=PARAMS["reaction_coefficient"],
        effective_salt_diffusivity=PARAMS["effective_salt_diffusivity"],
        beta=PARAMS["beta"],
        vapour_density=PARAMS["vapour_density"],
    )
    kwargs.update(overrides)
    return pc.CrustParameters(**kwargs)


def reference_film(**overrides) -> GasFilm:
    kwargs = dict(
        temperature=PARAMS["temperature"],
        relative_humidity=PARAMS["relative_humidity"],
        thickness=PARAMS["film_thickness"],
        atmospheric_pressure=PARAMS["atmospheric_pressure"],
    )
    kwargs.update(overrides)
    return GasFilm(**kwargs)


def reference_seed(**overrides) -> pc.SeedEvent:
    kwargs = dict(
        trigger_mass_fraction=PARAMS["equilibrium_mass_fraction"],
        seed_height=PARAMS["seed_height"],
        seed_porosity=PARAMS["crust_porosity"],
        declared_by="P2-A verification reference configuration",
    )
    kwargs.update(overrides)
    return pc.SeedEvent(**kwargs)


def reference_patch(
    interface: EvaporationInterface = EvaporationInterface.SOIL_TOP,
    law: pc.GrowthLaw = pc.GrowthLaw.FINITE_DA,
    activated: bool = False,
    height: float | None = None,
) -> SurfacePatch:
    patch = SurfacePatch(
        film=reference_film(),
        crust=reference_crust(),
        seed=reference_seed(),
        law=law,
        interface=interface,
    )
    if activated:
        patch.state.active = True
        patch.state.height = PARAMS["seed_height"] if height is None else height
        patch.state.solid_salt = (
            (1.0 - PARAMS["crust_porosity"])
            * PARAMS["crystal_density"]
            * patch.state.height
        )
    return patch
