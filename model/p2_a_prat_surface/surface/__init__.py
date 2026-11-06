"""P2-A external porous salt surface subgrid (standalone).

This package is an *isolated* module. It does not import, modify or depend on
any frozen P1 source, and it performs no DuMuX call. It supplies:

* `thermo`           - documented thermophysical correlations (source-traced)
* `water_activity`   - NaCl water activity (Pitzer osmotic coefficient)
* `film_evaporation` - physical gas-film evaporation mass flux
* `prat_crust`       - Prat (2026) porous salt growth law + reduced crust
                       transport + moving-interface balances + inventories
* `coupling`         - one-way diagnostic chain and conservative base-flux
                       reconstruction, with evaporation-interface options A/B

No scipy/torch dependency. numpy is used only for array bookkeeping.
"""

__all__ = [
    "thermo",
    "water_activity",
    "film_evaporation",
    "prat_crust",
    "coupling",
]
