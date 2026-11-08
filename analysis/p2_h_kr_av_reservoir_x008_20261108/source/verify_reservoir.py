"""Fail-closed checks for P2-H reservoir verification and the 50 h pilot gate."""
import csv
import json
import math
import os
import pathlib
import sys


root = pathlib.Path.cwd()
stage = sys.argv[1]
expected_salt_mole_fraction = float(os.environ.get("P2H_SALT_MOLE_FRACTION", "0.05"))
checks = []


def record(name, passed, detail):
    checks.append({"name": name, "passed": bool(passed), "detail": detail})


def rows(path):
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def finite(values):
    return all(math.isfinite(float(value)) for value in values)


def audit_case(path, expect_mode):
    outcome = json.loads((path / "outcome.json").read_text())
    record(f"{path.name}: exit", outcome["exit_code"] == 0, outcome)
    config = rows(path / "reservoir_config.csv")[0]
    record(f"{path.name}: boundary mode", config["bottom_mode"] == expect_mode, config)
    record(
        f"{path.name}: domain/water-level geometry",
        abs(float(config["domain_bottom_y_m"])) < 1e-15
        and abs(float(config["domain_top_y_m"]) - 0.2) < 1e-15
        and abs(float(config["water_level_parameter_m"]) - 0.01) < 1e-15
        and abs(float(config["water_table_y_m"]) + 0.01) < 1e-15,
        config,
    )
    record(
        f"{path.name}: initial/reservoir xNaCl",
        max(
            abs(float(config[key]) - expected_salt_mole_fraction)
            for key in (
                "configured_salt_mole_fraction",
                "reservoir_salt_mole_fraction",
                "initial_min_salt_mole_fraction",
                "initial_max_salt_mole_fraction",
            )
        )
        < 1e-14,
        config,
    )
    budget = rows(path / "open_boundary_budget.csv")
    record(f"{path.name}: finite boundary ledger", bool(budget) and finite(v for row in budget for v in row.values()), {"rows": len(budget)})
    reconstruction = []
    for row in budget[1:]:
        dt = float(row["dt_s"])
        if dt <= 0:
            continue
        liquid = (float(row["bottom_liquid_after_kg"]) - float(row["bottom_liquid_before_kg"])) / dt - float(row["cut_liquid_down_kg_s"])
        salt = (float(row["bottom_salt_after_kg"]) - float(row["bottom_salt_before_kg"])) / dt - float(row["cut_salt_down_kg_s"])
        reconstruction.append(max(abs(liquid - float(row["bottom_liquid_in_kg_s"])), abs(salt - float(row["bottom_salt_in_kg_s"]))))
    record(f"{path.name}: local face reconstruction", max(reconstruction, default=0.0) < 1e-16, {"max_abs_kg_s": max(reconstruction, default=0.0)})
    history = rows(path / "coupled_history.csv")
    s0 = abs(float(history[0]["soil_salt_kg"]))
    w0 = abs(float(history[0]["soil_water_kg"]))
    max_sr = max(abs(float(row["salt_residual_kg"])) for row in history)
    max_wr = max(abs(float(row["water_residual_kg"])) for row in history)
    record(f"{path.name}: salt balance", max_sr / max(s0, 1e-30) <= 1e-8, {"relative": max_sr / max(s0, 1e-30)})
    record(f"{path.name}: water balance", max_wr / max(w0, 1e-30) <= 1e-8, {"relative": max_wr / max(w0, 1e-30)})
    return history, budget


if stage == "verification":
    hydro_h, hydro_b = audit_case(root / stage / "hydrostatic_reservoir", "reservoir")
    evap_h, evap_b = audit_case(root / stage / "evaporating_reservoir", "reservoir")
    closed_h, closed_b = audit_case(root / stage / "no_flux_regression", "no_flux")
    hydro_changes = {
        key: abs(float(hydro_h[-1][key]) - float(hydro_h[0][key]))
        for key in ("soil_salt_kg", "soil_water_kg", "top_soil_mass_fraction", "top_soil_sw")
    }
    record("hydrostatic reservoir remains stationary", max(hydro_changes.values()) < 1e-11, hydro_changes)
    record(
        "hydrostatic reservoir reaction is zero",
        max(abs(float(row["bottom_liquid_in_kg_s"])) for row in hydro_b) < 1e-12
        and max(abs(float(row["bottom_salt_in_kg_s"])) for row in hydro_b) < 1e-12,
        {"last": hydro_b[-1]},
    )
    record(
        "no-flux ledger remains zero",
        max(abs(float(row["bottom_liquid_in_kg_s"])) for row in closed_b) == 0
        and max(abs(float(row["bottom_salt_in_kg_s"])) for row in closed_b) == 0,
        {"last": closed_b[-1]},
    )
    record(
        "evaporation case changes state",
        float(evap_h[-1]["vapour_kg"]) > 0 and abs(float(evap_h[-1]["soil_water_kg"]) - float(evap_h[0]["soil_water_kg"])) > 0,
        {"last": evap_h[-1]},
    )
    old_h = rows(root / stage / "no_flux_regression_old/coupled_history.csv")
    shared = (
        "time_s",
        "top_soil_mass_fraction",
        "top_soil_sw",
        "soil_salt_kg",
        "soil_water_kg",
        "vapour_kg",
    )
    regression_error = max(
        abs(float(new[key]) - float(old[key]))
        for new, old in zip(closed_h, old_h)
        for key in shared
    )
    record(
        "frozen no-flux executable regression",
        len(closed_h) == len(old_h) and regression_error < 1e-13,
        {"new_rows": len(closed_h), "old_rows": len(old_h), "max_abs": regression_error},
    )
elif stage == "pilot":
    path = root / stage / "pilot_ref_ale64"
    history, budget = audit_case(path, "reservoir")
    end = float(history[-1]["time_s"])
    active = [row for row in history if int(row["active"]) == 1]
    record("pilot reaches 50 h", abs(end - 180000.0) < 1e-8, {"end_s": end})
    record("pilot accepts ALE activation", bool(active), {"active_rows": len(active)})
    post = end - float(active[0]["activation_s"]) if active else 0.0
    # One complete registered 100 s field-output interval plus nonzero ALE
    # height/reaction evolution is the minimum nontrivial accepted window.
    record("pilot has nontrivial post-activation window", post >= 100.0, {"post_activation_s": post, "minimum_s": 100.0})
    bottom_water = float(budget[-1]["cumulative_bottom_water_kg"])
    record("reservoir supplies net water", bottom_water > 0, {"cumulative_bottom_water_kg": bottom_water})
    record(
        "nontrivial ALE fields",
        bool(active)
        and float(active[-1]["h_m"]) > 1e-7
        and float(active[-1]["internal_positive_kg_m2"]) > 0
        and float(active[-1]["epsilon_min"]) < 0.43,
        active[-1] if active else {},
    )
    attempts = rows(path / "donor_attempts.csv")
    accepted = [row for row in attempts if int(row["sufficient"]) == 1]
    transfer_rel = 0.0
    for row in accepted:
        transfer_rel = max(
            transfer_rel,
            abs(float(row["salt_transfer_residual_kg"])) / max(abs(float(row["initial_salt_kg"])), 1e-30),
            abs(float(row["water_transfer_residual_kg"])) / max(abs(float(row["initial_water_kg"])), 1e-30),
        )
    record("activation transfer conservation", bool(accepted) and transfer_rel <= 1e-12, {"relative": transfer_rel, "accepted_attempts": len(accepted)})
else:
    raise ValueError(stage)

result = {"stage": stage, "passed": all(item["passed"] for item in checks), "checks": checks}
print(json.dumps(result, indent=2))
raise SystemExit(0 if result["passed"] else 1)
