import csv
import json
import math
import pathlib
import sys


def read_rows(path):
    with path.open(newline="") as stream:
        return list(csv.DictReader(stream))


def f(row, key):
    return float(row[key])


def finite_or_none(value):
    return value if math.isfinite(value) else None


def summarize(run):
    budget = read_rows(run / "component_budget.csv")
    steps = read_rows(run / "p2f_accepted_steps.csv")
    if not budget or not steps:
        raise RuntimeError(f"{run}: missing accepted records")
    initial, final = budget[0], budget[-1]
    final_budget_time = f(final, "time_s")
    aligned = [row for row in steps
               if f(row, "time_s") <= final_budget_time + 1.0e-10
               and int(row["valid"])]
    if not aligned:
        raise RuntimeError(f"{run}: no crust record aligned to accepted soil state")
    last = aligned[-1]
    exit_status = int((run / "exit_status.txt").read_text().strip())
    width = 1.0e-3
    crust_salt = width * (f(last, "solid_salt_kg_m2")
                          + f(last, "dissolved_salt_kg_m2"))
    crust_water = width * f(last, "crust_water_kg_m2")
    evaporation = width * f(last, "evaporation_kg_m2_s") * f(last, "time_s")
    initial_salt = f(initial, "salt_mass_kg")
    initial_water = f(initial, "water_mass_kg")
    salt_residual = f(final, "salt_mass_kg") + crust_salt - initial_salt
    water_residual = f(final, "water_mass_kg") + crust_water + evaporation - initial_water
    salt_scale = max(abs(initial_salt), 1.0e-30)
    water_scale = max(abs(initial_water), 1.0e-30)
    active = [row for row in steps if int(row["active"])]
    activation = [row for row in steps if int(row["activation_event"])]
    invalid = [row for row in steps if not int(row["valid"])]
    ledger_residual = None
    ledger_relative = None
    ledger_path = run / "p2f_activation_ledger.csv"
    if ledger_path.exists():
        ledger = read_rows(ledger_path)[0]
        ledger_residual = max(abs(f(ledger, "salt_residual_kg")),
                              abs(f(ledger, "water_residual_kg")))
        ledger_scale = max(abs(f(ledger, "soil_salt_before_kg")),
                           abs(f(ledger, "soil_water_before_kg")), 1.0e-30)
        ledger_relative = ledger_residual / ledger_scale
    return {
        "run": run.name,
        "exit_status": exit_status,
        "accepted_time_s": final_budget_time,
        "accepted_steps": len(budget)-1,
        "domain_exit_time_s": (f(invalid[-1], "time_s") if invalid else None),
        "activation_time_s": (f(activation[0], "time_s") if activation else None),
        "height_m": f(last, "height_m"),
        "velocity_m_s": f(last, "velocity_m_s"),
        "interface_mass_fraction": f(last, "interface_mass_fraction"),
        "interface_saturation": f(last, "interface_saturation"),
        "solution_pe_max": max(abs(f(row, "solution_pe")) for row in steps),
        "moving_pe_max": max(abs(f(row, "moving_pe")) for row in steps),
        "salt_residual_kg": salt_residual,
        "water_residual_kg": water_residual,
        "salt_residual_relative": abs(salt_residual) / salt_scale,
        "water_residual_relative": abs(water_residual) / water_scale,
        "activation_ledger_absolute_residual_kg": ledger_residual,
        "activation_ledger_relative_residual": ledger_relative,
        "active_records": len(active),
        "invalid_records": len(invalid),
    }


root = pathlib.Path(sys.argv[1]).resolve()
rows = []
for path in sorted(root.iterdir()):
    if (path.is_dir() and path.name.startswith("grid")
            and (path / "component_budget.csv").exists()):
        rows.append(summarize(path))

by_name = {row["run"]: row for row in rows}
baseline = by_name.get("grid80_seed10_dt0025")
if baseline:
    for row in rows:
        if row["run"].startswith("grid") and row["activation_time_s"] is not None:
            row["height_relative_to_baseline"] = row["height_m"] / baseline["height_m"] - 1.0
            row["activation_relative_to_baseline"] = (
                row["activation_time_s"] / baseline["activation_time_s"] - 1.0)

grid = [by_name.get(name) for name in (
    "grid40_seed10_dt0025", "grid80_seed10_dt0025", "grid160_seed10_dt0025")]
grid_ready = all(row is not None and row["exit_status"] == 0
                 and row["activation_time_s"] is not None for row in grid)
if grid_ready:
    grid_activation_spread = max(row["activation_time_s"] for row in grid) / min(row["activation_time_s"] for row in grid) - 1.0
    grid_interface_spread = max(row["interface_mass_fraction"] for row in grid) / min(row["interface_mass_fraction"] for row in grid) - 1.0
    grid_height_spread = max(row["height_m"] for row in grid) / min(row["height_m"] for row in grid) - 1.0
else:
    grid_activation_spread = grid_interface_spread = grid_height_spread = math.inf

half_dt = by_name.get("grid80_seed10_dt00125")
dt_ready = (baseline is not None and half_dt is not None
            and baseline["exit_status"] == 0 and half_dt["exit_status"] == 0
            and baseline["activation_time_s"] is not None
            and half_dt["activation_time_s"] is not None)
dt_height = (abs(half_dt["height_m"]/baseline["height_m"]-1.0)
             if dt_ready else math.inf)
dt_activation = (abs(half_dt["activation_time_s"]
                     / baseline["activation_time_s"]-1.0)
                 if dt_ready else math.inf)

out = root / "coupled_summary.csv"
fields = sorted({key for row in rows for key in row})
with out.open("w", newline="") as stream:
    writer = csv.DictWriter(stream, fieldnames=fields)
    writer.writeheader()
    writer.writerows(rows)

criteria = json.loads((root.parents[1] / "source" / "model" / "p2_f_analytical_growth"
                       / "verification" / "acceptance_criteria.json").read_text())
tol = criteria["tolerances"]
checks = {
    "all_nominal_runs_activated": bool(rows) and all(row["active_records"] > 0 for row in rows),
    "all_nominal_runs_completed": bool(rows) and all(row["exit_status"] == 0 for row in rows),
    "combined_salt_budget": bool(rows) and all(row["salt_residual_relative"] <= tol["coupled_budget_relative"] for row in rows),
    "combined_water_budget": bool(rows) and all(row["water_residual_relative"] <= tol["coupled_budget_relative"] for row in rows),
    "activation_ledger": bool(rows) and all(row["activation_ledger_relative_residual"] is not None and row["activation_ledger_relative_residual"] <= tol["activation_budget_relative"] for row in rows),
    "peclet_validity": bool(rows) and all(row["solution_pe_max"] <= tol["maximum_solution_peclet"] and row["moving_pe_max"] <= tol["maximum_moving_peclet"] for row in rows),
    "at_least_one_adequate_growth_run": any(row["exit_status"] == 0 and row["active_records"] > 1 and row["invalid_records"] == 0 for row in rows),
    "grid_activation_time": grid_activation_spread <= tol["grid_activation_time_successive_relative"],
    "grid_interface_composition": grid_interface_spread <= tol["grid_top_concentration_successive_relative"],
    "grid_height": grid_height_spread <= tol["grid_evaporation_successive_relative"],
    "time_step_height": dt_height <= tol["timestep_successive_relative"],
    "time_step_activation": dt_activation <= tol["grid_activation_time_successive_relative"],
}
inadequate = root / "inadequate_supply"
inadequate_steps = read_rows(inadequate / "p2f_accepted_steps.csv")
inadequate_status = int((inadequate / "exit_status.txt").read_text().strip())
checks["insufficient_supply_flagged"] = (
    inadequate_status != 0 and any(not int(row["valid"]) for row in inadequate_steps))
(root / "coupled_acceptance.json").write_text(json.dumps({
    "checks": checks,
    "metrics": {
        "grid_activation_spread_relative": finite_or_none(grid_activation_spread),
        "grid_interface_spread_relative": finite_or_none(grid_interface_spread),
        "grid_height_spread_relative": finite_or_none(grid_height_spread),
        "dt_height_relative": finite_or_none(dt_height),
        "dt_activation_relative": finite_or_none(dt_activation),
    },
    "pass": all(checks.values()),
    "runs": rows,
}, indent=2) + "\n")
print(json.dumps(checks, indent=2))
if not all(checks.values()):
    raise SystemExit(2)
