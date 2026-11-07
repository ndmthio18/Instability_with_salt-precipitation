import csv
import json
import pathlib

import analytical_growth as ag


HERE = pathlib.Path(__file__).resolve().parent
OUT = HERE / "results"
OUT.mkdir(exist_ok=True)

p = ag.Parameters()
h0 = 1.0e-5
j = 1.8e-5
t_end = 2000.0
exact = ag.analytical_height(p, h0, j, t_end)

rows = []
for dt in (40.0, 20.0, 10.0, 5.0):
    state = ag.activation_trial(p, h0)
    for _ in range(round(t_end/dt)):
        trial = ag.evaluate_trial(state, p, j, dt)
        ag.commit(state, trial)
    rows.append({"dt_s": dt, "height_m": state.height,
                 "reference_height_m": exact,
                 "absolute_error_m": abs(state.height-exact),
                 "solid_salt_kg_m2": state.solid_salt,
                 "dissolved_salt_kg_m2": state.dissolved_salt,
                 "water_kg_m2": state.water,
                 "cumulative_evaporation_kg_m2": state.cumulative_evaporation})

with (OUT / "isolated_timestep_convergence.csv").open("w", newline="") as f:
    writer = csv.DictWriter(f, fieldnames=rows[0])
    writer.writeheader()
    writer.writerows(rows)

state = ag.activation_trial(p, h0)
with (OUT / "isolated_trajectory.csv").open("w", newline="") as f:
    fields = ["time_s", "evaporation_kg_m2_s", "height_m", "velocity_m_s",
              "top_mass_fraction", "solid_salt_kg_m2",
              "dissolved_salt_kg_m2", "water_kg_m2", "solution_peclet",
              "moving_peclet", "salt_flux_kg_m2_s", "water_flux_kg_m2_s"]
    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()
    for _ in range(200):
        trial = ag.evaluate_trial(state, p, j, 10.0)
        writer.writerow({
            "time_s": state.time, "evaporation_kg_m2_s": j,
            "height_m": state.height, "velocity_m_s": trial.velocity,
            "top_mass_fraction": trial.top_mass_fraction,
            "solid_salt_kg_m2": state.solid_salt,
            "dissolved_salt_kg_m2": state.dissolved_salt,
            "water_kg_m2": state.water,
            "solution_peclet": trial.solution_peclet,
            "moving_peclet": trial.moving_peclet,
            "salt_flux_kg_m2_s": trial.salt_flux,
            "water_flux_kg_m2_s": trial.water_flux,
        })
        ag.commit(state, trial)

summary = {
    "equilibrium_mass_fraction": p.equilibrium_mass_fraction,
    "solution_density_kg_m3": p.rho,
    "effective_diffusivity_m2_s": p.effective_diffusivity,
    "height_coefficient_per_m": p.height_coefficient,
    "reference_height_m_at_2000s": exact,
    "timestep_rows": rows,
}
(OUT / "isolated_summary.json").write_text(
    json.dumps(summary, indent=2) + "\n", encoding="utf-8")
