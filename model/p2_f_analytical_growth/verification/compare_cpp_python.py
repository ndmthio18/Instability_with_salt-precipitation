import json
import pathlib
import sys


model = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(model))
import analytical_growth as ag


line = next(line for line in pathlib.Path(sys.argv[1]).read_text().splitlines()
            if line.startswith("VECTOR "))
cpp = {token.split("=", 1)[0]: float(token.split("=", 1)[1])
       for token in line.split()[1:]}
p = ag.Parameters()
python = {
    "equilibrium_mass_fraction": p.equilibrium_mass_fraction,
    "solution_density": p.rho,
    "effective_diffusivity": p.effective_diffusivity,
    "velocity": ag.growth_velocity(p, 1.0e-5, 1.8e-5),
    "reference_height": ag.analytical_height(p, 1.0e-5, 1.8e-5, 2000.0),
}
relative = {key: abs(cpp[key]-value)/max(abs(value), 1.0e-300)
            for key, value in python.items()}
result = {"cpp": cpp, "python": python, "relative_error": relative,
          "maximum_relative_error": max(relative.values())}
pathlib.Path(sys.argv[2]).write_text(json.dumps(result, indent=2) + "\n")
if result["maximum_relative_error"] > 1.0e-12:
    raise SystemExit(2)
print(json.dumps(result, indent=2))
