"""Compare the new M1 history against the pinned pre-task M1 history."""
import csv
import hashlib
import json
import os
from pathlib import Path

if not os.environ.get("SLURM_JOB_ID") or os.uname().nodename.startswith(("ibcln", "login", "test")):
    raise RuntimeError("Slurm compute allocation required")

root = Path.cwd()
old_root = Path("/WORK/hulm_work/zhixin/prat_analytical_boundary_20261109")
old_binary = old_root / "build/soil_solver"
expected_old_binary = "69b2a5b7f5d580f66bc43a0ba7b619ae35d4906aaa43885faa72cdc65debbccf"
assert hashlib.sha256(old_binary.read_bytes()).hexdigest() == expected_old_binary
old_history = old_root / "runs/x009_analytical_kr0p001_n80_dt0p5/history.csv"
new_history = root / "runs/m1_dt0p5/history.csv"
old = list(csv.DictReader(old_history.open()))
new = list(csv.DictReader(new_history.open()))
assert len(old) == len(new), (len(old), len(new))
mapping = {
    "time_s": "time_s", "dt_s": "dt_s", "active": "active", "activation_s": "activation_s",
    "h_m": "h_m", "solid_kg_m2": "solid_kg_m2", "pore_salt_kg_m2": "pore_salt_kg_m2",
    "pore_water_kg_m2": "pore_water_kg_m2", "j_kg_m2_s": "j_kg_m2_s", "v_m_s": "v_m_s",
    "pe_u": "pe_u", "pe_move": "pe_move", "c_top": "c_top",
    "c_soil": "c_soil", "sw_soil": "sw_soil", "total_out_kg_m2_s": "total_out_kg_m2_s",
    "salt_out_kg_m2_s": "salt_out_kg_m2_s", "soil_salt_kg": "soil_salt_kg",
    "soil_water_kg": "soil_water_kg", "vapour_kg": "vapour_kg", "bottom_salt_kg": "bottom_salt_kg",
    "bottom_water_kg": "bottom_water_kg", "salt_residual_rel": "salt_residual_rel",
    "water_residual_rel": "water_residual_rel", "precip_flux": "precip_flux",
    "fill_liquid_flux": "fill_liquid_flux",
}
max_abs = {}
max_normalized = {}
for new_key, old_key in mapping.items():
    if new_key == "active":
        assert all(a[new_key] == b[old_key] for a, b in zip(new, old))
        max_abs[new_key] = 0.0
        max_normalized[new_key] = 0.0
    else:
        pairs = list(zip(new, old))
        differences = [abs(float(a[new_key]) - float(b[old_key])) for a, b in pairs]
        scale = max(max(abs(float(b[old_key])) for _, b in pairs), 1e-12)
        tolerance = 64*2.220446049250313e-16*scale
        max_abs[new_key] = max(differences)
        max_normalized[new_key] = max_abs[new_key]/tolerance
assert max(max_normalized.values()) <= 1.0, {"max_abs": max_abs, "max_normalized": max_normalized}
result = {
    "status": "PASS",
    "rows": len(new),
    "max_absolute_differences": max_abs,
    "max_roundoff_normalized_to_64eps": max_normalized,
    "diagnostic_exclusion": "legacy Da_top was logged at post-step height; the new audit logs the value actually used at pre-step height. Dynamics (h, v, fluxes and inventories) are compared pointwise.",
    "old_binary_sha256": expected_old_binary,
    "old_history_sha256": hashlib.sha256(old_history.read_bytes()).hexdigest(),
    "new_history_sha256": hashlib.sha256(new_history.read_bytes()).hexdigest(),
    "slurm_job_id": os.environ["SLURM_JOB_ID"],
    "hostname": os.uname().nodename,
}
(root / "unit/m1_regression.json").write_text(json.dumps(result, indent=2))
print(json.dumps(result, indent=2))
