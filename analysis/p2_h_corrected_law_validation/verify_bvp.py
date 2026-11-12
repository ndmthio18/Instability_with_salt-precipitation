"""Stage A: compare closed-form predictions with an independent solve_bvp."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

from analytical_model import (
    Parameters,
    analytical_profile,
    corrected_observables,
    solve_independent_bvp,
)


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"
RESULTS.mkdir(exist_ok=True)

DAS = (1.0, 4.0, 10.0, 20.0, 45.0, 100.0)
US = (0.0, 0.25, 0.5, 0.75, 0.9)


def relative_error(a: float, b: float, floor: float = 1e-30) -> float:
    return abs(a - b) / max(abs(a), abs(b), floor)


def main() -> None:
    p = Parameters()
    rows: list[dict[str, float]] = []
    profiles: list[pd.DataFrame] = []
    for da in DAS:
        for u in US:
            exact = corrected_observables(da, u, p)
            bvp = solve_independent_bvp(da, u, p)
            c_exact, _ = analytical_profile(bvp["z_m"], da, u, p)
            profile_abs = float(np.max(np.abs(c_exact - bvp["excess"])))
            profile_scale = max(float(np.max(np.abs(c_exact))), 1e-30)
            rows.append(
                {
                    "da": da,
                    "u": u,
                    "height_m": exact["height_m"],
                    "analytical_top_concentration": exact["top_concentration"],
                    "bvp_top_concentration": bvp["top_concentration"],
                    "top_concentration_abs_error": abs(
                        exact["top_concentration"] - bvp["top_concentration"]
                    ),
                    "analytical_top_flux_kg_m2_s": exact["top_flux_kg_m2_s"],
                    "bvp_top_flux_kg_m2_s": bvp["top_flux_kg_m2_s"],
                    "top_flux_relative_error": relative_error(
                        exact["top_flux_kg_m2_s"], bvp["top_flux_kg_m2_s"]
                    ),
                    "analytical_velocity_m_s": exact["velocity_m_s"],
                    "bvp_velocity_m_s": bvp["velocity_m_s"],
                    "velocity_relative_error": relative_error(
                        exact["velocity_m_s"], bvp["velocity_m_s"]
                    ),
                    "analytical_internal_signed_kg_m2_s": exact[
                        "internal_signed_kg_m2_s"
                    ],
                    "bvp_internal_signed_kg_m2_s": bvp[
                        "internal_signed_kg_m2_s"
                    ],
                    "internal_abs_error_kg_m2_s": abs(
                        exact["internal_signed_kg_m2_s"]
                        - bvp["internal_signed_kg_m2_s"]
                    ),
                    "profile_max_abs_error": profile_abs,
                    "profile_max_relative_error": profile_abs / profile_scale,
                    "bvp_nodes": bvp["nodes"],
                    "bvp_rms_residual_max": bvp["rms_residual_max"],
                }
            )
            if da in (1.0, 20.0, 100.0) and u in (0.0, 0.9):
                profiles.append(
                    pd.DataFrame(
                        {
                            "da": da,
                            "u": u,
                            "z_over_h": bvp["z_m"] / exact["height_m"],
                            "z_m": bvp["z_m"],
                            "analytical_concentration": p.ceq + c_exact,
                            "bvp_concentration": p.ceq + bvp["excess"],
                        }
                    )
                )
    frame = pd.DataFrame(rows)
    frame.to_csv(RESULTS / "bvp_verification.csv", index=False)
    pd.concat(profiles, ignore_index=True).to_csv(
        RESULTS / "bvp_profiles.csv", index=False
    )
    summary = {
        "status": "PASS"
        if frame["velocity_relative_error"].max() < 1e-8
        and frame["profile_max_relative_error"].max() < 1e-8
        else "FAIL",
        "cases": len(frame),
        "max_velocity_relative_error": float(frame["velocity_relative_error"].max()),
        "max_top_flux_relative_error": float(frame["top_flux_relative_error"].max()),
        "max_top_concentration_abs_error": float(
            frame["top_concentration_abs_error"].max()
        ),
        "max_profile_relative_error": float(
            frame["profile_max_relative_error"].max()
        ),
        "max_internal_abs_error_kg_m2_s": float(
            frame["internal_abs_error_kg_m2_s"].max()
        ),
        "parameters": p.to_dict(),
        "boundary_conditions": {
            "bottom": "rho*eps*Deff*dC/dz = Bsoil",
            "top": "Sref-rho*mathcalD*dC/dz = rho*ks*(C-Ceq)",
        },
    }
    (RESULTS / "bvp_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2))
    if summary["status"] != "PASS":
        raise SystemExit("independent BVP gate failed")


if __name__ == "__main__":
    main()
