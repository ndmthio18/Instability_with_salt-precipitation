"""Create the preregistered controlled ALE benchmark matrix."""

from __future__ import annotations

import csv
import json
from pathlib import Path

from analytical_model import Parameters, thickness_for_da


HERE = Path(__file__).resolve().parent
CONFIGS = HERE / "configs"
CONFIGS.mkdir(exist_ok=True)
DAS = (1.0, 4.0, 10.0, 20.0, 45.0, 100.0)
US = (0.0, 0.25, 0.5, 0.75, 0.9)


def row(case_id, suite, da, u, cells, power, rtol, max_step_div, mode):
    return {
        "case_id": case_id,
        "suite": suite,
        "da": da,
        "u": u,
        "cells": cells,
        "power": power,
        "rtol": rtol,
        "max_step_div": max_step_div,
        "end_tau": 20,
        "outputs_per_tau": 20,
        "boundary_mode": mode,
    }


def main() -> None:
    rows = []
    for da in DAS:
        for u in US:
            rows.append(
                row(f"main_da{da:g}_u{u:g}", "main", da, u, 64, 1.5, 1e-8, 50, "neumann")
            )
    representatives = ((1.0, 0.0), (1.0, 0.9), (20.0, 0.0), (20.0, 0.9),
                       (100.0, 0.0), (100.0, 0.9))
    for da, u in representatives:
        for cells in (32, 64, 128):
            for power in (1.0, 1.5, 2.0, 3.0):
                rows.append(
                    row(
                        f"space_da{da:g}_u{u:g}_n{cells}_p{power:g}",
                        "spatial", da, u, cells, power, 1e-8, 50, "neumann",
                    )
                )
        for label, rtol, step_div in (
            ("coarse", 1e-7, 25), ("baseline", 1e-8, 50),
            ("strict", 1e-9, 50), ("smallstep", 1e-8, 100),
        ):
            rows.append(
                row(
                    f"time_da{da:g}_u{u:g}_{label}",
                    "temporal", da, u, 64, 1.5, rtol, step_div, "neumann",
                )
            )
        rows.append(
            row(
                f"boundary_da{da:g}_u{u:g}_eqin", "boundary", da, u,
                64, 1.5, 1e-8, 50, "equilibrium_inflow",
            )
        )
    for da in (1.0, 20.0, 100.0):
        rows.append(
            row(
                f"boundary_da{da:g}_dirichlet", "boundary", da, 0.0,
                64, 1.5, 1e-8, 50, "dirichlet",
            )
        )
    with (CONFIGS / "cases.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    p = Parameters()
    lengths = [
        {
            "da": da,
            "height_m": thickness_for_da(da, p),
            "height_um": thickness_for_da(da, p) * 1e6,
            "height_over_ell": da**0.5,
            "ell_m": p.ell,
            "grain_m": p.grain,
            "height_over_grain": thickness_for_da(da, p) / p.grain,
        }
        for da in DAS
    ]
    manifest = {
        "case_count": len(rows),
        "main_cases": len(DAS) * len(US),
        "spatial_cases": 6 * 3 * 4,
        "temporal_cases": 6 * 4,
        "boundary_cases": 9,
        "parameters": p.to_dict(),
        "physical_lengths": lengths,
        "main_mesh": {"cells": 64, "power": 1.5},
        "main_time": {
            "duration": "20 h^2/Deff",
            "output_dt": "(h^2/Deff)/20",
            "max_step": "(h^2/Deff)/50",
        },
        "comparison_window_candidate": "10 <= t/(h^2/Deff) <= 20",
    }
    (CONFIGS / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
