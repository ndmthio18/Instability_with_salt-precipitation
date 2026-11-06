#!/usr/bin/env python3
"""Analyze a matched P1 B0/B1 pair without redefining onset terminology."""

from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def read_budget(path: Path) -> list[dict[str, float]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return [{key: float(value) for key, value in row.items()}
                for row in csv.DictReader(handle)]


def slope(xs: list[float], ys: list[float]) -> float | None:
    if len(xs) < 3:
        return None
    xm = sum(xs) / len(xs)
    ym = sum(ys) / len(ys)
    denominator = sum((x - xm) ** 2 for x in xs)
    return sum((x - xm) * (y - ym) for x, y in zip(xs, ys)) / denominator if denominator else None


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("b0", type=Path)
    parser.add_argument("b1", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)

    b0 = read_budget(args.b0 / "component_budget.csv")
    b1 = read_budget(args.b1 / "component_budget.csv")
    times_b0 = [row["time_s"] for row in b0]
    times = [row["time_s"] for row in b1]
    amplitudes = [row["mode_amplitude"] for row in b1]
    positive = [(time, amplitude) for time, amplitude in zip(times, amplitudes) if amplitude > 0.0]
    minimum_index = min(range(len(amplitudes)), key=amplitudes.__getitem__)
    detected_minimum_time = times[minimum_index]
    fit = [(time, math.log(amplitude)) for time, amplitude in positive
           if time >= detected_minimum_time and time <= times[-1]]
    fitted_growth_rate = slope([item[0] for item in fit], [item[1] for item in fit])

    components = ("water", "air", "salt")
    result = {
        "time_end_s": times[-1],
        "accepted_steps": len(times) - 1,
        "initial_mode_amplitude": amplitudes[0],
        "minimum_mode_amplitude": amplitudes[minimum_index],
        "new_model_mode_minimum_time_s": detected_minimum_time,
        "post_minimum_log_slope_1_s": fitted_growth_rate,
        "reference_times_s": {
            "lsa_global_neutral_crossing": 31409.1577927,
            "lsa_observable_mode_crossing": 45474.3545908,
            "richards_detected_onset": 50394.3336331
        },
        "max_abs_budget_residual_mol_m2": {
            label: max(abs(row[f"{component}_residual_mol_m2"])
                       for row in data for component in components)
            for label, data in (("B0", b0), ("B1", b1))
        },
        "final": {
            label: {key: data[-1][key] for key in
                    ("top_salt_mean", "top_liquid_density_kg_m3", "top_sw",
                     "mode_amplitude", "water_inventory_mol_m2", "salt_inventory_mol_m2")}
            for label, data in (("B0", b0), ("B1", b1))
        },
        "interpretation_limit": "The mode-amplitude minimum is a numerical detection diagnostic, not the LSA global neutral crossing or a nonlinear transport-onset time."
    }
    (args.output / "benchmark_summary.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")

    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print(json.dumps(result, indent=2))
        return 0

    hours_b0 = [time / 3600.0 for time in times_b0]
    hours = [time / 3600.0 for time in times]
    fig, axes = plt.subplots(3, 1, figsize=(7.2, 8.0), sharex=True, constrained_layout=True)
    axes[0].plot(hours_b0, [row["top_salt_mean"] for row in b0], label="B0")
    axes[0].plot(hours, [row["top_salt_mean"] for row in b1], "--", label="B1")
    axes[0].set_ylabel("top NaCl mole fraction")
    axes[0].legend()
    axes[1].semilogy(hours, amplitudes, color="tab:red", label="B1 mode 28")
    axes[1].axvline(31409.1577927 / 3600.0, color="0.5", ls=":", label="LSA global crossing")
    axes[1].axvline(45474.3545908 / 3600.0, color="tab:blue", ls=":", label="LSA mode-28 crossing")
    axes[1].axvline(50394.3336331 / 3600.0, color="tab:green", ls=":", label="Richards detection")
    axes[1].set_ylabel("surface mode amplitude")
    axes[1].legend(fontsize=8)
    axes[2].plot(hours_b0, [row["salt_inventory_mol_m2"] for row in b0], label="B0 salt")
    axes[2].plot(hours, [row["salt_inventory_mol_m2"] for row in b1], "--", label="B1 salt")
    axes[2].set_ylabel("NaCl inventory (mol m$^{-2}$)")
    axes[2].set_xlabel("time (h)")
    axes[2].legend()
    fig.savefig(args.output / "B0_B1_COMPARISON.png", dpi=200)
    plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
