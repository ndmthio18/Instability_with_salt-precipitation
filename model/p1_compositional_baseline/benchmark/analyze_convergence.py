#!/usr/bin/env python3
"""Evaluate the frozen one-factor P1 convergence campaign.

An endpoint mode-amplitude minimum is reported as an unresolved onset/growth
diagnostic, not converted into a detection time or fitted growth rate.
"""

from __future__ import annotations

import argparse
import csv
from datetime import datetime, timezone
import hashlib
import json
import math
import platform
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np


COMPONENTS = ("water", "air", "salt")
PROFILE_FIELDS = ("p_liq", "S_liq", "x^NaCl_liq", "rho_liq")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def read_budget(path: Path) -> list[dict[str, float]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return [
            {key: float(value) for key, value in row.items()}
            for row in csv.DictReader(handle)
        ]


def linear_slope(x: np.ndarray, y: np.ndarray) -> float | None:
    if x.size < 3:
        return None
    centred = x - np.mean(x)
    denominator = float(np.dot(centred, centred))
    return float(np.dot(centred, y - np.mean(y)) / denominator) if denominator else None


def summarize_budget(case_dir: Path) -> tuple[dict[str, object], list[dict[str, float]]]:
    rows = read_budget(case_dir / "component_budget.csv")
    amplitude = np.asarray([row["mode_amplitude"] for row in rows])
    times = np.asarray([row["time_s"] for row in rows])
    minimum_index = int(np.argmin(amplitude))
    endpoint = minimum_index == len(rows) - 1
    positive = amplitude > 0.0
    fit_mask = positive & (times >= times[minimum_index]) if not endpoint else np.zeros_like(positive)
    growth_rate = linear_slope(times[fit_mask], np.log(amplitude[fit_mask]))
    maximum_residual = max(
        abs(row[f"{component}_residual_mol_m2"])
        for row in rows
        for component in COMPONENTS
    )
    summary = {
        "time_end_s": rows[-1]["time_s"],
        "accepted_steps": len(rows) - 1,
        "initial_mode_amplitude": amplitude[0],
        "final_mode_amplitude": amplitude[-1],
        "minimum_mode_amplitude": amplitude[minimum_index],
        "mode_minimum_time_s": times[minimum_index],
        "mode_minimum_is_endpoint": endpoint,
        "mode_monotonic_nonincreasing": bool(np.all(np.diff(amplitude) <= 0.0)),
        "numerical_detection_time_s": None if endpoint else times[minimum_index],
        "post_minimum_growth_rate_1_s": growth_rate,
        "final_top_salt_mean": rows[-1]["top_salt_mean"],
        "final_top_sw": rows[-1]["top_sw"],
        "final_salt_inventory_mol_m2": rows[-1]["salt_inventory_mol_m2"],
        "max_abs_component_budget_residual_mol_m2": maximum_residual,
    }
    return summary, rows


def read_vtu_profile(path: Path) -> dict[str, np.ndarray]:
    arrays: dict[str, np.ndarray] = {}
    for item in ET.parse(path).iter("DataArray"):
        name = item.attrib.get("Name")
        if name:
            dtype = int if item.attrib.get("type", "").startswith(("Int", "UInt")) else float
            arrays[name] = np.fromstring(item.text or "", sep=" ", dtype=dtype)
    coordinates = arrays["Coordinates"].reshape(-1, 3)
    connectivity = arrays["connectivity"]
    offsets = arrays["offsets"]
    centres = []
    start = 0
    for stop in offsets:
        centres.append(coordinates[connectivity[start:stop]].mean(axis=0))
        start = int(stop)
    centres = np.asarray(centres)
    z_values = np.unique(centres[:, 1])
    profile: dict[str, list[float]] = {"z_m": []}
    profile.update({field: [] for field in PROFILE_FIELDS})
    for z_value in z_values:
        selected = np.isclose(centres[:, 1], z_value, rtol=0.0, atol=1.0e-12)
        profile["z_m"].append(float(z_value))
        for field in PROFILE_FIELDS:
            profile[field].append(float(np.mean(arrays[field][selected])))
    return {key: np.asarray(value) for key, value in profile.items()}


def final_vtu(case_dir: Path) -> Path:
    paths = sorted(case_dir.glob("*.vtu"))
    if not paths:
        raise FileNotFoundError(f"no VTU output in {case_dir}")
    return paths[-1]


def relative_profile_norms(
    reference: dict[str, np.ndarray], candidate: dict[str, np.ndarray]
) -> dict[str, dict[str, float]]:
    result = {}
    z_ref = reference["z_m"]
    for field in PROFILE_FIELDS:
        ref = reference[field]
        candidate_on_ref = np.interp(z_ref, candidate["z_m"], candidate[field])
        difference = candidate_on_ref - ref
        result[field] = {
            "relative_L1": float(np.mean(np.abs(difference)) / max(np.mean(np.abs(ref)), 1.0e-30)),
            "relative_L2": float(np.linalg.norm(difference) / max(np.linalg.norm(ref), 1.0e-30)),
            "relative_Linf": float(np.max(np.abs(difference)) / max(np.max(np.abs(ref)), 1.0e-30)),
        }
    return result


def relative_difference(candidate: float, reference: float) -> float:
    return abs(candidate - reference) / max(abs(candidate), abs(reference), 1.0e-30)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base", type=Path, required=True)
    parser.add_argument("--horizontal-coarse", type=Path, required=True)
    parser.add_argument("--vertical-coarse", type=Path, required=True)
    parser.add_argument("--time-half", type=Path, required=True)
    parser.add_argument("--plan", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--binary-sha256", required=True)
    parser.add_argument("--remote-result", required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    plan = json.loads(args.plan.read_text(encoding="utf-8"))
    case_paths = {
        "base": args.base,
        "horizontal_coarse": args.horizontal_coarse,
        "vertical_coarse": args.vertical_coarse,
        "time_half": args.time_half,
    }
    summaries = {}
    histories = {}
    profiles = {}
    input_hashes = {}
    for label, path in case_paths.items():
        summaries[label], histories[label] = summarize_budget(path)
        vtu = final_vtu(path)
        profiles[label] = read_vtu_profile(vtu)
        input_hashes[label] = {
            "component_budget_sha256": sha256(path / "component_budget.csv"),
            "final_vtu": vtu.name,
            "final_vtu_sha256": sha256(vtu),
        }
    base = summaries["base"]
    acceptance = plan["acceptance"]
    comparisons = {}
    for label in ("horizontal_coarse", "vertical_coarse", "time_half"):
        candidate = summaries[label]
        slope_available = (
            base["post_minimum_growth_rate_1_s"] is not None
            and candidate["post_minimum_growth_rate_1_s"] is not None
        )
        slope_difference = (
            relative_difference(
                float(candidate["post_minimum_growth_rate_1_s"]),
                float(base["post_minimum_growth_rate_1_s"]),
            )
            if slope_available else None
        )
        checks = {
            "budget": float(candidate["max_abs_component_budget_residual_mol_m2"])
            <= float(acceptance["inventory_absolute_budget_residual_mol_m2"]),
            "surface_salinity": relative_difference(
                float(candidate["final_top_salt_mean"]), float(base["final_top_salt_mean"])
            ) <= float(acceptance["surface_salinity_relative_difference"]),
            "mode_amplitude": relative_difference(
                float(candidate["final_mode_amplitude"]), float(base["final_mode_amplitude"])
            ) <= float(acceptance["mode_amplitude_relative_difference_when_resolved"]),
            "post_minimum_growth_rate": (
                slope_available
                and slope_difference <= float(acceptance["post_minimum_growth_rate_relative_difference"])
            ),
        }
        comparisons[label] = {
            "final_surface_salinity_relative_difference": relative_difference(
                float(candidate["final_top_salt_mean"]), float(base["final_top_salt_mean"])
            ),
            "final_mode_amplitude_relative_difference": relative_difference(
                float(candidate["final_mode_amplitude"]), float(base["final_mode_amplitude"])
            ),
            "post_minimum_growth_rate_relative_difference": slope_difference,
            "growth_rate_comparison_available": slope_available,
            "profile_norms_after_vertical_interpolation": relative_profile_norms(
                profiles["base"], profiles[label]
            ),
            "checks": checks,
            "pair_pass": all(checks.values()),
        }
    overall = all(value["pair_pass"] for value in comparisons.values())
    result = {
        "overall_pass": overall,
        "failure_semantics": (
            None if overall else
            "One or more registered checks failed. In particular, an endpoint mode minimum leaves no post-minimum growth rate, so the onset-sensitive gate cannot pass even when terminal observables are insensitive."
        ),
        "plan_sha256": sha256(args.plan),
        "summaries": summaries,
        "comparisons_to_base": comparisons,
        "input_hashes": input_hashes,
    }
    (args.output / "convergence_assessment.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )

    with (args.output / "profile_norms.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=("comparison", "field", "relative_L1", "relative_L2", "relative_Linf"),
        )
        writer.writeheader()
        for comparison, values in comparisons.items():
            for field, norms in values["profile_norms_after_vertical_interpolation"].items():
                writer.writerow({"comparison": comparison, "field": field, **norms})

    import matplotlib
    import matplotlib.pyplot as plt
    styles = {
        "base": ("#000000", "-"),
        "horizontal_coarse": ("#0072B2", "--"),
        "vertical_coarse": ("#D55E00", "-."),
        "time_half": ("#009E73", ":"),
    }
    labels = {
        "base": "base (280 × 100, Δt)",
        "horizontal_coarse": "horizontal coarse (140 × 100, Δt)",
        "vertical_coarse": "vertical coarse (280 × 50, Δt)",
        "time_half": "time refined (280 × 100, Δt/2)",
    }
    fig, axes = plt.subplots(3, 1, figsize=(7.2, 8.0), sharex=True, constrained_layout=True)
    for label, history in histories.items():
        color, line = styles[label]
        time_h = np.asarray([row["time_s"] for row in history]) / 3600.0
        axes[0].plot(time_h, [row["top_salt_mean"] for row in history], color=color, ls=line, label=labels[label])
        axes[1].semilogy(time_h, [row["mode_amplitude"] for row in history], color=color, ls=line, label=labels[label])
        residual = [max(abs(row[f"{component}_residual_mol_m2"]) for component in COMPONENTS) for row in history]
        axes[2].semilogy(time_h, np.maximum(residual, 1.0e-18), color=color, ls=line, label=labels[label])
    axes[0].set_ylabel("top NaCl mole fraction")
    axes[1].set_ylabel("mode-28 amplitude")
    axes[2].set_ylabel("max budget residual\n(mol m$^{-2}$)")
    axes[2].set_xlabel("time (h)")
    axes[0].legend(fontsize=8)
    for axis in axes:
        axis.grid(alpha=0.25)
    for suffix in ("png", "svg"):
        fig.savefig(args.output / f"CONVERGENCE_TRAJECTORIES.{suffix}", dpi=300, facecolor="white")
    plt.close(fig)

    fig, axes = plt.subplots(1, 2, figsize=(7.6, 4.0), sharey=True, constrained_layout=True)
    for label, profile in profiles.items():
        color, line = styles[label]
        axes[0].plot(profile["x^NaCl_liq"], profile["z_m"], color=color, ls=line, label=labels[label])
        axes[1].plot(profile["S_liq"], profile["z_m"], color=color, ls=line, label=labels[label])
    axes[0].set(xlabel="NaCl mole fraction", ylabel="z (m)")
    axes[1].set_xlabel("liquid saturation")
    axes[0].legend(fontsize=7)
    for axis in axes:
        axis.grid(alpha=0.25)
    for suffix in ("png", "svg"):
        fig.savefig(args.output / f"CONVERGENCE_PROFILES.{suffix}", dpi=300, facecolor="white")
    plt.close(fig)

    generated_files = sorted(
        path for path in args.output.iterdir()
        if path.is_file() and path.name not in {"manifest.json", "SHA256SUMS.txt"}
    )
    manifest = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "generator": str(Path(__file__).resolve()),
        "generator_sha256": sha256(Path(__file__)),
        "python": platform.python_version(),
        "platform": platform.platform(),
        "numpy": np.__version__,
        "matplotlib": matplotlib.__version__,
        "accepted_binary_sha256": args.binary_sha256,
        "remote_result_root": args.remote_result,
        "plan": str(args.plan.resolve()),
        "plan_sha256": sha256(args.plan),
        "inputs": input_hashes,
        "outputs": {path.name: sha256(path) for path in generated_files},
    }
    (args.output / "manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    checksum_files = sorted(
        path for path in args.output.iterdir()
        if path.is_file() and path.name != "SHA256SUMS.txt"
    )
    (args.output / "SHA256SUMS.txt").write_text(
        "".join(f"{sha256(path)}  {path.name}\n" for path in checksum_files),
        encoding="utf-8",
    )
    print(json.dumps(result, indent=2))
    return 0 if overall else 2


if __name__ == "__main__":
    raise SystemExit(main())
