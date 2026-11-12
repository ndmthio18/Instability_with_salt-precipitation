#!/usr/bin/env python3
"""Extract the accepted x0=0.08 post-precipitation one-way forcing history."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
INPUT = HERE / "inputs"
OUT = HERE / "results" / "x008_forcing.csv"


def read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def interpolate(rows: list[dict[str, float]], time_s: float) -> dict[str, float]:
    if not (rows[0]["absolute_time_s"] <= time_s <= rows[-1]["absolute_time_s"]):
        raise ValueError("interpolation time outside source support")
    for left, right in zip(rows, rows[1:]):
        if left["absolute_time_s"] <= time_s <= right["absolute_time_s"]:
            span = right["absolute_time_s"] - left["absolute_time_s"]
            weight = 0.0 if span == 0.0 else (time_s - left["absolute_time_s"]) / span
            return {
                key: left[key] + weight * (right[key] - left[key])
                for key in left
            }
    return dict(rows[-1])


def main() -> None:
    flux_path = INPUT / "flux_decomposition.csv"
    time_path = INPUT / "x008_timeseries.csv"
    summary_path = INPUT / "case_summary.csv"
    summary = next(
        row for row in read(summary_path)
        if abs(float(row["x0_mole_fraction"]) - 0.08) < 1e-12
    )
    onset = float(summary["t_precip_s"])
    evaporation = {float(row["time_s"]): float(row["j_evap_water_kg_m2_s"])
                   for row in read(time_path)}
    rows = []
    for row in read(flux_path):
        if abs(float(row["x0_mole_fraction"]) - 0.08) > 1e-12:
            continue
        absolute = float(row["time_s"])
        if absolute not in evaporation:
            raise RuntimeError(f"missing evaporation sample at {absolute}")
        diffusive = float(row["j_s_diff_kg_m2_s"])
        rows.append({
            "absolute_time_s": absolute,
            "evaporation_water_kg_m2_s": evaporation[absolute],
            "b_soil_kg_m2_s": -diffusive,
            "soil_advective_nacl_kg_m2_s": float(row["j_s_adv_kg_m2_s"]),
            "soil_net_nacl_kg_m2_s": float(row["j_s_net_kg_m2_s"]),
            "surface_phase_flux_kg_m2_s": float(row["j_s_precip_boundary_kg_m2_s"]),
        })
    rows.sort(key=lambda row: row["absolute_time_s"])
    if any(row["b_soil_kg_m2_s"] < 0.0 for row in rows):
        raise RuntimeError("x0=0.08 source contains upward diffusive samples")
    event = interpolate(rows, onset)
    selected = [event] + [row for row in rows if row["absolute_time_s"] > onset]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fields = ["elapsed_time_s", *event.keys()]
    with OUT.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields)
        writer.writeheader()
        for row in selected:
            writer.writerow({"elapsed_time_s": row["absolute_time_s"] - onset, **row})
    manifest = {
        "case": "x0=0.08",
        "precipitation_onset_absolute_s": onset,
        "forcing_duration_s": selected[-1]["absolute_time_s"] - onset,
        "source_sampling_s": 300.0,
        "interpolation": "linear only at t_precip; ALE callbacks linearly interpolate stored accepted states",
        "b_soil_definition": "-j_s_diff at the first internal soil SCVF; positive downward",
        "coupling_scope": "one-way prescribed soil history; no crust feedback to soil",
        "inputs": {path.name: sha256(path) for path in (flux_path, time_path, summary_path)},
        "output": {OUT.name: sha256(OUT)},
        "rows": len(selected),
    }
    (OUT.parent / "forcing_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
