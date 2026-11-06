#!/usr/bin/env python3
"""Extract horizontally averaged P1 fields from ASCII VTU snapshots."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path
import xml.etree.ElementTree as ET

import numpy as np


FIELDS = ("p_liq", "S_liq", "x^NaCl_liq", "rho_liq")


def read_vtu(path: Path):
    arrays = {}
    for item in ET.parse(path).iter("DataArray"):
        name = item.attrib.get("Name")
        if name:
            dtype = int if item.attrib.get("type", "").startswith(("Int", "UInt")) else float
            arrays[name] = np.fromstring(item.text or "", sep=" ", dtype=dtype)
    coordinates = arrays["Coordinates"].reshape(-1, 3)
    connectivity = arrays["connectivity"]
    offsets = arrays["offsets"]
    centers = []
    start = 0
    for stop in offsets:
        centers.append(coordinates[connectivity[start:stop]].mean(axis=0))
        start = stop
    return np.asarray(centers), arrays


def profile(path: Path):
    centers, arrays = read_vtu(path)
    z_values = np.unique(centers[:, 1])
    result = []
    velocity = arrays["velocity_liq (m/s)"].reshape(-1, 3)
    for z in z_values:
        mask = np.isclose(centers[:, 1], z, rtol=0.0, atol=1e-12)
        row = {"z_m": float(z)}
        for field in FIELDS:
            row[field] = float(arrays[field][mask].mean())
        row["velocity_liq_z_mean_m_s"] = float(velocity[mask, 1].mean())
        row["velocity_liq_norm_max_m_s"] = float(np.linalg.norm(velocity[mask], axis=1).max())
        result.append(row)
    return result


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("b0_vtu", type=Path)
    parser.add_argument("b1_vtu", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=False)
    profiles = {"B0": profile(args.b0_vtu), "B1": profile(args.b1_vtu)}

    keys = tuple(profiles["B0"][0])
    with (args.output / "vertical_profiles.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=("case",) + keys)
        writer.writeheader()
        for label, data in profiles.items():
            for row in data:
                writer.writerow({"case": label, **row})

    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(2, 2, figsize=(8.0, 7.0), sharey=True, constrained_layout=True)
    labels = (("p_liq", "liquid pressure (Pa)"),
              ("S_liq", "liquid saturation"),
              ("x^NaCl_liq", "NaCl mole fraction"),
              ("rho_liq", "liquid density (kg m$^{-3}$)"))
    for axis, (field, xlabel) in zip(axes.ravel(), labels):
        for case, data in profiles.items():
            axis.plot([row[field] for row in data], [row["z_m"] for row in data], label=case)
        axis.set_xlabel(xlabel)
        axis.set_ylabel("z (m)")
        axis.legend()
    fig.savefig(args.output / "B0_B1_VERTICAL_PROFILES.png", dpi=200)
    plt.close(fig)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
