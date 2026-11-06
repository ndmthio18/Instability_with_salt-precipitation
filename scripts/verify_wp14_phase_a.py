"""Independent acceptance checks for WP14 Phase A delivery."""
from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path


PAPER2 = Path(__file__).resolve().parents[1]
BASE = PAPER2 / "analysis" / "wp14_nonlinear_timestep_20261106" / "phase_a"
FIGURES = BASE / "figures"
GRIDS = ("R24", "R48", "R96")
AMPLITUDES = (1e-6, 1e-5, 1e-4)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle))


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    checks = []

    def check(name: str, condition: bool, detail=None) -> None:
        checks.append({"name": name, "passed": bool(condition), "detail": detail})

    manifest = json.loads((BASE / "source_manifest.json").read_text())
    check("source manifest has 54 frozen files", len(manifest) == 54, len(manifest))
    for item in manifest:
        path = PAPER2 / item["path_relative_to_paper2"]
        check(f"source exists: {item['path_relative_to_paper2']}", path.is_file())
        if path.is_file():
            check(f"source size: {item['path_relative_to_paper2']}", path.stat().st_size == item["bytes"])
            check(f"source hash: {item['path_relative_to_paper2']}", sha256(path) == item["sha256"])

    normalized = read_csv(BASE / "normalized_harmonic_timeseries.csv")
    milestones = read_csv(BASE / "normalized_milestones.csv")
    modes = read_csv(BASE / "mode_checkpoints.csv")
    transport = read_csv(BASE / "transport_pulse_metrics.csv")
    grids = read_csv(BASE / "grid_difference_by_amplitude.csv")
    check("normalized row count", len(normalized) == 1305, len(normalized))
    check("milestone row count", len(milestones) == 27, len(milestones))
    check("mode checkpoint row count", len(modes) == 63, len(modes))
    check("transport pulse row count", len(transport) == 36, len(transport))
    check("grid difference row count", len(grids) == 60, len(grids))
    check("all normalized ratios valid", all(int(row["ratio_valid"]) == 1 for row in normalized))
    check("all predefined milestones resolved", all(row["crossing_T"] != "" for row in milestones))

    for grid in GRIDS:
        ref = next(row for row in milestones if row["grid"] == grid and float(row["amplitude"]) == 1e-6
                   and float(row["threshold"]) == 0.1)
        enhanced = next(row for row in milestones if row["grid"] == grid and float(row["amplitude"]) == 1e-4
                        and float(row["threshold"]) == 0.1)
        check(f"normalized acceleration persists {grid}", float(enhanced["crossing_T"]) < float(ref["crossing_T"]))
        normalized_fraction = 1 - float(enhanced["crossing_T"]) / float(ref["crossing_T"])
        absolute_fraction = 1 - float(enhanced["absolute_H23_milestone_T"]) / float(ref["absolute_H23_milestone_T"])
        check(f"absolute H23 exaggerates acceleration {grid}", absolute_fraction > normalized_fraction,
              {"absolute": absolute_fraction, "normalized": normalized_fraction})

    normalized_states = [row for row in modes if row["coordinate"] == "normalized_Rharm_0.1"]
    check("all normalized-milestone states remain m2 dominated",
          all(int(row["dominant_mode"]) == 2 and float(row["m2_fraction"]) > 0.98 for row in normalized_states))

    for depth in (5, 10):
        prefix = f"{depth}mm"
        row = next(item for item in grids if float(item["amplitude"]) == 1e-4
                   and item["metric"] == f"{prefix}_positive_advective_integral_m")
        check(f"A1e-4 R48/R96 positive integral within 5% at {depth} mm",
              float(row["relative_difference_R48_R96"]) < 0.05,
              float(row["relative_difference_R48_R96"]))
        row = next(item for item in grids if float(item["amplitude"]) == 1e-4
                   and item["metric"] == f"{prefix}_principal_pulse_peak_T")
        check(f"A1e-4 R48/R96 principal pulse timing within 1% at {depth} mm",
              float(row["relative_difference_R48_R96"]) < 0.01,
              float(row["relative_difference_R48_R96"]))

    for grid in GRIDS:
        max_odd = {}
        for amplitude in AMPLITUDES:
            selected = [float(row["odd_mode_fraction"]) for row in normalized
                        if row["grid"] == grid and float(row["amplitude"]) == amplitude]
            max_odd[amplitude] = max(selected)
        check(f"A1e-4 has smallest max odd fraction {grid}",
              max_odd[1e-4] < max_odd[1e-5] < max_odd[1e-6], max_odd)

    figure_names = (
        "Figure_1_normalized_harmonic_development",
        "Figure_2_physical_and_normalized_mode_selection",
        "Figure_3_principal_positive_advective_pulses",
        "Figure_4_grid_sensitivity_by_amplitude",
    )
    for name in figure_names:
        for suffix in (".pdf", ".svg", ".png", ".alignment.json", ".collision_audit.json", ".pdf_text_audit.json"):
            check(f"figure artifact {name}{suffix}", (FIGURES / f"{name}{suffix}").is_file())
        alignment = json.loads((FIGURES / f"{name}.alignment.json").read_text())
        collision = json.loads((FIGURES / f"{name}.collision_audit.json").read_text())
        text_audit = json.loads((FIGURES / f"{name}.pdf_text_audit.json").read_text())
        check(f"alignment pass {name}", alignment.get("verdict") == "PASS", alignment.get("verdict"))
        check(f"collision pass {name}", collision.get("verdict") == "PASS", collision.get("verdict"))
        font_pass = bool(text_audit.get("auditable")) and int(text_audit.get("below_minimum_count", -1)) == 0
        check(f"font audit pass {name}", font_pass,
              {"minimum_found_pt": text_audit.get("minimum_found_pt"),
               "below_minimum_count": text_audit.get("below_minimum_count")})

    for required in ("REPORT_PHASE_A_ZH.md", "FIGURE_CONTRACT.md", "QA_NOTES.md", "REPRODUCE.md",
                     "PHASE_A_DECISION.json", "phase_a_protocol_and_summary.json"):
        check(f"delivery file {required}", (BASE / required).is_file())

    failures = [item for item in checks if not item["passed"]]
    result = {
        "phase": "WP14 Phase A",
        "delivery_passed": not failures,
        "checks_count": len(checks),
        "failures": failures,
    }
    (BASE / "phase_a_acceptance.json").write_text(json.dumps(result, indent=2))
    print(json.dumps(result, indent=2))
    if failures:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
