#!/usr/bin/env python3
"""Fail if the P1 source accidentally acquires precipitation/mineral model code."""

from pathlib import Path


ROOT = Path(__file__).resolve().parents[1] / "source" / "appl" / "p1"
FORBIDDEN = ("2pncmin", "mineralization", "solidnacl", "precipitationsource")


def main() -> int:
    findings = []
    for path in ROOT.glob("*.*"):
        text = path.read_text(encoding="utf-8").lower().replace("_", "").replace(" ", "")
        for token in FORBIDDEN:
            if token in text:
                findings.append(f"{path.name}: {token}")
    if findings:
        raise SystemExit("P1 scope violation: " + "; ".join(findings))
    print("PASS: no mineral/precipitation model token in P1 application source")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
