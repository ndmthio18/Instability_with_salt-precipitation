"""Prepare the explicitly authorized A=1e-5 three-grid campaign."""
from pathlib import Path

import wp12_prepare_amplitude as base


base.OUT = base.PAPER2 / "analysis" / "wp13_amplitude_1e5_20261106"
base.AMPLITUDE = 1.0e-5
base.CAMPAIGN_TAG = "WP13"
base.CAMPAIGN_NAME = "WP13 intermediate-amplitude three-grid convergence test"
for grid, specification in base.SPECS.items():
    specification["name"] = f"WP13-{grid}-A1e5-6T"


if __name__ == "__main__":
    base.main()
