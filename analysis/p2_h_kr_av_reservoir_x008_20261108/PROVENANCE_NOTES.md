# Provenance notes

`provenance.json` records commit `cdf9a790...`, the commit present when the
diagnostic package and figures were generated. The physical pilot itself ran
the registered scenario at `b9815894...`; the only intervening commits added
diagnostic capture and censored postprocessing and did not alter the physical
model or pilot inputs.

The solver executable was reused from the immediately preceding reservoir
revision because the x=0.08 change is a runtime parameter. Its SHA-256 was
verified as `ad3188a97326887b1f980bfb96c40205697bc444e9bdf7028040bf761216e088`.

