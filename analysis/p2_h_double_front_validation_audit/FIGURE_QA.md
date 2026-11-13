# Figure evidence and QA

All figures use actual accepted numerical outputs; analytical curves are explicitly distinguished. Deterministic computations have no statistical sampling interval. The uncertainty evidence is grid/time refinement, not an invented error bar. Source values and failed/censored prefixes are retained.

Font is Times New Roman, following the user's typography preference. The generic static checker only recognizes a sans-serif whitelist and flags this family; this is a documented typography override, not a claim that its whitelist check passed. Actual PDF font embedding, rendered glyph sizes, panel alignment and collision geometry are separately audited. PNG/PDF/editable SVG are the requested deliverables; the static checker's optional TIFF warning is inapplicable. Positive-time/positive-value guards are explicit on logarithmic axes.

## Panel contract

| Figure | Panels and distinct role | Data/inference limits |
|---|---|---|
| 1 | a accepted step; b weighted component; c Jacobian count; d smallest liquid cell | Internal time, not requested output time. Crosses mark actual failures. A large component is not itself a failed RMS norm. |
| 2 | a–b crust diffusion; c–d liquid diffusion; e–f contact concentration; left u=.6, right x008 | Matched mesh/time curves, zero time omitted on log axes; no one-hour exclusion. |
| 2b | a–b top, c–d bottom, e–f net speed in first two seconds | Explicit 0–2 s display subset; complete early history remains in Figure 2 and source tables. |
| 2c | a–b seed-gap sensitivity; c–d uniform versus stationary-concentration initialization | Different physical initial conditions, not a convergence comparison; initial inventories differ. |
| 3 | a top; b bottom; c actual thickness | Same final implementation, genuine moving geometry; no inventory-to-theory conversion. |
| 4 | a top speeds; b net speeds; c top discrepancy; d net discrepancy | Solid ALE and dashed analytical; full transient retained, stated analytical assumptions apply. |
| 5 | a u=.6; b u=.9; c x008: signed top, diagnostic-bottom and net discrepancies | Diagnostic inferred bottom law is not independently validated. |
| 6 | rows u=.6/.9; columns bottom speed, bottom position, actual thickness | kappa=.5/1/2 at identical resolution/time; no parameter selected by fit. |
| 7 | a volume mean; b core mean; c final profile; d bottom trace | Absent early core is NaN; theoretical constant never inputs the ALE. |
| 8 | top, bottom, net, thickness grid differences | 128→256 and 256→512, common saved times, nonzero velocity normalization. |
| 8b | top, bottom, net, thickness time differences | Two tolerances and additional 512-cell check. Exactly zero differences are omitted only on log ordinate; tables retain them. |

All panels require the final render-time alignment JSON, PDF collision JSON, PDF text-size audit and visual review. Earlier previews are superseded by the final render. Long labels and short-page clipping found in the first pass were corrected; no opaque annotation boxes conceal plotted data. Final QA summary is saved as `results/figure_qa.json`.

Final visual review on 2026-11-13: all 11 figures inspected in a contact sheet, with Figures 3, 4 and 7 additionally inspected at full resolution. No missing labels, clipped legends or misplaced units observed. Overlapping curves denote numerically close results; complete discrepancies are retained in tables. Figure 3c includes integrated analytical thickness from the same initial thickness. All 11 final PDFs passed glyph-size and collision checks; 7 numerical test groups and 241 delivery/provenance assertions passed.
