# Phase A numerical and figure QA

- Nine frozen cases loaded; each has 7,201 accepted scalar/flux states and 145 saved full-field states through 6T.
- Source manifest contains 54 files with SHA-256, sizes and repository-relative paths.
- Generated tables contain 1,305 normalized field rows, 27 registered milestone rows, 63 modal checkpoints, 36 depth-specific pulse summaries and 60 amplitude-separated grid-difference rows.
- All four multi-panel figures passed the 1.5 pt rendered alignment gate.
- All four PDFs passed the 5 pt glyph-floor audit; minima are 6.0, 6.0, 6.0 and 5.6 pt.
- All four final PDFs passed the strict rendered collision audit with zero FAIL and zero WARN.
- Panel-by-panel visual inspection confirmed readable labels, consistent grid colors, explicit blank/unresolved pulse panels and no clipped annotations.
- Static source validation reports one non-blocking uncertainty warning because the scientific terms “seed mode” and “surface mean” trigger an aggregate-data heuristic. These are deterministic single trajectories and spatial means, not seed/fold replicate averages; no uncertainty interval exists or was omitted.
- The 600 dpi TIFF files are reproducible intermediates. PDF/SVG/PNG and the plotting source are the lightweight tracked delivery.
