# WP14 Phase B QA notes

- The plotting source passed the strict static preflight with 21 PASS, 0 WARN and 0 FAIL.
- Figures 5--7 passed the 1.5 pt panel-alignment gate.
- PDF effective text-size audits passed at the 5 pt floor. The minimum effective sizes were 5.6 pt for Figures 5 and 7 and 8.0 pt for Figure 6.
- Strict rendered collision audits passed with 0 FAIL and 0 WARN for all three figures.
- All three final PNGs were inspected at high resolution after the final render. Curves, panel labels, legends/direct labels, scales and axis text are legible; no visible clipping or unintended overlap was found.
- The two high relative-difference points in Figure 7a are the mirror-symmetry and odd-mode diagnostics. Their absolute values remain extremely small; the final report explicitly prevents their large relative percentages from being read as a macroscopic modal or transport branch change.
- The shaded `0.25 <= Delta_t/Delta_x <= 4` region is a registered interpretation band for observed finite-resolution differences. It is not a confidence interval or uncertainty distribution.
- Transient 600-dpi TIFF sources and collision-overlay PDFs were generated for QA and then removed. PDF, SVG, PNG and machine-readable audit records are retained.
