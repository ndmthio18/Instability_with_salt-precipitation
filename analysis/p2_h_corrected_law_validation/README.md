# Independent validation of the corrected finite-Da growth law

This directory implements task 25. Its primary observable is the ALE top
boundary velocity, not an inventory-equivalent height.

## Frozen baseline

- Branch point: eaa5156420000fffe35d854ece86b52c0bbbe44a
- Material parameters: configs/manifest.json
- Constant evaporation: 1.78e-5 kg m-2 s-1
- Main matrix: six target \(Da\) values by five \(u\) values
- \(Da\) is varied only through initial crust thickness; \(Da^*\) is fixed.

## Boundary modes

- **neumann**: prescribe the downward diffusive salt flux, reconstruct the
  bottom face concentration from that flux, and use the reconstructed
  concentration for advective inflow.
- **equilibrium_inflow**: preserve the previous x008 convention
  \(J_{s,b}=C_{\mathrm{eq}}F_b-B_{\mathrm{soil}}\) while separately reporting
  the Neumann-reconstructed face concentration.
- **dirichlet**: use the original finite-volume \(C_b=C_{\mathrm{eq}}\)
  reservoir boundary and diagnose its resulting diffusive flux.

The main validation uses **neumann**. The other modes are boundary-condition
audits and must not be mixed with the main normalized-velocity comparison.

## Local Stage A

    python -m unittest discover -s analysis\p2_h_corrected_law_validation -p "test_*.py" -v
    python analysis\p2_h_corrected_law_validation\verify_bvp.py

Postprocessing and all result-level acceptance tests:

    cd analysis\p2_h_corrected_law_validation
    python analyze_validation.py
    python -m unittest discover -p "test_*.py" -v

Local Python versions and the remote compiler/SUNDIALS environment are frozen
in `execution_record.json`. The exact 135-row matrix is in
`configs/cases.csv`; derived values and physical lengths are in
`configs/manifest.json`.

The BVP solver works in \(\zeta=z/\ell\) to avoid dimensional residual scaling.
It starts from a zero profile and is not calibrated against the closed form.

## HPCC stages

The production run directory is:

    /WORK/hulm_work/zhixin/p2_h_corrected_law_validation_20261112

The build script compiles and tests the operator. The array script executes one
row of configs/cases.csv per Slurm array task. Both scripts refuse computation
outside Slurm or on login/test nodes.

Successful jobs were build `16569420`, 135-case array `16570915`, and pack
`16571056`. Pre-solver CRLF/wrapper failures and the failed first pack are
retained in the remote run directory and listed in `execution_record.json`.

## Outputs

- `VALIDATION_REPORT.md`: scientific findings and the seven required answers.
- `results/`: compact BVP, ALE, convergence, boundary and x008 tables.
- `figures/`: six figures in editable SVG, PDF and PNG plus source manifests.
- `analyze_validation.py`: deterministic table/figure regeneration.
- `test_validation_results.py`: numerical acceptance tests.

The raw 135-case histories and profiles are intentionally not committed. They
remain in the checksummed remote archive and the ignored local `raw_remote/`
directory; the archive SHA-256 is recorded in `execution_record.json`.

## Evidence policy

Raw case histories and profiles are retained remotely and transferred as a
checksummed archive. Compact derived tables, scripts, tests, reports, and
publication figures are committed. Failed or physically inadmissible cases are
kept with their status and stderr; they are not silently replaced.

## Interpretation boundary

Agreement of the analytical BVP with solve_bvp establishes mathematical
correctness only. Agreement with the transient ALE top velocity additionally
requires a quasi-steady window, boundary consistency, and spatial/time
convergence. Cases outside those conditions remain outside demonstrated
validity even if their numerical solver completes.
