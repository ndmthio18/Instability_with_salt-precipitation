# Task 28 — double-front numerical and physical audit

Entry point: [FINAL_REVIEW.md](FINAL_REVIEW.md). Machine-readable conclusions: `results/validation_status.json`. This is an independent geometric comparison: ALE thickness is zt−zb, never crystal inventory divided by the analytical core fraction.

## Scope and baseline

Branch: `research/p2-h-double-front-validation-audit`; exact parent: `04f9bb2d856512531ebcb6265e4b46ef35cc43b0`; remote: `instability`. All Task27 code/results remain intact. The soil flux is an external input from the previously computed x008 column, not a new coupled-soil calculation. x008 output is restricted to 160951.5 s (44.70875 h) after the forcing origin. At CVODE internal overshoot the endpoint forcing is held constant for interpolation; this is recorded rather than described as a new soil prediction.

## Reading order

1. `NUMERICAL_FAILURE_DIAGNOSIS.md`: reproduced failures, exact conserved-variable transformation, error control, Jacobian and nonlinear-solve ablations.
2. `INTERFACE_FLUX_CONVERGENCE.md`: actual early relaxation, internal steps, grid and initial-condition comparisons, earliest sampled reliable window.
3. `BOTTOM_KINETICS_AUDIT.md`: independent moving-frame jump and region-inventory budgets; additional exposed-area assumption.
4. `TOP_AND_NET_GROWTH_VALIDATION.md`: separate top/net hypotheses, error cancellation, geometric identity.
5. `CORE_SOLID_FRACTION_AUDIT.md`: birth/core/mean/trace distinctions and limits of causal attribution.
6. `FINAL_REVIEW.md`: quantified conclusions, eight gates and remaining limitations.

`PROTOCOL.md` states numerical targets and comparison rules. `FIGURE_CONTRACT.md` and `FIGURE_QA.md` document plot semantics, typography and rendered checks. Full accepted histories are stored per case in `results/histories` as lossless CSV.gz, including failed prefixes; summary tables live in `results`; figures are PNG/PDF/editable-text SVG in `figures`.

## Reproduction

All C++ builds and simulations run on HPCC Slurm compute nodes. Use the registered dependency installation and GCC 12.2.0 in the checked-in sbatch files. These scripts contain the actual execution directory; choose a new directory and distinct run IDs for an independent reproduction, never overwrite previous raw evidence.

`reproduce.sbatch` recreates the original failures. `reproduce_stress.sbatch` independently recreates the intermediate double/extended-precision stress failures with distinct diagnostic names. `build_setup.sbatch` verifies the corrected solver; `build_sqrt.sbatch` verifies the arithmetic optimization. Final runs use frozen `source_sqrt_v7`, `run_case_sqrt_v7` and `FRESH_JAC=setup`. The only v6→v7 change evaluates eps^(3/2) as eps*sqrt(eps), including the matching derivative; it changes no physical equation or tolerance. The implementation bridge is explicitly tested.

`run_fast_p01.sbatch` performs the complete P0/P1 replay. `gate_fast.sbatch` must pass before `run_final.sbatch`; the latter also requires the previously audited P2 balance/sensitivity gate. `final_cases.csv` pins all primary material/grid/tolerance settings. `pack_final.sbatch` archives hashes and selected output only after the final array exits. Full raw profiles and failed attempts remain remotely archived; local source-data tables and selected profiles are traceable to their SHA-256 manifests.

For full recomputation, first retrieve `task28_complete_raw.tar.gz` from the remote directory given in `provenance/README.md`, verify its SHA-256 against `provenance/task28_complete_raw.sha256`, and extract its `raw/` directory here. The archive is intentionally outside Git. Already computed histories in `results/histories/` and all delivered figures can be inspected without the raw archive. Exact local package versions are in `provenance/local_analysis_environment.json`.

Local postprocessing uses Python with numpy, pandas, matplotlib and PyMuPDF:

```powershell
python analyze.py
python analyze_final.py
python test_audit.py
python make_figures.py
python supplement_early.py
python supplement_time.py
python write_final_report.py
python qa_figures.py
python audit_delivery.py
```

The plotting source imports the included panel-alignment auditor. Numerical tables and final reports must be regenerated after changing any source data. Conservation checks, numerical convergence, kinetic sensitivity and physical validity remain separate conclusions.
