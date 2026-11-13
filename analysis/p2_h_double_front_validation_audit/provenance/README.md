# Provenance

Primary remote directory: `/WORK/hulm_work/zhixin/p2_h_double_front_audit_20261112`. All executable builds and simulations use Slurm compute nodes, GCC 12.2.0 and the existing SUNDIALS 6.7 dependency installation. Slurm allocation/state/exit/node/time records accompany raw source/data SHA-256 manifests.

The exact Task 27 reproduction is in source_legacy. Development versions are separate frozen source directories and executable names: smoke_v0, double_v1, longdouble_v2, early_v3, nonlinear_v5 and setup_v6. Output IDs identify each run. Retried outputs never replace previous failed prefixes.

The intermediate every-step fresh-Jacobian experiment did not receive a separate frozen v4 source directory before the next mode-selection branch was added. Its logs and executable-run evidence are retained; it is not advertised as an exactly pinned final implementation. Its FRESH_JAC=1 execution path is also represented by the later source, but that is a source-equivalence statement, not proof of an identical binary. The final setup_v6 source/executable hash pair is pinned separately and was used for the first matched production and kinetic comparisons. The final primary matrix was replayed with the separately pinned serial sqrt_v7 source/executable; see the implementation-bridge table and primary case registry.

Failed CVODE attempts have status, diagnostics and error-weighted state snapshots. The separate high-D stress test has a different physical test problem and must not be combined with reference-material production rows. An intermediate SKIP_STRESS run tests only remaining checks; successful final setup verification runs the entire original suite.

Large full raw archives remain on HPCC with SHA-256 manifests. Local derived compressed tables and selected profile files are traceable to those raw hashes. Partial copies and extraction workspaces are ignored by Git; they are not published as complete evidence.


Two early high-D stress failure snapshots originally shared the generic diagnostic basename. Their original console logs and frozen source versions were retained, and both stress cases were reproduced under distinct `stress_repro_0`/`stress_repro_1` IDs to preserve separate weighted-state evidence. Reproductions are labelled as new runs, not falsely represented as the original overwritten snapshots.

The optional parallel-v8 source and logs are retained separately. Its short 1/8-thread and v7 histories were byte-identical; its duplicate 512-cell trajectory was intentionally stopped after a common-prefix equality audit, with a distinct censored status. All primary physical comparisons use the serial sqrt_v7 executable.


The first stage pack reported `tar: raw: file changed as we read it` because another diagnostic run added files to that directory. It is a development snapshot, not a certified complete archive. Completed primary files are checked against final per-file hashes; the final quiescent full-raw archive supersedes this snapshot. A later SHA-256 of the unchanged stage archive verifies transfer integrity only, not original completeness. The early kinetic checksum file had a basename-label typo; its digest verified the explicitly selected kinetic archive, and the script label is corrected for future packs.
