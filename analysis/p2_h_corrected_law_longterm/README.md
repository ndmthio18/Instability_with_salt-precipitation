# Long-term corrected-law applicability study

Task 26; branch research/p2-h-corrected-law-longterm; baseline
c6fb6277c4cdf1689d80f4bbdda7a421060d874a.

Run `python prepare.py` to reproduce the registered P0 matrix. Production
build and simulations use build.sbatch and run.sbatch under
/WORK/hulm_work/zhixin/p2_h_corrected_law_longterm_20261112 on HPCC cnmix.
The shared operators are staged from model/p2_h_ale_coupling; source and
executable hashes are recorded in build_sha256.txt.

Dependencies: GCC 12.2.0; SUNDIALS 6.7.0; local Python 3.12.4,
NumPy 2.0.0, SciPy 1.14.0, pandas 2.2.2, Matplotlib 3.9.2.
With raw histories/profiles/status files retrieved under raw/, run
`python analyze.py`. It writes compact tables and eight PNG/PDF/SVG figures.
Raw histories are ignored by Git; the prior x008 forcing CSV is a read-only
dependency already in the baseline. The task never overwrites prior results.

## Completed campaign

70 registered trajectories: 33 P0, 24 P1, 8 P1 event diagnostics and 5 P1
adaptive spatial checks. These counts include the eight original state-bound
failures, retained separately from the event reruns. A successful scheduler
exit at a physical screen does not mean that the requested target was reached.
The scientific findings and eight required answers are in
LONGTERM_VALIDATION_REPORT.md; final per-case decisions are in results/.

Generate the matrix with prepare.py, prepare_p1.py, prepare_events.py and
prepare_refine.py. Execute only on compute nodes with the matching array
scripts. The first P0 array used rows 0--29; the additional P0 refinement
used rows 30--32 via `sbatch --array=30-32 run.sbatch`. Build versions and
job mappings are in execution_record.json. Historical driver versions are
preserved in provenance/ and full remote source snapshots.

After retrieving histories, selected profile rows and status files into raw/:

```text
python soil_audit.py
python analyze.py
python summarize_physics.py
python audit_results.py
python -m unittest -v test_models.py
```

Run these commands from this directory. Committed history.csv.gz and
profiles.csv.gz are lossless compressed CSV tables (pandas reads them directly).
history contains all outputs through 3600 s and every 600 s afterward;
profiles contain exact selected raw rows. Full scalar histories and all
production profile outputs remain under the HPCC run directory and in
evidence_final.tar.gz. Compact transport archives are ignored by Git;
their checksums and the full archive checksum are committed. final_sha256.txt
refers to FULL remote originals, not locally selected profile files.

Figures use editable Times New Roman SVG text, embedded-font PDF and 250 dpi
PNG. Captions and grid choices are described in the report. Figure 6 represents
sampled endpoints, not a continuous universal validity boundary.

No shared numerical operator, earlier study, or forcing history was changed.
No phi_core is used to convert inventory to a claimed geometric height.
