# Failure and refusal register

No numerical failure was selectively removed or automatically rerun.
See `deliverables/audit/case_summary.csv` for every case and
`deliverables/audit/verification.json` for accepted times and reasons.

| Class | Cases | Outcome |
|---|---|---|
| Budget stop before activation | All twelve `new_n320_h{1,2,5,10}_d{2,5,10}`, `old_n320`, `pre_n320`, `dt_half_n320` | Exit 1; no attempted activation; no final 2500 s state. Boundary-area attribution is in REPORT_ZH.md. |
| Adaptive depth refusal | `new_n40_h10_d2`, `new_n80_h10_d2` | Exit 0 at 2500 s as precrust only; insufficient admissible liquid; not coupling success. |
| Old saturation refusal | `old_n160`, `old_n80_h5`, `old_n80_h10` | Exit 0 at 2500 s as precrust only. |
| Diagnostic precrust controls | `pre_n40`, `pre_n80`, `pre_n160` | Deliberately no activation. Their completion is not a coupling pass. |
| Deliberate depth refusal | Controlled `insufficient_depth`; real `refusal/depth_refused` | No soil/crust transfer. Real transport history equals `refusal/precrust_reference`. |
| Injected trial rejection | `rollback_activation`, `rollback_coupled` | One rejected trial each; accepted states conserve inventory; both continue to 2500 s. |

Primary job **16371535 TIMEOUT** interrupted seven unfinished runs. The exact
archived paths, accepted times, original commands and replay order are in
`deliverables/resume_metadata.json`. Job **16371664 COMPLETED** replayed only
those interrupted runs, from the identical initial state with the same binary
and command. Numerical exit-1 results were not retried. Archived prefixes are
in HPCC `failures/16371535/` and the hashed full-history archive.

Audit job **16371564 FAILED** before final verification because a frozen-source
manifest expected a project metadata file absent from the historical staging.
Its script and output are retained at HPCC `failures/audit_16371564/`; its log
is also in `deliverables/audit-16371564.out`. The corrected audit hashes the
complete actual soil application plus protected ALE/BDF/film/Prat dependencies.
Job **16371746 COMPLETED** with no integrity errors. This operational correction
does not turn the failed scientific convergence gates into passes.

Figure QA jobs **16372772** and **16372774** only rerendered plots. They did not
run the solver or modify numerical histories. The second revision also uses
common donor-mass axes and removes overlapping initial-saturation ticks.

The 97 MB full-history archive remains on HPCC with a SHA-256 in provenance.
An optional local download was stopped after the compact evidence had been
verified; the local `.partial` is not an archive and is not used as evidence.
