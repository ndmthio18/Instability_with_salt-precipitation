# Running-case snapshots

This calculation was still running when the first research commit was prepared.
The registered input, plan, screening and launch provenance are tracked directly.
Automatically refreshed local reports and diagnostics are excluded from Git at
their live paths; immutable snapshots under `snapshots/` retain the captured
report, metrics, diagnostics and hashes. A snapshot is not a completed result.

After completion, collect and verify the final evidence, create a new snapshot,
and update `PROGRESS.md`. Do not overwrite the previous snapshot. The live
monitor and solver continue independently of Git operations.
