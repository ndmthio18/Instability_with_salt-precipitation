# WP11 pre-solver submission failure record

## Scope

The first Slurm submissions for WP11 were:

- R48-M1: job `16305507`
- R96-M1: job `16305508`

Both jobs ran on `ibc11b11n34` and failed with exit code `1:0` after six seconds on
2026-11-05. The solver executable was never invoked: both
`production_started.json` sentinels remained zero bytes, both progress records remained
empty, and neither run produced an accepted state.

## Root cause

The immutable A28 profile was staged once at the WP11 durable root as
`eigenprofile.dat`. The worker incorrectly attempted to copy it from the grid-specific
`R48/` and `R96/` directories. Both jobs therefore stopped at the pre-solver staging
step with `FileNotFoundError`.

The worker has been corrected to copy `params.input` from the grid directory and the
shared `eigenprofile.dat` from the WP11 root. No production resubmission is authorized
by this correction itself. In accordance with `command/18.md`, a new submission must
wait for explicit user authorization; there is no restartable numerical state.

## Preserved evidence

The subdirectories `R48/` and `R96/` contain the scheduler stdout, stderr and worker
failure JSON copied from durable storage. The zero-byte stdout files are preserved as
evidence that the failure occurred before solver output.
