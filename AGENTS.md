# Paper 2 — independent research repository

This directory is a standalone Git repository. The parent `22Surrogate_model`
repository owns Paper 1 and the shared baseline, and ignores this entire directory.
Run Git commands from this directory or with an explicit `git -C` path.
Do not stage, commit or push changes in the parent as part of ordinary Paper 2 work.

## Scope and scientific rules

- Read `Option_A_Research_Plan.md`, `README.md`, `PAPER1_REUSE.md` and `PROGRESS.md`.
- Current scope: Phase 0 source audit, then Pilot-0 feasibility. Repository setup
  does not establish feasibility or authorize a production campaign.
- Paper 2 studies effective salt redistribution. Onset information inherited from
  Paper 1 is an input/reference, not evidence that effective transport occurs.
- Follow `analysis/onset_semantics_20261102/CONTRACT_ZH.md`: LSA `t_crit` is the global critical-mode crossing; `t_seed(0.10)` is the growth-qualified observable mode's own crossing (theta=.10, floor=.03, fallback enabled), never an amplitude-reaching time. Keep DuMuX `t_detect` and undefined transport `t_eff` separate. Use `scripts/onset_semantics.py` or normalized views; preserve raw source fields and censoring.
- Preserve correlated soil joint support, constitutive branches, perturbation
  protocol, units, failed runs and censoring. A solubility threshold crossing is
  not a precipitation-kinetics prediction.
- Review assumptions and evidence–claim alignment critically. Ask about material
  ambiguity before dependent work; state any necessary assumptions.

## Source and output ownership

- Treat `paper1_sources.lock.json` as the frozen dependency contract. Verify hashes
  before using data. A changed upstream HEAD does not silently update the lock.
- Treat Paper 1 data and source code as read-only references. New Paper 2 code,
  manifests and results belong here. Before modifying a shared model, establish
  an isolated versioned copy/checkout and record its upstream commit and patch.
- Never modify the frozen saturated DuMuX source for unsaturated work. The
  unsaturated baseline provenance is in its own `CHANGES.md` and patch.
- No scipy/torch dependency is assumed; do not introduce either merely to simplify
  the verified numpy solvers.
- Runs lasting more than a couple of minutes belong on a compute host. Pin BLAS
  threads to one for process pools, retain run provenance and pull results back.
- Raw VTK fields, databases and generated/external data stay outside Git. Keep
  small explicit manifests, diagnostics, code and final scientific deliverables
  in version control. New committed result directories should be explicitly scoped.

## Progress and Git

- Update this repository's `PROGRESS.md` after meaningful research changes.
- Do not call Notion. Do not change the parent Paper 1 progress percentage merely
  because Paper 2 starts. Status stays `In progress` until the user decides otherwise.
- Commit messages must not include AI attribution or `Co-Authored-By` trailers.
- Paper 2's user-specified remote is `git@github.com:ndmthio18/surrogate_model.git`.
  Use its own `origin` and `codex/paper2-main` branch; do not change Paper 1's remote.
