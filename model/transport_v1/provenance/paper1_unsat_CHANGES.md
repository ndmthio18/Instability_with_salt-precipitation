# kiemle2023b_unsat — the UNSATURATED fork (provenance & change record)

## Why this folder exists

The repository previously ran both papers' cases through **one** executable
(`code/dumux_kiemle2023b/.../test_model_salt`), with `Problem.FullySaturated` switching
between the Bringedal 2022 saturated model and the Bringedal 2025 partially saturated one.
That left one real coupling surface: any change made for the unsaturated reproduction
(instrumentation, the gravity parameter, constitutive switches) touched the same binary
that produced — and re-verifies — the frozen saturated golden results.

This folder decouples them:

| folder | role | status |
|---|---|---|
| `code/dumux_kiemle2023b/` | **saturated 2022** work; produced the frozen Tier-2 golden | **frozen** — do not modify for unsaturated needs |
| `code/dumux_kiemle2023b_unsat/` (this) | **partially saturated 2025** reproduction and all future unsaturated development | active |

The fork deliberately does **not** contain the `FullySaturated` switch or any of the
2022-specific physics (linear EOS, constant D, Sw=1 shortcut). It cannot run the
saturated model at all — misuse becomes an unknown-parameter no-op instead of a silently
wrong result.

## Origin (unmodified)

- Source: DaRUS **doi:10.18419/darus-4610** `kiemle2023b-main.tar.gz`
  (md5 `06da7b6c13c74f3e6e7948cfa76d9012`, verified against the DaRUS metadata; the
  tarball itself lives once in `../dumux_kiemle2023b/` and is shared).
- Built against the DUNE/DuMuX stack in `~/Dumux/dumux_salt_REV/dumux/` (referenced,
  not copied), CMAKE_BUILD_TYPE=Release, same compiler and flags as the sibling build.

## Our modifications (all in `kiemle2023b-main/appl/simpleRichardsnc/`)

Full unified diff: **`kiemle2023b_unsat_modifications.patch`** (applies cleanly to the
pristine tarball; verified). Only `main.cc` and `model.hh` change. The set is the
"minimal unsaturated reproduction kit":

1. **`Problem.GravitationalAcceleration`** — one parameter used *consistently* by the
   hydrostatic initial/boundary profiles and the momentum balance. Default **9.81** =
   pristine upstream = what produced the paper's numerical tables. (The paper's own two
   legs disagree: Table 5.1 prints 9.8, the LSA driver implies 9.80015, the DuMuX code
   uses 9.81. See `../paper_comparison/UNSAT_NOTES.md` §7.1.)
2. **`Problem.KrwMualemExponent`** (default `false` = upstream `L`-as-factor). `true`
   gives the paper's eq. (2.5) `S_e^L`-as-exponent form, which is what the authors' own
   LSA uses. Measured effect on the base Table 5.9 case: onset 33250 → 19150 s (−42 %),
   flipping the DuMuX/LSA ratio across 1. See UNSAT_NOTES §7.5.
3. **`PerturbationType = Eigenprofile`** + `Problem.EigenprofilePath` — runtime tabulated
   A(y) replacing upstream's recompile-per-case commented-out polynomials (which are kept
   verbatim in the source as the validation anchor, main.cc "Option C").
4. **`Problem.RandomSeed`** (default 1 = upstream) and `PerturbationHeight` honoured by
   the `Random` type (0.0 recovers upstream whole-domain behaviour exactly).
5. **Onset instrumentation** — `<name>_sigmatop.csv` written every step (validated to
   3e-12 against VTU extraction) and `TimeLoop.OutputInterval` VTU throttling (default 1
   = upstream). Without throttling a Table 5.9 run writes ~56 GB of VTU.
6. **Runtime similarity parameters** — `Problem.SaltDiffusivity`,
   `ReferenceViscosity`, `PureWaterDensity`, `DensityLinearCoefficient`,
   `DensityQuadraticCoefficient`, `MolarMassH2O`, and `MolarMassNaCl`. Their defaults are
   the pristine author-code values, so old reproduction inputs are unchanged. The
   `surrogate_consistent_v1` campaign writes the Bringedal-equation values explicitly,
   making the DuMuX density, `R_E`, `gamma`, and diffusion scaling agree with the LSA.

Removed relative to the sibling folder (i.e. never added here): the `FullySaturated`
switch and its four gated branches, and the member flag.

## Verification of the fork

- **Fork vs shared binary, same unsaturated input, 20 steps**: `x_top_mean` bit-identical;
  `sigma_top` max relative difference **1.6e-10**; first VTU identical, later VTUs differ
  at ulp level. Cause: deleting the (dead) `FullySaturated` branches changes code layout,
  and `-O3 -march=native` re-schedules floating-point ops (FMA contraction). Compiler
  effect, not physics.
- **Gold standard (completed 2026-08-23)**: the fork reproduces the paper's Table 5.9
  base-case eigenprofile onset **exactly — 33250.0 s vs the paper's 33250 s** (authors'
  A(y), g = 9.81, TEnd = 60000 s), and its full sigma_top trajectory matches the shared
  binary's to **2.4e-10** relative over the common window. The shorter TEnd leaves only a
  1.04x rise past the minimum, but the full-length shared-binary run already established
  the same minimum with an 8.7x rise, and the two trajectories are identical to float
  noise, so the identification carries over.
- **Runtime-parameter extension (2026-08-24)**: the tracked patch applies cleanly to the
  pristine tarball; the fork recompiles; a generated input varying every hydraulic and
  similarity parameter completes a real 8x4-cell DuMuX smoke run and writes the onset
  diagnostic. Python-side evaluation of the same input matches all target dimensionless
  groups and the full initial saturation profile below `2.3e-16`. An old input that omits
  every new runtime key reproduces the saved first 100 s of the paper trajectory
  bit-for-bit, verifying backward-compatible defaults.

## How to reproduce / revert

- From scratch: extract `../dumux_kiemle2023b/kiemle2023b-main.tar.gz` here, apply
  `kiemle2023b_unsat_modifications.patch` (`patch -p0` from inside the extracted root's
  parent), then configure with
  `cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_PREFIX_PATH="<all ~/Dumux/dumux_salt_REV/dumux/*/build-cmake, ;-separated>"`
  and `make test_model_salt` in `build-cmake/appl/simpleRichardsnc/`.
- The driver `code/paper_comparison_2025/run_case_2025.py` points here by default;
  `KIEMLE_EXE` overrides.
