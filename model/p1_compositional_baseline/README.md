# P1 conservative compositional baseline

This isolated DUNE module implements the P1-only isothermal H2O--air--NaCl baseline. It uses DuMuX `2pnc` with molar balances and contains no mineral, solid-NaCl or crust state.

## Build on the audited `sal` stack

Configure with the DUNE module build directories in `CMAKE_PREFIX_PATH`, `-DCMAKE_BUILD_TYPE=Release` and `-DCMAKE_DISABLE_FIND_PACKAGE_SuperLU=TRUE`, then build `p1_compositional_1d` and `p1_compositional_2d` with one compile job. Exact commits and executable hashes are in `manifests/dependencies_20261106.json`.

On the audited host, set the following before execution because the default OpenMPI/UCX path deadlocked in `MPI_Init`:

```bash
export OMPI_MCA_pml=ob1
export OMPI_MCA_btl=self,tcp
export OMPI_MCA_mtl=^ofi
export UCX_TLS=tcp,self
export OMP_NUM_THREADS=1
export TBB_NUM_THREADS=1
```

Run V1--V7 with:

```bash
python3 verification/run_verification.py \
  --binary-1d /path/to/p1_compositional_1d \
  --binary-2d /path/to/p1_compositional_2d \
  --output /new/versioned/output
python3 verification/analyze_verification.py /new/versioned/output
```

`run_verification.py` refuses to reuse an output directory. The acceptance criteria are stored separately and hashed before execution.

## Formal Case 36877 inputs

`benchmark/case36877/B0_production.input` and `B1_production.input` preserve the 0.60 m by 0.20 m domain and 600 by 200 target grid. The B1 input references the already frozen repository eigenprofile whose SHA-256 is recorded in `reference.json`.

## Accepted evidence

- `results/p1_verification_20261106/`: final build-07 V1--V7 assessment and provenance.
- `results/p1_benchmark_20261106/`: bounded B0/B1 exploratory comparison, figures and checksums.
- `results/p1_convergence_20261106/`: registered one-factor spatial/time-step assessment, figures and checksums.
- `../../docs/P1_GO_NO_GO.md`: scientific decision and stopping rule.

The benchmark result is a right-censored non-detection of instability, not an onset at 60,000 s. The formal 600 by 200 inputs are preserved but were not executed in P1; the bounded benchmark used 280 by 100 cells. Consult the decision report before reusing these results.
