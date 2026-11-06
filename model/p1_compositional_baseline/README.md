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
