# Reproduction instructions

## Repository state

- Branch: `research/p2-h-interface-concentration-audit`
- Base: `1fe9d1f7f9fee8094b671367992d463a8e8b96c9`
- HPCC root: `/WORK/hulm_work/zhixin/p2_h_interface_concentration_audit_20261108`
- Reused tested DuMuX solver SHA-256:
  `ad3188a97326887b1f980bfb96c40205697bc444e9bdf7028040bf761216e088`

Raw VTK/CSV fields remain on HPCC and are not committed because the project
`AGENTS.md` requires large fields to remain outside Git. The delivery contains
checksums, compact diagnostic tables, source, Slurm logs and figures.

## Registered jobs

- `16438804`: compile standalone ALE counterfactual driver.
- `16438812_[0-7]`: fixed-bottom-concentration and initialization counterfactuals.
- `16438805_[0-6]`: matched four-case, ALE-grid and soil-grid coupled audit.
- `16438821_7`: exact 0.1 micrometre activation replay.
- The postprocessing job ID is recorded in `provenance.json` after execution.

All numerical kernels, compilation and plotting ran on Slurm compute nodes.
Login-node actions were limited to staging, submission and read-only monitoring.

## Coupled cases

From the HPCC root:

```bash
sbatch interface_audit_coupled.slurm
sbatch --array=7 interface_audit_coupled.slurm
```

The case registry is in `run_interface_audit.py`. It freezes
`k_r=1e-3 m/s`, `alpha_av=1`, `epsilon_birth=0.43`, `Sw_min=0.99`,
`MaxDonorDepth=5 mm`, `dt=0.5 s` and the 3.5 micrometre activation threshold,
except for the separately named 0.1 micrometre historical replay.

## Standalone counterfactuals

```bash
build=$(sbatch --parsable interface_audit_build_counterfactual.slurm)
sbatch --dependency=afterok:${build} interface_audit_counterfactual.slurm
```

The driver writes every accepted 0.5 s state, the first rejected field if any,
the bottom salt-source ledger, and the mass-conservative initialization
adjustment ledger.

## Postprocessing

After all coupled cases terminate:

```bash
sbatch interface_audit_postprocess.slurm
```

The script independently reconstructs the ALE bottom operator and reconstructs
the DuMuX top flux from soil storage change plus bottom reservoir influx. It
then generates the required CSV files, figures, verification record and SHA-256
manifests. Failed runs are never interpolated beyond their last accepted time.
