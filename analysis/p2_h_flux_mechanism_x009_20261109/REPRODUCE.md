# 复现

## 分支与工作目录

- Git branch: `research/p2-h-flux-mechanism-x009`
- Verified result commit: `82bfef1f747dc84e8db63c0ce719cbccbdcf773b`
- HPCC: `/WORK/hulm_work/zhixin/p2_h_flux_mechanism_x009_20261109`

## Slurm 命令

非交互 SSH 需先设置：

```bash
export PATH=/usr/local/bin:/usr/bin:/bin:/rmprog/slurm/v22.05.7/bin
cd /WORK/hulm_work/zhixin/p2_h_flux_mechanism_x009_20261109
```

实际提交链：

```bash
sbatch --parsable source/flux_x009_build.slurm
sbatch --parsable --dependency=afterok:16445808 source/flux_x009_runs.slurm
sbatch --parsable --dependency=afterok:16445813 repo/model/p2_h_ale_coupling/flux_x009_postprocess.slurm
```

对应作业为 `16445808`、`16445813`、`16445845`。输入、源文件、二进制和输出校验和分别保存在远端：

- `input_sha256_build4.txt`
- `run_input_sha256.txt`
- `unit/source_binary.sha256`
- `post_input_sha256_updated.txt`
- `analysis_output/output_sha256.txt`

修订后的最终后处理与紧凑打包作业分别为 `16446287` 和 `16446308`。紧凑归档 SHA-256 为 `e9f9481d1913e3bf2c8d1f5b2a73d7b1efa99d69fc1a4aa48c35aa32786437c3`。

## 后处理入口

```bash
python3 repo/model/p2_h_ale_coupling/analyze_flux_mechanism.py \
  --x008-progressive runs/x008_progressive_reservoir_n80_a64_dt0p5_diagnostic \
  --x008-adaptive /WORK/hulm_work/zhixin/p2_h_reservoir_x008_h3p5_20261108/diagnostic_bound \
  --x009 runs/x009_progressive_reservoir_n80_a64_dt0p5 \
  --x009-half runs/x009_progressive_reservoir_n80_a64_dt0p25_timecheck \
  --out analysis_output
```

完整原始场留在 HPCC。Git 只保存紧凑 CSV、图、报告、脚本、校验和与原始结果索引，不提交大体积逐步原始文件。
