# 复现说明

代码分支：`research/prat-analytical-growth-law-comparison`

基线：`c6f24b979999a9ae87602cfa1c6135597c82cbcb`

HPCC 工作目录：`/WORK/hulm_work/zhixin/prat_growth_law_comparison_20261109`

所有编译、数值计算和绘图均通过 Slurm 计算节点执行。典型命令为：

```bash
cd /WORK/hulm_work/zhixin/prat_growth_law_comparison_20261109
sbatch repo/model/prat_analytical_boundary/build.slurm
sbatch repo/model/prat_analytical_boundary/growth_comparison.slurm M1 0.5
sbatch repo/model/prat_analytical_boundary/growth_comparison.slurm M2 0.5
sbatch repo/model/prat_analytical_boundary/growth_comparison.slurm M3 0.5
sbatch repo/model/prat_analytical_boundary/regression.slurm
```

只有 M2/M3 的首轮计算满足条件时，才提交 `dt=0.25 s`：

```bash
sbatch repo/model/prat_analytical_boundary/growth_comparison.slurm M2 0.25
sbatch repo/model/prat_analytical_boundary/growth_comparison.slurm M3 0.25
```

最后在所有必需结果就绪后运行：

```bash
sbatch repo/model/prat_analytical_boundary/growth_analysis.slurm
```

实际 Job ID、节点、命令参数、退出码、二进制及原始历史校验和见 `records/`、`provenance.json` 与 `checksums.json`。完整接受步历史和求解日志保留在上述 HPCC 工作目录；失败曲线不外推。

本次实际 Job ID：构建 `16466140`；M1/M2/M3 (`dt=0.5 s`) 为 `16466165/16466166/16466167`；M2/M3 (`dt=0.25 s`) 为 `16466175/16466176`；最终 M1 回归 `16466174`；后处理 `16466190`。失败的回归审计尝试 `16466168/16466173` 保留在 HPCC 日志中。
