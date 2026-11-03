# command/15 交付索引：唯一一次连续修正严格背景（0 → 4T）

**结果：运行正常结束，4800/4800 接受步到达 4T，48 项强制验收全部通过，
4T 内未触及局部溶解度阈值（最大 x = 0.056969，占 0.1 的 57.0%）。**
修正求解器通过了旧 v1 的失败窗口（第 4414 步，3 次迭代接受）。
未重试、未重启、未回放、未中途调整容差；R96、追加二维、多波长与 8T 均未启动。

| 文档 | 内容 |
|---|---|
| [执行与强制验收](EXECUTION_AND_ACCEPTANCE_ZH.md) | 启动前准入、运行成本、48 项验收逐项数值、失败预案 |
| [与历史严格背景的对照](HISTORICAL_COMPARISON_ZH.md) | 剖面、七面库存、带符号通量、累计与闭合、Newton 诊断、旧失败窗口 |
| [科学审查](SCIENTIFIC_REVIEW_ZH.md) | R24/R48 的 M0/M1 相对背景的响应、七面再分配、表层均值与局部最大、网格敏感性、能说与不能说 |
| [图 A](Figure_A_background_validation.png) | 背景验证：与历史差异、Newton 成本、旧失败窗口、逐步守恒闭合 |
| [图 B](Figure_B_background_referenced_response.png) | 科学结果：表层均值/局部最大、零模态响应、网格敏感性、七面再分配 |

数据产物：`run/`（输入、日志、事件、七面预算、收敛记录、归档回执）、
`validation.json`、`historical_comparison.json`、`checkpoint_comparison.csv`、
`scientific_analysis.json`、`background_referenced_responses.csv`、
`seven_plane_redistribution.csv.gz`、`surface_mean_versus_local_max.csv`、
`grid_sensitivity_vs_background.csv`、`cost_audit.json`、`acceptance_receipt.json`、
`delivery_verification.json`、`provenance.json`、`numerical_tests.txt`。

原始场（4801 帧，压缩 59.1 MiB）按仓库约定不入 Git，留在 frodo
`/temp/zhixin/paper2_wp8_full4T_20261103/run/fields.tar.gz` 与本地忽略缓存
`tmp/wp8_full4T_20261103/run/`，两端逐帧 SHA-256 已核对（4801/4801）。

## 三个一句话结论

1. **背景可用。** 与历史接受前缀 4414 个共有时刻的全部差异都在舍入量级
   （max|Δx| = 9.02e-16，max|Δp| = 3.06e-10 Pa），低于每一项原门槛数个量级。
2. **零模态响应在 0→3T 上网格不敏感。** R24 与 R48 推出的响应之比为 1.003/1.009/1.015
   （T/2T/3T），显著高于 9.02e-16 的可分辨下限。
3. **4T 的响应幅度尚未收敛，不可引用。** 3T 之后比值升至 2.74（3.6775T）与 382（4T），
   且本应等于背景的无扰动 M0 自身在 4T 偏离 1.24e-06（R24）对 2.90e-09（R48）。

## 复算

```bash
cd paper/2/scripts
python validate_wp8_full4t.py
python compare_wp8_full4t.py
python analyze_wp8_full4t.py
python plot_wp8_full4t.py
python report_wp8_full4t.py
python verify_wp8_full4t_delivery.py
```

以上脚本都不启动求解器。`wp8_full4t_campaign.py` 的 `prepare`/`launch` 对已存在的
运行目录会断言失败，因此不提供第二次 4T 的启动入口。
