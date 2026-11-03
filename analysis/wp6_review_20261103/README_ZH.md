# WP6 command/8 复核交付

日期：2026-11-03。基线提交：`6e8e03428779dd1f64134cc5878c9f9ebb03e927`。

**结论：本轮复核完成，当前 4T 为 NO-GO。** 异步停止分析已经修复并验证；原 2T 配对结果逐字段重现。背景差异已确认为随时间发展、受时间步影响的数值扰动，但单个端点指标无法保证全场对称性，尚不能唯一分离求解器、边界离散和噪声放大的贡献。原失败筛查保留，未放宽阈值。

只新增了一条获准的 M0 dt/2、2T 对照：frodo，1 核，2,400 接受步，求解墙钟 310.24 s；独立收支通过，2,401 帧回传并核验 SHA。没有启动 4T、8T、完整域或 nx60，也没有更改 Paper 1、WP5 或 command/7 基准。

| 工作 | 交付与状态 |
|---|---|
| 背景诊断 | [BACKGROUND_ZH.md](BACKGROUND_ZH.md)，原 1,201 帧及对照 2,401 帧全场、时间序列、空间与模态分析 |
| 异步停止 | [ASYNC_STOPPING_ZH.md](ASYNC_STOPPING_ZH.md)，22 项测试及 9,608 行原结果零差异回归 |
| 七面输运 | [TRANSPORT_ZH.md](TRANSPORT_ZH.md)，[seven_plane_review.csv](seven_plane_review.csv)，[paired_transport.csv](paired_transport.csv) |
| 下一步决策 | [4T_AUTHORIZATION_ZH.md](4T_AUTHORIZATION_ZH.md)，条件式授权草案与成本；当前 NO-GO |
| 图 | Figure_A_background：时间/模态；Figure_B_spatial：浅层噪声；Figure_C_transport：带符号的浅层输运；均有 PNG/PDF |
| 验收与溯源 | `review_delivery_checks.json`、`artifact_manifest.json`、`control_acceptance_and_cost.json`、`amendment_20261103_B01.json` |

本轮 1–2 mm 累计差异是**持续、超过收支残差、但极小且未做完整配对收敛验证**的早期信号。不能表述为零输运，也不能表述为已证实的有效盐输运或推迟溶解度达到。两条基准在 2T 均删失；`t_eff` 未定义。

## 复现

在本 Paper 2 独立仓库根目录执行：

```text
python scripts/test_wp6_paired_stopping.py
python scripts/analyze_wp6_review.py
python scripts/build_wp6_review_products.py
python scripts/verify_wp6_review.py
```

这些命令不启动求解器。原场留在远端，并压缩归档于忽略目录 `tmp/wp6_review_20261103`；绘图缓存亦在该目录。重新取得背景缓存可用 `python scripts/wp6_review_background.py baseline` / `half`，仅在 frodo 后处理已经完成的命名运行。`wp6_review_campaign.py collect dt_half` 和 `archive dt_half` 只收集既有对照，不应重复 stage/launch。

运行前的 B01 仅登记唯一 dt/2 对照与原阈值不变，没有追认新对称性门槛。准备过程原样保存于 `provenance/prepare_control.py`，是历史执行记录，不应重跑。
