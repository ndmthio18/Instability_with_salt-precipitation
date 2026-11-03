# WP6 command/9 交付索引

**已完成唯一 M1 dt/2 计算、严格细步配对和七面时间敏感性分析。A 支持提出有界探索性 4T；B 定量输运及溶解度延迟仍 NO-GO。未启动 4T/8T。**

- [结果报告](RESULT_ZH.md)：运行验收、增长/检测、七面 U_dt、时钟差界、浅层全时程与背景噪声解释。
- [4T 决策和预算](DECISION_4T_ZH.md)：A/B 分开，具名 dt/2 探索分支的具体范围与限制。
- [方法与证据边界](METHOD_AND_SCOPE_ZH.md)：严格配对、初始库存修正、时间差处理及预登记工程判据。
- `amendment_20261103_C01.json`：新增 M1 前锁定的运行范围和工程判据；没有更改原门槛。
- `seven_plane_timestep_summary.csv`：七面终点及最新严格共同时刻的敏感性；`paired_timestep_sensitivity.csv`：每个实际比较时刻的 U_dt 和误差界。
- `refined_pair/`：2,401 个严格共同接受时刻，8 控制体共 19,208 行；各面通量、修正库存、残差及事件删失。
- `growth_endpoint_comparison.csv`、`growth_and_detection.json`、`vertical_profile_comparison.csv`：增长、检测、表面浓度和垂向结构。
- Figures A/B：增长与剖面、配对输运敏感性（PNG/PDF）。原始 VTU/NPZ 留在忽略目录，不提交 Git。
- `delivery_checks.json`、`test_results.json`、`artifact_manifest.json`：原基准与源码哈希保护、31 项测试、归档核验及来源清单。

最关键的限制：1–2 mm 末态信号同号、细化后分别下降约 6.8%/4.5%，但 **1 mm 在严格相同的约 1.02T 时刻差约 68%**。因此末态接近不能当作全时程输运已收敛。U_z、U_x、U_domain 仍 unknown，两种溶解度事件均删失，t_eff 未定义。
