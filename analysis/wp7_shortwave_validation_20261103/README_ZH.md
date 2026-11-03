# WP7 交付索引

- `R48_PLAN_AND_ADMISSION_ZH.md`：冻结协议、用户授权覆盖、R48 执行与 R96 成本审计。
- `R24_R48_COMPARISON_ZH.md`：共同物理时间的增长、模态、深度、七面输运和阈值比较。
- `BACKGROUND_1D_VALIDATION_ZH.md`：一单元背景与二维 M0 的部分验证及失败门槛。
- `LOCAL_SCVF_DIAGNOSTIC_ZH.md`：局部守恒面通量能力、隔离编译和一步桥接；不冒充晚期机制证据。
- `FINAL_REVIEW_ZH.md`：最终科学判断及 R96/多波长/8T 的 GO/NO-GO。
- `r48_comparison.json`、`decision.json`、`r96_cost_estimate.json`：机器可读结论与成本。
- `R24_R48_checkpoints.csv`、`final_mode_energy_R24_R48.csv`、`paired_transport_R24.csv`、`paired_transport_R48.csv`、`transport_R24_R48_final.csv`：可复算表格。
- `Figure_A_R24_R48_growth.*` 至 `Figure_D_R24_R48_transport.*`：PNG/PDF 证据图，均已视觉核查。
- `sal_collect_r48.json`、`sal_collect_spatial.json`、`sal_collect_r96_cost.json`：归档和逐场哈希回执。
- `runs/`：输入、计划、在线诊断、独立复算和小型证据；完整 VTU 归档位于忽略目录 `tmp/wp7_shortwave_validation_20261103/`。
- `failed_attempts/`：早期迁移包绝对路径导致的两次零步失败，保留而未覆盖。

`r48_preflight_summary.json` 是用户资源授权前的历史快照，后续执行状态以 `sal_collect_r48.json` 和 `r48_comparison.json` 为准。

本轮实际完成 R48 M0/M1 4T 配对和条件性 R96 十步成本测试；没有启动 R96 长时、多波长或 8T。
