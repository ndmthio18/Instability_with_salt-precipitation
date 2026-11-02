# Onset 命名和来源一致性审计

结论：Case 36877 的历史 `critical_nx=28` 来自指定模态的 crossing 接口，不能作为全局临界模式。
Paper 1 原始数据、Paper 2 已冻结输入和模拟结果均保留；新增字段和报告采用 [统一定义](CONTRACT_ZH.md)。
本轮没有重算 LSA、重新提取 DuMuX onset 或启动模拟。

| Case 36877 项目 | 冻结值 |
|---|---:|
| 全局 critical mode n_c | 5 |
| t_crit | 31409.157792677746 s |
| theta=0.10 observable mode | 28 |
| t_seed(0.10) | 45474.35459078279 s |
| selection / sharpness | preferred / 0.10104692701116125 |
| 历史 nx28 manifest 拟合 t_detect | 50394.33363312911 s |
| 同运行 handoff 拟合 t_detect | 50394.33363314537 s |
| 采样最小值时间 | 50476.53359577 s |
| 检测 plateau 区间 | [44110.12395306, 56615.57146552] s |
| t_eff | 未定义，不赋值 |

两个来源的拟合 t_detect 差约 1.63e-8 s，远小于采样尺度。各自精度原样保留，没有通过舍入覆盖差异。
上述检测值来自 nx28 主协议运行，不能移植为同一案例其他模态的检测时间。

## 生成路径与证据强度

审阅代码来自锁定的 Paper 1 提交 `158395f3d25a29680dea5bcd0d3a24e52a28b4cd` 的 Git blob，
不是当前工作副本的隐式导入。逐文件 hash 和行号见 `source_provenance.json`、`SOURCE_EXCERPTS.md`。

1. `dominant_mode_cli.py` 的 nx_map 将 observable 模式写为 `nx`，全局模式另存 `critical_nx`。
2. `run_design_campaign_parallel.py:174–176` 读 nx_map 的 `nx`，构造 `nx_by_case`。
3. `pipeline.py:445` 将 case 的模式传作 `nx`；`prepare_and_run_case` 将它传入 LSA/分辨率审核。
4. `solvers.py:325` 保存 `requested_nx`；`run_lsa` 的扫描、二分及最终 `critical_rayleigh` 都传递 `nx=nx`。
5. `lsa_unsaturated.py:578`：仅 nx=None 才生成 1..nx_max；给定 nx 时搜索数组只含指定模式。
6. 返回的 `critical_nx` 是这个搜索数组的最小值位置。`pipeline.py:336` 把整个 lsa 对象写入 manifest，未重命名。

该 manifest 同时记录 requested_nx=28、critical_nx=28、eigenprofile_t_hat=0.0017177719703013337，
以及自身 crossing 秒值，与冻结 observable 行一致。因此 nx_max=30 不能证明它做过全局扫描。
`least_stable_critical_nx` 同样处在受限搜索上下文，不能当作全局字段。
原 manifest 没有完整的历史启动 argv 和源码提交身份：已确认字段语义与存储路径，
但没有声称追溯出当时具体使用的 nx_map 文件，或独立证明历史二进制身份。

## 存量映射检查

- 436 条旧运行和 84 条 WP1 补充记录分别标准化，不把两批相加当作独立物理案例数。
- 两组 Pilot 候选、209 行候选筛选、129 行高渗透率筛选及 14 条密集时序共用同一字段约定。
- 全部冻结 theta=0.10 行共 10,985 条：10,718 条 accepted，其中 preferred 10,586 条、fallback 132 条；267 条 no_sharp_mode 保持空值。
  这是冻结 seeded_modes 文件总体，不是 Paper 1 回归筛选后的 9,354 个案例，二者统计单位不能混用。
- 代码用 `t_nx > t_stop` 排除越界，PDF Eq.(10) 使用严格 `<`。全部 10,718 条 accepted 存量记录均严格早于时间限制，
  没有等于或超过边界的案例；本轮记录这一代码/文字边界差异，不修改 Paper 1 或重选模式。
- sharpness、fallback 和时间限制检查使用存储值；候选最优化/特征值未重算，不能据此声称重新证明每个模式最优。
- 密集时序原 `t_seed_s` 是案例 canonical 值。标准字段将其与模态 30/60 自身 crossing 分开，避免错误计算延迟或比值。
- weak-interval 等检测类别继续保留。其候选拟合时间仍可追溯，但不填入 resolved-only 的 `dumux_t_detect_s`。

旧证据文件通过哈希保护，新的命名视图保留每个原始字段。最终检查详见 `verification.json`。
研究结论不变：现有数据尚未建立相对配对背景的有效盐输运，t_eff 仍不能定义或用于训练。
