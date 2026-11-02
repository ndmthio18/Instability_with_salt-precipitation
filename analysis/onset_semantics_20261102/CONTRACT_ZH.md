# Paper 2 onset 定义与字段约定 v1

本轮按 `command/2.md` 和用户随后指定的合作稿
`D:/3OutputFile/1Research/Article/22-1-surrogate_model/69e50893f5d9b55adfb56494/AWR_onstime/onstime.pdf`
执行。PDF SHA-256 为 `e57fc093419d388ca9bdf005ac5a5b08bb73ae5a402ead7b08a33c1a39bfb498`。
依据是该 PDF 的 Section 2.2–2.4，文件页码 8–14（印刷页码 7–13）；
不是桌面文件，也不是凭文件名推定的 `onstime(5).pdf`。
对应文字提取见 `paper1_sections_2_2_to_2_4.txt`，身份见 `source_provenance.json`。

## 四个符号、三类物理问题

| 类别 | 符号 | 严格含义 | 标准字段（秒） |
|---|---|---|---|
| LSA 状态稳定性 | t_crit | 搜索集 nx=1,…,30 上全局 critical mode n_c 的首次 neutral crossing，Rc(t)=R_E | `lsa_t_crit_s` |
| LSA 可检测性筛选 | t_seed(0.10) | growth-qualified observable mode 自身的 neutral crossing | `lsa_t_seed_theta010_s` |
| 数值检测 | t_detect | DuMuX 表面盐摩尔分数标准差最小值的协议相关检测时间 | `dumux_t_detect_s` |
| 输运有效性 | t_eff | 相对严格配对背景的输运差异达到经验证标准的物理时间；目前尚未定义 | `transport_t_eff_s`，始终空值 |

前两行均属于 LSA onset；数值检测是第二类问题，有效输运是独立第三类。
“global”限于 Paper 1 规定且已筛选的模态搜索，不宣称对无限连续波数证明全局最优。
两种 LSA 时间均用 t=hat(t) phi H/E 转换；标准化直接读取冻结秒值与无量纲值，不重新求解。
`t_seed` 不是达到某个有限扰动幅值的时刻，也不是 t_detect 的别名。
初始化在 t=0 施加所选模态的本征函数，不是在 t_seed 才施加扰动。

## Observable 模态的选择

对每个候选 nx 定义自身 crossing hat(t_c(nx))，并计算
`S(nx)=sigma(nx,2 hat(t_c(nx))) hat(t_c(nx))`。
候选为 `{n_c} ∪ {2,…,30,40,60,80,120}`，只考虑 nx≥n_c 且
`hat(t_c(nx)) < min(hat(t_precip),hat(t_90))` 的模式。

1. 选满足 S≥theta=0.10 的最小 nx，标记 `preferred`。
2. 若没有，选上述时间可行模式中 S 最大者；若 S≥S_floor=0.03 则接受，标记 `fallback_best`。
3. 否则不赋予 observable 模态或 t_seed。默认启用 fallback；不得将缺失补为 t_crit。

`lsa_observable_selection_theta010`、`lsa_observable_status_theta010`、sharpness、阈值、floor 和协议 hash 必须随时间保留。
主 first-event 分类始终基于 t_crit，不改用 t_seed。
0.10 和 0.03 是 sharpness 筛选值，不是盐浓度或扰动幅值阈值。

## 检测类别与实际初始化模态

`initialized_mode_nx` 是实际运行模态，`lsa_t_initialized_mode_crossing_s` 是该模式自身的 LSA crossing。
它只有在该运行确实初始化 observable 模态时才对应该运行的 t_seed 参考。
模态 30/60 的运行不能因为 case_id 相同就被赋予模态 28 的 crossing 作为自身参考。

Paper 1 的检测使用表面 std 最小值及其 0.1% plateau，plateau 不超过记录长度的 1/4 且随后上升才归为 resolved。
`dumux_t_detect_s` 只填 `resolved_onset` 的点值。`dumux_t_detect_estimate_s` 保留旧提取器的候选数值，
`dumux_t_detect_sample_s` 是离散采样最小值；low/high 是最小值 plateau 的时间区间，不是概率置信区间。
weak-interval、immediate-growth、no-onset/censored 状态原样保留，不提升为已分辨的点事件。
无扰动试算的持续可见波动判据（例如 sigma>1e-8 连续 5 步）另属 visibility，既不是 Paper 1 的 t_detect，也不是 t_eff。

## 原字段映射：必须先识别来源

| 来源及原字段 | 标准解释 |
|---|---|
| handoff `observations.csv`: `critical_nx`, `t_onset_s`, `t_hat_onset` | `lsa_global_critical_mode_nx`, `lsa_t_crit_s`, `lsa_t_hat_crit`，连同分类和 censoring |
| handoff `seeded_modes.csv`, sharpness_min=.10: `nx_seeded`, `t_onset_seeded_s` | `lsa_observable_mode_nx_theta010`, `lsa_t_seed_theta010_s`，仅 status=ok |
| handoff `seeded_modes.csv`: `critical_nx` | 全局 n_c；与 observations 核对。no_sharp_mode 行可为空，不能反过来清空全局记录 |
| handoff `dumux_onsets.csv`: `critical_nx` | 从案例级全局记录连接的 n_c；不同于下行同名字段 |
| 历史 `case_manifest.json`: `lsa.requested_nx=28`, `lsa.critical_nx=28` | 指定模式 28 的限制搜索结果，不是全局 n_c |
| 同 manifest: `lsa.t_onset_lsa_s`, `eigenprofile_t_hat` | 指定模式自身 crossing，后者也是生成 eigenprofile 的参考时间 |
| handoff `seeded_nx`, `t_lsa_mode_s` | `initialized_mode_nx`, `lsa_t_initialized_mode_crossing_s`；保留 `t_lsa_mode_source` 来源 |
| handoff `onset_s`, `onset_sample_s` | 拟合检测时间、采样最小值；结合 `onset_class` 使用 |
| manifest `dumux.onset.onset_refined_s`, `onset_s` | 拟合检测时间、采样最小值；**manifest 的 onset_s 与 handoff 的 onset_s 含义不同** |
| Phase 0 `lsa_t_seed_010_s`；dense `t_seed_s` | 案例 canonical t_seed(0.10)，不一定是该运行实际模式的 crossing |
| WP1 `t_lsa_mode_s`, `seeded_nx` | 指定模态 manifest 的秒值与 requested_nx；原 manifest SHA 逐行保留 |
| 旧背景计划 `legacy_detect_s` / `old_detect_s`，结果 `legacy_detection_s` | 历史扰动运行的 t_detect 参考，不是无扰动运行检测事件；原舍入值保留 |
| `first_sustained_visible_s` / visibility criterion | 描述性波动可见时间，独立保留，不映射为任一种 onset |
| `t_eff_defined=false` | `transport_t_eff_status=undefined_not_a_label`；不得由任何 LSA/检测时间填充 |

## 表、报告和接口的使用方式

`scripts/onset_semantics.py` 是统一读取接口，只读经过 SHA 校验的冻结 handoff。
`FrozenOnsets.case(case_id)` 读案例级定义，`run(raw)` / `run_by_key(run_key)` 附加运行级定义。
不导入求解器，不选择新模态，不改变原始来源值。

历史 Phase 0 表、remote snapshot、manifest、注册运行计划和命令文档是不可改写的来源证据。
新增的 `views/run_audit.csv`（436 行）、`wp1_supplement_audit.csv`（84 行）、
`pilot_screening.csv`（209 行）、`highk_screening.csv`（129 行）保留原表每一列，并追加标准字段；
新分析应读这些视图或使用上述接口。`table_mapping.json` 记录原表到新视图的对应关系。
`views/pilot_candidates.json` 保留原嵌套对象并增加 `normalized_onset`。
`views/case_onsets.csv` 是这些运行/候选涉及的案例级表，不能作为全库人口统计。

`dense_history/run_assessment.csv` 原列保留并直接增加标准字段；重新生成的 17 页诊断 PDF 和代表时间线
分别标出 t_crit、canonical t_seed 和实际初始化模态的 crossing，粉色区间明确属于 DuMuX 检测。
库存/场/预算/阈值表没有新增 onset 值，原物理时间、积分和阈值事件保持不变。
各报告的统一命名注记适用于其历史 onset 用语；旧日志和冻结快照按此映射阅读，不追改历史证据。
研究方案中有关 t_eff 的阈值和时间线仍是待验证的候选设计，不是已定义的标签。

运行计划的旧字段保留兼容；未来生成/收集背景记录的脚本附加 `legacy_onset_reference`，明确它指向历史扰动运行。
成对输运接口显式输出空的 `transport_t_eff_s`。本轮未执行任何 launch/collect 远程动作。

## 复现与范围

在 Paper 2 根目录依次执行：

```powershell
python scripts/prepare_paper1_inputs.py --verify-only
python scripts/audit_onset_semantics.py --pdf 'D:/3OutputFile/1Research/Article/22-1-surrogate_model/69e50893f5d9b55adfb56494/AWR_onstime/onstime.pdf'
python scripts/analyze_dense_history.py
python scripts/write_dense_history_report.py
python scripts/check_onset_semantics.py
```

本轮是术语、来源和存量数值映射审核，不是重新验证全部 LSA 特征值或 DuMuX 物理模型。
