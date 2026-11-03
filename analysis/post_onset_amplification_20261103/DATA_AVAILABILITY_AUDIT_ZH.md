# WP4 数据可用性与审计边界

2026-11-03。执行 `command/4.md`；本轮仅分析已有数据，没有运行 DuMuX、重算 LSA、修改远端原始文件或 Paper 1。

## 统计单位与完整性

主集 436 条冻结运行、253 个物理案例，344 resolved、78 weak-interval、6 immediate-growth、8 right-censored。12 条未达到请求终点（9 resolved、3 censored），incomplete 是正交的完整性标记，不能作为第五种互斥 onset 类别相加。主集 359 primary 标记不是 359 个独立物理案例。

WP1 为独立辅助来源：84 条、48 个案例，70 resolved、14 weak；其中 7 个案例与主集重叠。不能把 520 条当作独立土壤样本。`run_metrics.csv` 用 source_group 明确分组，`case_summary.csv` 每来源每 case 一行。病例代表行优先 is_primary，之后按 run_key 排序，选择不依赖增长大小；另列每案例最大 resolved G，明确这是有利方向的筛选统计。

已从 frodo 只读读取 404 份 sigmatop（主集 320、206 个案例；WP1 84）。每份 SHA256 与既有 Phase0 快照逐字节匹配。404 次重新检测均与冻结检测类别一致。主集剩余 116 条在 cxhpc：主机密钥验证失败，未绕过验证；其冻结端点放大和类别仍可用，G_max、拟合、相对幅度不能由端点反推。访问证据见 `cxhpc_access.json`。这 116 条中有 110 条 resolved，故全量 344 的端点统计与原始时序 234 的统计必须区别。

主集已知 VTU 保存数：292 条两帧、14 条 101 帧、12 条一帧、2 条零帧、116 条未核验。14 条密集运行覆盖 10 个案例。复用此前 1,414 帧的已验哈希产品和逐帧时间校正，不重读原始大场。新选 10 条稀疏记录，读取 18 个不同初/末场；其中 30858、29422 只有初始场。未选择的可访问两帧记录标记 available_not_selected，不能写成数据缺失。

320 条主集原始时序中 308 条达到请求终点；WP1 84/84 达到。未到终点的原因按原审计 stop_interpretation、returncode、record_complete 原样保留，无法区分故障/人工终止的记录不猜测。43 条主集表面均值曾达到 0.1 mol/mol，意味着局部阈值至迟此时已经达到；均值未达到不证明所有节点均有效。历史按模式 crossing 倍数设定的终点不等于物理饱和停止。初末场不能排除两帧之间更早的阈值穿越。

## 来源与定义

固定 Paper1 提交 `158395f3d25a`，`prepare_paper1_inputs.py --verify-only` 通过：12 个来源文件、归档和 10 个 handoff 文件。遵循 `analysis/onset_semantics_20261102/CONTRACT_ZH.md`，分别保存全局 t_crit、canonical t_seed(0.10)、实际初始化模态 crossing、DuMuX t_detect。t_eff 保持 undefined。

主集 106 份及 WP1 84 份 sigmatop 还保存 fourier_top 和 mode_proj，已与 sigma 对照。其他时序不伪造模态信息。物理参数来自各运行已哈希输入；协议、网格、振幅、初始模态、LSA 来源、阈值、终点和保存帧数都保留在 `run_metrics.csv`。

新增场审计隔离 case35733、19674、38852 的具体 hf1 运行：输入预测密度与场密度相差约 1.80–2.01 kg/m³，而正常场误差约 5e-12 kg/m³。记录诊断值仍可描述，但不混入经来源核验的物理证据。没有把当前可执行文件视为历史二进制证明。

## 复现与文件

执行顺序：`collect_post_onset.py` → `analyze_post_onset.py` → `collect_post_onset_fields.py` → `analyze_post_onset_fields.py` → `summarize_post_onset.py` → `plot_post_onset.py` → `report_post_onset.py` → `test_post_onset.py` → `verify_post_onset.py`，脚本均在本项目 scripts/。采集有哈希缓存；远端只读。原始时序和选中场的本地缓存位于忽略的 `tmp/post_onset_amplification_20261103/`；不把原始大场提交 Git。

`series_source_manifest.json` 记录 404 份 CSV 的远端路径/原始哈希/本地缓存；`sparse_field_evidence.json` 记录18帧哈希及已有LSA元数据；`dense_reuse_hashes.json` 校验70个既有产品；`source_manifest.json` 冻结所有本轮输入/脚本；`artifact_manifest.json` 校验交付文件。缺失的 cxhpc 原始时序、部分终场、历史 ELF 全证据和有限增长率输出保持 missing/unknown，不填0。
