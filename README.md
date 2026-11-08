# Paper 2 — effective salt redistribution

`research/p2-h-crystal-boundary` 将默认激活阈值设为3.5 μm，保持零初始固体、
真实质量初始化及epsilon_birth=0.43，并在隔离应用中修复顶部侧面误分类。
HPCC四网格均成功激活；40/80/160层完成2500 s，320层及dt/2在激活后快速
进入ALE状态保护边界。联合预算通过，完整四网格和定量物理准入仍为NO-GO。
见[本轮报告及七图](analysis/p2_h_crystal_boundary_20261108/REPORT_ZH.md)和
[复算记录](analysis/p2_h_crystal_boundary_20261108/REPRODUCE.md)。

`research/p2-h-adaptive-donor` 增加了守恒、限深的多层孔隙液初始化，保留
`top_only` 旧模式；ALE、BDF和土壤物理方程保持不变。新方法解决80/160层的
局部单BOX供液限制，但320层的高度容差把侧面误判为顶面并施加蒸发，
额外面积0.03617%解释了激活前的水预算失败；本轮未改动冻结边界。完整四网格验证
仍为NO-GO。2 mm限深还会使40/80层10 um案例拒绝激活，不能隐去这些结果。
见[本轮中文报告](analysis/p2_h_adaptive_donor_20261108/REPORT_ZH.md)、
[复算说明](analysis/p2_h_adaptive_donor_20261108/REPRODUCE.md)和
[初始化假设](model/p2_h_ale_coupling/ADAPTIVE_FORMULATION.md)。

此前 `task/9.md` 的 **P2-H 条件 DuMuX–ALE 水盐耦合**已完成独立分支验证：
C++隐式库存积分复现冻结M3，从零固体守恒激活，激活后真实双向交换且不再重复
执行P2-G析盐扣除。HPCC最终8组实验均到2500 s，其中5组形成ALE、3组按供液
筛查保持预盐壳；联合水盐误差低于1.11e-13。盐壳网格与时间步门通过，但土壤
加密后激活失败，完整空间收敛及水力物理有效性仍NO-GO。见
[中文报告和图件](analysis/p2_h_ale_coupling_20261107/REPORT_ZH.md)、
[模型与复算入口](model/p2_h_ale_coupling/README.md)、
[方程及边界](docs/P2_H_FORMULATION.md)。

2026-11-07 的 **P2-G 最简表面析盐**：保留 Richards–salt，
通过表面溶解度投影守恒地累积外部固体，并独立对照 Prat 增高需求；不使用人工
晶种或盐层内部反应求解。HPCC 七条 1000 s 基准已完成，联合水盐预算通过，
Prat 需盐量约为同期析盐量的 73%；80/160 层析盐量差 2.081%，未过 2% 门槛。
代码链条可运行，定量物理高度和二维失稳结论尚未建立。见
[模型入口](model/p2_g_surface_precipitation/README.md)、
[公式与边界](docs/P2_G_FORMULATION.md)、
[验收报告](analysis/p2_g_surface_precipitation_20261107/REPORT_ZH.md)。
以下保留历史研究结果，不代表本轮启动了其他算例。

`command/16.md` 的 WP9 已完成：R96-M1 在 sal 单线程、R96-M0 在用户单独授权后于 hpcc
四线程完成，均为 9600/9600 接受步并通过守恒、时间轴、逐帧哈希和独立 VTU 再积分。
R48/R96 的预注册 G1–G5 判定为 G1/G2 通过、G3/G4/G5 失败，最终 **NO-GO**；4T 的
能量加权波数相差 46.54%，顶面背景参照响应比为 0.285，1 mm 初始校正输运量相差
65.83%。小扰动窗内 LSA–DuMuX 增长率 11 个准入点的相对差中位数 5.44%、最大 9.78%。
R96 同时减半 `dx` 和 `dt`，`nz=200` 未变，盒内又缺少全局临界 `nx=5` 模态，故不能
声称纯空间收敛或唯一因果归因。`t_eff`、结晶延迟及训练标签均未定义；R192、多波长、
全域和 8T 未启动。42 项数值测试与 1020 项交付检查全过。详见
[交付索引](analysis/wp9_r96_convergence_20261104/README_ZH.md)、
[三分辨率比较](analysis/wp9_r96_convergence_20261104/R24_R48_R96_COMPARISON_ZH.md)和
[最终判定](analysis/wp9_r96_convergence_20261104/DECISION_ZH.md)。

`command/15.md` 已完成：唯一一次连续修正严格背景从 t=0 走到 4T，4800/4800 接受步正常结束，48 项强制验收与 657 项交付检查全部通过，4801 帧逐帧哈希归档。与历史接受前缀 4414 个共有接受时刻的差异全在舍入量级（max|Δx| = 9.02e-16，max|Δp| = 3.06e-10 Pa），旧失败步 4414 以 3 次迭代被接受，Newton 最大迭代由 18 降至 5。有了这条背景后，R24/R48 的 M0/M1 零模态响应可在 0→4T 上参照：3T 之前两种横向分辨率之比为 1.003–1.015，4T 升至 382，且本应等于背景的无扰动 M0 自身在 4T 偏离 1.24e-06（R24）对 2.90e-09（R48），因此 4T 幅度尚未收敛、不可引用。全程未触及局部溶解度阈值（最大 x = 0.056969，占 0.1 的 57.0%）；未定义 t_eff，未启动 R96、追加二维、多波长或 8T。详见 [交付索引](analysis/wp8_full4T_20261103/README_ZH.md)、[执行与验收](analysis/wp8_full4T_20261103/EXECUTION_AND_ACCEPTANCE_ZH.md)、[历史对照](analysis/wp8_full4T_20261103/HISTORICAL_COMPARISON_ZH.md)和[科学审查](analysis/wp8_full4T_20261103/SCIENTIFIC_REVIEW_ZH.md)。

`command/3.md` 基准工作包已完成：sal M0/M1 均正常结束并归档；三深度额外输运尚未分辨，Phase 2 暂不具备启动条件。六项交付与四类图见 [结果报告](analysis/pilot0_integrated_20261102/PILOT_RESULTS_ZH.md)和[下一阶段决策](analysis/pilot0_integrated_20261102/NEXT_STAGE_DECISION_ZH.md)。本轮未追加模拟，自动轮询已结束；重现步骤见 [执行交接](analysis/pilot0_integrated_20261102/EXECUTION_HANDOFF_ZH.md)。

已完成 `command/2.md` 的 [onset 来源审计](analysis/onset_semantics_20261102/AUDIT_ZH.md)与[字段约定](analysis/onset_semantics_20261102/CONTRACT_ZH.md)。新分析使用 `analysis/onset_semantics_20261102/views/` 中保留原列并附加标准字段的表；原始来源表和注册计划保持冻结。Case 36877：全局 nx=5、t_crit=31409.16 s；observable nx=28、t_seed(0.10)=45474.35 s；历史 t_detect=50394.33 s。t_eff 尚未定义，本轮未重算或启动模拟。

本目录已独立为一个 Git 仓库。研究方案见
[Option_A_Research_Plan.md](Option_A_Research_Plan.md)，当前任务是 Phase 0
审计与 Pilot-0 准备，已完成小规模背景及诊断验证，完整 Pilot A/B 尚未启动。

2026-11-02 已完成可访问范围的 Phase 0 审计：见
[审计报告](analysis/phase0_20261102/PHASE0_AUDIT_ZH.md)、
[输运诊断约定](analysis/phase0_20261102/TRANSPORT_DIAGNOSTIC_CONTRACT.md) 和
[Pilot-0 候选案例](analysis/phase0_20261102/PILOT0_CANDIDATES_ZH.md)。
运行表重现 436 条 / 253 个案例 / 293 条 resolved primary；frodo 的 320 条及
后续 84 条已检查，cxhpc 的 116 条文件仍未核实。新增发现包括归一化守恒量、
VTK 速度缺少 krw 因子，以及 case 38852 旧场与当前输入的密度公式冲突。
候选为 case 36877 / 38852；尚未证明有效输运，完整 A/B 配对尚未启动。

已建立 [隔离输运模型](model/transport_v1/README.md)，实现 `phi Sw x` 预算、
共享离散面通量、表面/全域浓度阈值及停止记录。静态检查与 frodo 编译通过，
详见 [编译报告](analysis/transport_build_20261102/BUILD_REPORT_ZH.md)。
随后在 frodo 完成 8 项小规模背景/事件配置，结果见
[验证报告](analysis/transport_validation_20261102/VALIDATION_REPORT_ZH.md)。
均匀二维背景、归一化盐水收支及事件记录通过预设门槛；近零通量使用了绝对误差下限，
不能宣称所有微弱输运信号均已分辨。垂向网格收敛、真正 1D/2D 一致性和 Pilot 可行性仍待验证。

进一步完成 [14 条历史密集时序分析](analysis/dense_history_20261102/REPORT_ZH.md)：
覆盖 10 个物理案例、1,414 帧。输入—字段一致性和时间轴修复通过；首次观测越界前，
三个上部控制体库存均随保存时间增加。case 36877 模态 60 虽有约 3.9 倍 std 放大，
末帧 sigma 仅为表面均值的 0.0066%。现有数据支持早期富集/扰动诊断，尚不足以识别
相对于无扰动背景的额外再分配；本次只做后处理，没有启动新模拟。项目另行进行的
600×200 无外加扰动试算已正常结束，见 `analysis/unperturbed_36877_20261102/RESULT_ZH.md`。
高渗透率 case 42929 的延长无扰动试算正在进行；注册输入和阶段快照位于
`analysis/highk_unperturbed_42929_20261102/`，快照不代表最终结果。

已完成 [case 36877 严格配对准备审计](analysis/matched_pair_36877_20261102/AUDIT_ZH.md)及
[最小运行计划](analysis/matched_pair_36877_20261102/MINIMAL_RUN_PLAN_ZH.md)：现有 M0 可复用，
同一 Paper 2 ELF 的 nx=28、幅值 1e-6 M1 尚缺；待运行输入和原始 eigenprofile 已冻结。
配对诊断经预算重算、自身相减和代数/拒绝条件检查；未定义 t_eff，未启动新增模拟。

## 仓库边界

| 内容 | 位置与管理方式 |
|---|---|
| Paper 1、公共基线代码和冻结数据归档 | 父级 `22Surrogate_model` 仓库，既有远程与分支保持不变 |
| Paper 2 研究方案、代码、文稿、分析记录 | 本仓库，分支 `codex/paper2-main` |
| Paper 1 数据版本和文件校验值 | `paper1_sources.lock.json` |
| Paper 2 使用的本地外部数据副本 | `data/external/paper1/`，已忽略，不重复提交 |
| Paper 2 进展 | `PROGRESS.md` |
| 正式合作稿和 11 月 1 日报告源文件 | 独立 Overleaf 仓库，定位见 `PAPER1_REUSE.md` |

本仓库不是父仓库的分支或子模块。父仓库忽略整个 `paper/2/`，因此 Paper 2
的修改不会进入 Paper 1 的提交。按用户指定，本仓库的独立 GitHub 远程为
`git@github.com:ndmthio18/surrogate_model.git`，远程名 `origin`，工作分支为
`codex/paper2-main`。推送 Paper 1 不会备份 Paper 2；本仓库需单独推送。

## 准备 Paper 1 输入数据

在本目录执行（仅 Python 标准库）：

```powershell
python scripts/prepare_paper1_inputs.py --verify-only
python scripts/prepare_paper1_inputs.py
git status
```

默认读取本目录上两级的 Paper 1 仓库。若本仓库移到别处：

```powershell
python scripts/prepare_paper1_inputs.py --paper1-root 'D:/path/to/22Surrogate_model'
```

脚本核验固定提交存在、选定源文件及 LFS 归档的 SHA-256，再从归档中按白名单
读取 handoff 文件，写到以基线提交号区分的数据目录。已有同名文件必须通过
校验才会复用；不覆盖不匹配文件。数据目录还保存独立 provenance 清单。
源码校验允许 Windows CRLF 与 LF 的纯换行差异；归档和 handoff 数据仍进行
逐字节 SHA-256 校验。

数据准备不代表 Phase 0 已完成：原始二维场、离散盐通量、停止原因和 onset
之后的覆盖长度仍需逐运行审计。方案中的 436 次运行、253 个物理案例和 293
条配对观测是待追溯的报告统计，不能由旧 handoff 直接替代。

## 日常工作

```powershell
git -C 'D:/3OutputFile/1Research/Article/22Surrogate_model/paper/2' status
```

新代码与分析写在本仓库。不要在父仓库现有模型上直接试改 Paper 2 功能；
需要改变数值模型时，先建立可追溯的隔离源码副本并记录上游版本及补丁。
不要把新的输运结果写回 Paper 1 的结果目录。
