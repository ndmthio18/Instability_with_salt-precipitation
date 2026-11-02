# command/3.md 持续执行交接

状态：**sal M0 running / M1 running, authorized parallel pair**，不是已完成的 Pilot 结果。
用户已授权 command/3 全工作包，并授权使用 cxhpc 或 sal。已迁移至 sal；用户最新明确授权 M0/M1 各一个核并行，覆盖原串行限制。

## 当前执行身份

- 本仓库：`D:/3OutputFile/1Research/Article/22Surrogate_model/paper/2`，分支 `codex/paper2-main`。
- sal：`/tmp/zhixin_paper2_pilot0_20261102_r01`。
- **M0 Slurm job 19160，M1 Slurm job 19161，均已 RUNNING**；原 afterok 依赖已解除，没有新提交或重复模拟。
- 两条腿同时运行，各 1 rank、1 个物理核心、BLAS/OMP/MKL=1、墙钟上限 12 h；M1 启动检查 CPU 亲和性与 M0 不重叠。
- frodo worker 2453288 已在启动 M1 前取消；`frodo_queue_cancelled.json` 留存证据。禁止恢复该队列。
- frodo 原 case42929 未修改，继续运行。sal 不必等待它。
- `paper-2-pilot-0` heartbeat 已改为每 30 min 检查 sal；无变化静默，完成、失败或决策时通知。
- 正式比较用新 sal M0/M1；原 frodo M0 仅作迁移桥接，避免主机变化混入微弱物理信号。

## 运行来源与限制

原 ELF `d78ce5122cf69e4dee5a4d12bcc726779959f304b343246009cb45bb37e894cc`、21 项源码、原输入及 nx28/1e-6 profile 未变。
原 frodo 动态加载器及依赖库随包部署；204 项文件通过 SHA 校验。`sal_migration/` 保存部署与输入记录。
MPI 为单 rank singleton；仅求解器子进程移除 Slurm/PMIx 自动发现变量，仍在 Slurm 分配和 CPU 约束内。
首次加载探针发现 MPI 路径/启动方式不兼容，已补齐 OpenMPI/PMIx 组件、共享数据并通过加载，随后才启动 M0。

sal home 配额已满；独立计算目录位于 `/tmp` tmpfs。输入在本地有备份，完成后必须回传小文件并归档 VTU，核验 SHA，保留远端证据。
Slurm 站点 RealMemory=1M 无法申请正常内存；使用单核分配、进程地址空间 8 GiB 限额，预检实际可用内存约 215 GiB。
cxhpc 当前主机密钥检查失败，未绕过，未在该机启动任务。

## 必须使用的操作

```powershell
python scripts/pilot0_sal_campaign.py collect
```

读取 `sal_migration/latest_status.json`，不要再 prepare/submit。失败先查 Slurm、wrapper 与 solver 日志，不盲目重提。
本轮 checkpoint/restart 未验证，VTK 不能直接作为重启状态。

1. M0 完成后 worker 比较原 frodo M0 的完整预算/表面历史，生成 `bridge_validation.json`。
   这是主机桥接测试，容差不是物理有效输运阈值。并行不跳过最终门槛；失败时不得发布配对物理结论，并评估停止尚在运行的 M1。
2. M1 第一接受步后自动检查实际初始 VTK：600×200 坐标、与新 M0 的压力/Sw、nx28/1e-6 公式及三个 BOX 区域初始库存。
   `initial_field_check.json` 失败即停止。本次启动已通过实际网格、压力、Sw、扰动公式及初始库存6项检查；快照见 `sal_migration/M1_initial_field_check_snapshot.json`。
3. 完成后 collect 回传小文件；执行下列归档将完成腿的原始场备份到本地忽略目录：

```powershell
python scripts/pilot0_sal_campaign.py archive
```

4. 正式分析只用新 sal 配对：

```powershell
python scripts/analyze_matched_pair.py --m0 analysis/pilot0_integrated_20261102/sal_migration/runs/M0 --m1 analysis/pilot0_integrated_20261102/sal_migration/runs/M1 --out analysis/pilot0_integrated_20261102/sal_baseline_pair
```

未完成或 wrapper 失败必须拒绝；不得用 diagnostic-only 替代验收。接受时刻不同可显式用 `--common-accepted-times`，只比较精确共有物理时刻的原始累计量，保留未覆盖时段，不能按行号相减。

最初 5 步桥接输出相同，仅为局部检查，完整结果 pending。早期每步约 20–35 s，600 步的 M0 粗估 3–6 h；M1 启动后实测估计，不能套用原 frodo 成本保证完成。

并行变更来源：`sal_migration/parallel_activation.json`、`parallel_resource_preflight.json`、`parallel_gate_checks.json`。检查时62个空闲逻辑CPU、约217 GiB可用内存。原 M0 正在执行的包装器未改动，M1 使用留存旧版本后的新包装器；原 ELF、输入、profile 不变。`enable_parallel_remote.py` 是已执行记录，禁止重复执行。正式分析接口强制核验 M0 bridge 和 M1 initial check 及其输出哈希。

## 继续工作顺序

- 完成 Task 2 后自动继续 Tasks 3–4：完整 C0/C1、Delta C、初始修正 D_I、adv/diff、表面均值/极值、扰动幅值时序和实际初始库存核验。
- 每接受步 sigma/Fourier/mode 已由同 ELF 的 `*_sigmatop.csv` 输出，不能把 M0 Constant 分支的内置投影当作 nx28 投影；共同基底比较需从场重投影。
- 现有 adv/diff 是总守恒面通量分解。`Delta adv` 不等于独立的扰动协方差对流通量，不能用 VTK velocity*x 替代守恒面通量。若需后者，须明确诊断缺口及版本/桥接成本。
- 从新 M1 的 sigma 重新提取 t_detect，使用 Paper 1 冻结提取器（只读复制并记录 hash）；不能套用 50394.334 s 历史检测值。
- 基础预算通过且差异值得分辨时，Task 5A 顺序运行同物理窗口的 M0/M1 dt/2；先重估墙钟。Task 5B 先做嵌套分级网格和实际切面预检；Task 5C 有条件开展。
- 若扰动弱、无可分辨额外输运，不机械加密；按 Task 6 的实际 max(x)、增长率和成本决定有界双方延窗，绝不增幅值。超预算先给成本方案。
- Task 7 在基准评价后做 case38852 WP1 设计；仅满足 command/3 的条件时运行最小配对。不能用密度冲突旧场，也不能据两案例宣称单参数因果。
- 完成能合理推进的 Tasks 8–9，补齐六项交付与 Figures A–D，逐项写 completed/pending/failed/censored/inconclusive。
- 每个实质里程碑更新本项目 `PROGRESS.md`，提交推送；完成本轮并汇报后暂停 heartbeat `paper-2-pilot-0`，不要归档用户聊天。

## 来源问题与已完成检查

`preflight_attempt01_failure.json` 保留首次登记拒绝。旧准备计划引用的是本地 M0 输入/plan 文本哈希；远端 LF 输入 hash 与实际 M0 plan 一致。
本轮使用 `M0_actual_params.input` 的精确输入 hash，并比较 plan 的 JSON 语义一致性，实际远端 plan hash另录在 `preflight.json`。
注册 M1 plan 中 `frozen_preparation_M0_*` 和继承的 `M0_plan_sha256` 都指旧本地准备引用，实际 M0 字节以新 `runs/M0/` 和 `preflight.json` 为准。
这是来源换行处理，不是模拟物理/数值设置更改。旧文件没有改写。

`execution_checks.json`：13 项执行门槛、正式输入接口、实际 M0 预算和自身相减检查通过。
未完成的 U_dt/U_z/U_cut/U_x 均为 unknown，不填 0；不注册任何 1%/3%/5% t_eff 阈值。
