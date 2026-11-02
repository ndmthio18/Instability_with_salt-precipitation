# command/3.md 持续执行交接

状态：**in progress / queued_resource_gate**。本目录是执行中的研究工作包，不能当作已完成的 Pilot 结果。
用户已通过 `command/3.md` 授权本轮 M1、小规模验证、条件性双方延窗及提交推送；原 command/1 计划中的“待批准”已被本授权取代。

## 当前执行身份

- 本仓库：`D:/3OutputFile/1Research/Article/22Surrogate_model/paper/2`，分支 `codex/paper2-main`。
- frodo：`/temp/zhixin/paper2_transport_20261102/pilot0_integrated_20261102_r01`。
- 已提交 worker PID：2453288；**这是排队包装器，不是求解器 PID**。
- 本聊天 heartbeat：`paper-2-pilot-0`，每 30 分钟继续检查和推进。无状态变化时静默；完成、失败或需要决策时汇报。
- 只允许该 worker 执行已登记的唯一 M1，没有自动重试、网格批次或延窗分支。
- 参数/profile 与旧 `planned_M1` 字节完全相同，ELF 为 `d78ce5122cf69e4dee5a4d12bcc726779959f304b343246009cb45bb37e894cc`，21 项源码 hash 已核验。

## 为什么尚未开始求解

资源检查时 frodo 有原 case 42929 在算，PID 2416984，约 27% 完成，负载约 25/24。
其他研究任务也在占用计算资源；不终止或修改它们。M1 等待全部其他 `paper2_transport` 进程结束，
并要求 1 min 负载≤20、5 min 负载≤22、可用内存≥4 GiB、磁盘空余≥10 GiB。
每 60 s 检查；排队最多 12 h，超时成为 `pending_resource_queue_expired`，不无限重排。
资源门槛通过后，M1 自身最多执行 12 h。单 rank、BLAS/OMP/MKL=1、nice=10。
若停止/失败，保留输入、日志、最后接受状态、返回码和所有输出哈希。

原 M0 2492.69 s 墙钟只提供成本参考；旧 nx28 不同 ELF 的 19122 s 也不能保证新 M1 成本。
当前负载下原 case 42929 可能还需数小时。不得为了赶进度并行再开大型 Paper 2 求解。

## 必须使用的操作

```powershell
python scripts/pilot0_campaign.py collect
```

读取 `latest_status.json`。**不要再次执行 prepare 或 submit**。远端独立目录、launch claim 和 flock 防重复。
`worker_failure.json` 优先于任何旧 worker_state；若有失败，先查日志与实际进程，不能盲目重提。
本轮尚未验证 checkpoint/restart，不能把 VTK 当作可直接重启状态。

## 运行中及完成后的门槛

1. 第一接受步后，worker 对实际新 M1 的 t=0 VTK 进行独立核查：600×200 坐标、与 M0 的压力/Sw、原 profile 的 nx28/1e-6 公式、三个 BOX 区域初始库存。
   输出 `initial_field_check.json`。失败即停止该新进程；不使用历史 nx28 初始库存替代。
2. 完成后 `collect` 下载全部小型实际输出及哈希；VTK 留在 frodo，后处理应远端进行，再回传指标和代表图/压缩场。
3. 正式分析命令（M1 尚未完成时必须拒绝）：

```powershell
python scripts/analyze_matched_pair.py --m0 analysis/pilot0_integrated_20261102/runs/M0 --m1 analysis/pilot0_integrated_20261102/runs/M1 --out analysis/pilot0_integrated_20261102/baseline_pair
```

不得用 `--diagnostic-only` 代替真实配对验收。`validate_pair_inputs` 已通过，但 `validate_pair` 的完整输出验收尚 pending。
如果实际接受步不同，现接口拒绝错配；可显式加 `--common-accepted-times` 使用已通过9项检查的精确共有接受时刻交集，保留原始累计量，明确未覆盖时段。
只有必要时才开发经过测试的守恒一致时间对齐，不按行或步号直接相减。

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
