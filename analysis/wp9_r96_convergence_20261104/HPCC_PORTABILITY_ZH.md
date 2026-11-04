# hpcc（探索1000）可移植性验证

**范围：只验证基础设施，不产生任何科学结果，也不授权任何生产算例。**
本轮在 sal 上进行的唯一一条 R96-M1 未受影响，仍按 `command/16` 的冻结协议运行到底。

## 1 为什么要做这件事

sal 的计算节点此刻负载 100/128，R96 的单步耗时从 1.3 s 漂到 2.7 s。要把后续算例
（dt 对齐的 R96、R96-M0 或 nx=192，均需单独授权）放到别的机器上，前提是**机器本身
不构成协议差异**——这正是父仓库对 frodo/cxhpc 的既有要求。本节用一次十步基准把这件事
测出来，而不是假设。

## 2 环境与可达性（只读巡检）

| 项 | 值 |
|---|---|
| 登录节点 | `ibcln01`，CentOS 7，内核 3.10.0-1160，宿主 glibc **2.17** |
| 工作区 | `/WORK/hulm_work/zhixin/`（Lustre，7.0 PB，已用 27%） |
| 空间配额 | 1.86 TB / 2.05 TB |
| **文件数配额** | **2 057 468 / 2 000 000，已超，宽限剩约 1 天 14 小时** |
| 余额 | 43 190 |
| 分区 | `cnmix` 15 个全空闲节点（840 核）+ 262 核混合空闲；`cnall` 多数 alloc/drain |
| 账户既有作业 | 1 个运行中（56 核），2 个 PENDING（依赖） |

**文件数配额已超是一个现实约束。** 一条完整 R96 要写 4801 个 VTU，只占超额量的 0.23%，
但在宽限期内增加文件不是好做法。真要在这里跑长算例，应先与账户持有人确认清理，
或改用更大的 `OutputInterval` 并及时打包。本次只增加了 361 个文件（bundle 350 + 11 帧）。

## 3 可执行文件为什么能跨过 CentOS 7

sal 的 `paper2_transport` 链接到 glibc 2.38，而 hpcc 的宿主 glibc 只有 2.17。
之所以仍能运行，是因为这个 bundle 是**自包含**的：它自带
`lib/ld-linux-x86-64.so.2` 与 `lib/libc.so.6`，而 `launch_solver.py` 直接用自带 loader
启动：

```
os.execve(loader,[loader,'--library-path',str(r/'lib'),str(r/'bin/paper2_transport')]+argv,env)
```

所以宿主的动态链接环境完全不参与，只要内核版本满足自带 glibc 的下限即可。
CentOS 7 的 3.10 内核满足。

## 4 十步基准：结果

Slurm 作业 `16286679`，分区 `cnmix`，1 节点 1 任务 1 核，`COMPLETED`、`ExitCode 0:0`，
节点 `ibc11b11n16`。作业脚本带有守卫：无 `SLURM_JOB_ID` 或落在登录/测试节点时拒绝计算。

| 项 | sal | hpcc |
|---|---:|---:|
| ELF SHA-256 | `1120767…df81` | **同一哈希** |
| 接受步 | 10 | 10 |
| 停止原因 | `requested_horizon` | `requested_horizon` |
| 末接受时刻 / s | 1136.8588647695694 | 同 |
| 求解器墙钟 / s | 25.647 | **15.570**（求解器自报 12.786） |

**逐值比较（四个诊断文件，2 948 个数值单元）：**

| 文件 | 行数 | 数值单元 | 不同的单元 | 最大相对差 | 逐字节相同 |
|---|---:|---:|---:|---:|:--:|
| `surface_events.csv` | 12 | 132 | **0** | 0.0 | 是 |
| `transport_budget.csv` | 89 | 2 640 | **0** | 0.0 | 是 |
| `weighted_online.csv` | 12 | 110 | **0** | 0.0 | 是 |
| `*_sigmatop.csv` | 12 | 66 | **0** | 0.0 | 是 |

**结论：同一可执行文件在 hpcc 上逐字节复现 sal 的结果，单核快约 1.65 倍（按墙钟）
或 2.0 倍（按求解器自报 CPU 时间）。因此主机不是协议差异**，与父仓库对
frodo/cxhpc 的既有判定一致，不应进入任何协议哈希。

## 5 这条结论能用来做什么、不能做什么

**能：** 后续任何单独授权的算例（dt 对齐的 R96、R96-M0、nx=192）可以放在 hpcc 上，
而不引入主机相关的数值差异；按 `cnmix` 的空闲核数，单核作业可立即开始。

**不能：**

- 不能用来改变 sal 上这条 R96。它已跑过一半，`command/16` 明确禁止重启与重试，
  且本轮的单次授权只覆盖那一条。
- 不构成对任何算例的授权。十步基准是基础设施检查，不是科学结果。
- 不能据此推断长算例的速度：十步包含启动开销，且两台机器的瞬时负载都会变。
- 文件数配额仍在宽限期，长算例前必须重新核对。

## 6 复算

```bash
cd paper/2/scripts
python wp9_hpcc_portability.py stage     # 需要 sal 上取下的 bundle（哈希已登记）
python wp9_hpcc_portability.py submit
python wp9_hpcc_portability.py status
python wp9_hpcc_portability.py collect
```

证据：`hpcc_portability/stage.json`、`job.json`、`environment.txt`、`sacct.txt`、
`solver.log`、`stop_record.json`、四个诊断 CSV 与 `portability.json`。
远端目录 `/WORK/hulm_work/zhixin/paper2_wp9_portability_20261104/`。

---

## 7 线程扩展：同一可执行文件本来就支持共享内存多核

**更正：** 本文件第 5 节原先只检查了 MPI（bundle 里没有 `mpirun`/`mpiexec`/`orterun`，
`launch_solver.py` 又写死 `OMPI_MCA_ess=singleton`），由此说"多核做不到"是不完整的。
该二进制**编译进了 TBB 后端的共享内存多线程**：符号表里有
`Multithreading::ExecutionBackends::TBB`、`DUMUX_NUM_THREADS` 和 `Assembly.Multithreading`，
并链接了自带的 `libtbb.so.12`；日志开头的 `Colored 19200 elements with 4 colors` 就是
着色装配被启用的证据。被并行化的区域包括 `FVAssembler::assembleJacobianAndResidual`、
`assembleResidual`、`CVFEGridVolumeVariables::update` 和 `CVFEGridFluxVariablesCache::update`。

以往所有算例都把它压在单线程上（worker 用 `sched_setaffinity` 钉一个 CPU，
`cpus-per-task=1`），所以这部分能力从未被使用。

### 实测（十步，`cnmix`，作业 16288117/18/19/20）

| 线程 | 节点 | 求解器 / s | s/步 | 加速 | 末步装配 / s | 末步线性求解 / s |
|---:|---|---:|---:|---:|---:|---:|
| 1 | ibc11b11n01 | 16.037 | 1.604 | 1.00 | 0.40 | 0.82 |
| 2 | ibc11b11n01 | 9.447 | 0.945 | **1.70** | 0.21 | 0.45 |
| 4 | ibc11b11n01 | 8.418 | 0.842 | **1.91** | 0.12 | 0.45 |
| 8 | ibc11b11n04 | 22.075 | 2.208 | 0.73 | 0.09 | 1.60 |

八线程那条落在另一个繁忙节点上，装配仍继续下降（0.09 s）而线性求解被挤到 1.60 s，
因此**不可比**，不能据此说八线程变慢。1/2/4 三条同节点。

### 结果逐位相同

`surface_events.csv`、`transport_budget.csv`、`weighted_online.csv` 在
串行参照与 1/2/4/8 线程五个版本之间**哈希完全一致**，且与 sal 的十步基准一致。
原因是着色装配让每个自由度收到贡献的顺序由颜色决定，与线程调度无关。

**所以多线程是零数值代价的加速：不改源码、不改输入、不改 ELF 哈希，4 线程约 1.9 倍。**

### 瓶颈转移

1→4 线程时装配从 0.40 s 降到 0.12 s（约 3.3 倍），而 ILU-BiCGSTAB 线性求解
只从 0.82 s 降到 0.45 s 就饱和；四线程时**线性求解占一步的约 79%**。
再往上加线程收益很小，下一个杠杆在线性求解器本身——
`main.cc:848` 有一行被注释掉的 `UMFPackIstlSolver`，对 38 994 个未知量的二维问题
直接法通常很快。那需要重新编译，会改变 ELF 哈希，属于需要单独登记和验证的协议分支。
