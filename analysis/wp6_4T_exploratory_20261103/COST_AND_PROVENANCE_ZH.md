# 执行成本与证据保全

仅执行 command/10 授权的两条新物理模拟：`WP6-M0-explore4T`、`WP6-M1-explore4T`。frodo 管理器 PID 2889702，实际亲和性分别为物理核 3、9，其 SMT 同胞为 15、21；启动前检查同胞占用，未与已有忙任务共用物理核。启动资源记录：可用内存约 21.9 GB、磁盘约 233.7 GB，1 分钟负载 6.74。

## 实测成本

| 项目 | M0 | M1 |
|---|---:|---:|
| 接受步数 | 4800 | 4800 |
| 求解墙钟 / s | 740.985 | 711.296 |
| 求解 CPU / s | 709.169 | 681.390 |
| 求解峰值 RSS / KiB | 45492 | 45664 |
| 后处理及压缩墙钟 / s | 689.193 | 647.714 |
| 整体受限 worker 墙钟 / s | **1446.419** | **1374.408** |
| 整体 worker 硬限 / s | 3600 | 3600 |
| 求解输出 / byte | 4316155814 | 4453378121 |
| 原始场压缩归档 / byte | 702270232 | 980508795 |

整体 worker 包含启动核验、求解、逐帧后处理、压缩及结束清单；两腿并行约 24.1 分钟全部结束，均正常返回、未超时。合计占核墙钟 **0.783563 核时**，低于 2 核时预算；求解 CPU 合计 **0.386266 核时**。不得以 RSS 一列代替预留内存需求，也不把求解墙钟与整体 worker 墙钟混为一谈。

随后对已经完成的每腿各执行一次独立空间统计（不调用求解器），仍固定在各自释放的物理核，每次设置 900 s 超时。未采集这两次空间统计的独立 CPU 计费值，保留为 unknown；即按两次各 900 s 的保守占核上界计入，总远端占核墙钟上界为 **1.283563 核时**，仍在 2 核时内。本地分析、绘图与文件传输另计，不属于物理求解。

较长的后处理/压缩受磁盘 I/O 影响，不能仅按约 12 分钟求解时长预算更大域。下一阶段不自动申请、提交或扩展任务，建议和成本使用原则见 `NEXT_STAGE_ZH.md`。

## 来源、归档和哈希

- 源基线提交：`c9241f2c1f769b7e7fa263b1488d1fa16bc2db21`，分支 `codex/paper2-main`。
- D01 在启动前登记；每腿输入与对应旧细步输入只有名称、结束时间变化。配置、profile、D01 和 ELF 哈希保存在各腿 `plan.json`、`started.json`。
- 已验证 ELF SHA256：`112076759fdff0942000307f729c2bd3c7e96963229b43c1b67c1f12bbf9df81`。本轮不修改或重编译物理源码。
- 启动器、worker、管理器和原后处理远端字节与本地脚本逐一核对，见 `remote_execution_source_hashes.json`。
- 远端根：`/temp/zhixin/paper2_wp6_4T_exploratory_20261103_r01`。原始输入、日志、VTU、压缩文件和后处理保留，不清理远端证据。
- 本地完整原始场归档：`tmp/wp6_4T_exploratory_20261103/WP6-M0-explore4T_fields.tar.gz` 及 `WP6-M1-explore4T_fields.tar.gz`，共约 1.683 GB；两腿各 4801 个 VTU，共 **9602 帧**，校验压缩文件 SHA256 及每个成员的 SHA256。校验收据为 `explore4T_archive.json`。
- 日志和各个结果文件 SHA256 在 `finished.json`，逐 VTU 清单在 `post/field_hashes.json`。附加空间统计重新读取全部场并再次核验，`spatial/` 和 `spatial_M0/` 的 audit/transfer 保存源脚本、源场与小产品哈希。
- 原始 VTU/PVD、NPZ 和压缩包保持 gitignored；Git 保存输入、报告、图、结果表、日志与必要清单。`artifact_manifest.json` 是本轮交付文件的 SHA256 索引，`delivery_checks.json` 汇总核验。

## 复核入口（不会启动物理模拟）

在 Paper 2 根目录运行：

```text
python scripts/analyze_wp6_4t.py
python scripts/plot_wp6_4t.py
python scripts/verify_wp6_4t.py
```

绘图需要本地忽略缓存中的 `spatial/selected_fields.npz` 与 `spatial_M0/selected_fields.npz`。完整交付校验还需要归档和每腿忽略缓存；仅克隆 Git 不会获得这些大型证据。必要时从保留的远端执行 collect/archive，而不是再次 stage/launch。新建配对结果须使用 `analyze_wp6_review.py` 的 `--m0 --m1 --out` 接口和一个新输出目录；默认路径属于旧结果，不能覆盖旧历史。
