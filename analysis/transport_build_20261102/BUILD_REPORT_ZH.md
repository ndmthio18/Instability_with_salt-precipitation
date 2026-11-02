# Paper 2 隔离源码、输运诊断与编译报告

> 2026-11-02 命名核对：LSA `t_crit` 是全局临界模态首次 crossing；`t_seed(0.10)` 是增长筛选（含 fallback）选中的 observable 模态自身 crossing；`t_detect` 是 DuMuX 检测时间。`t_eff` 独立且尚未定义。原始字段按[统一映射](../onset_semantics_20261102/CONTRACT_ZH.md)阅读，不能把 `t_seed` 解释为达到指定幅值的时间。

2026-11-02。**完成实现、静态检查、配置和编译；没有运行求解器或任何模拟。** 编译成功不代表数值守恒、网格收敛或有效输运已验证。

## 交付位置

- 源码与说明：`../../model/transport_v1/README.md`。
- 主诊断实现：`../../model/transport_v1/source/appl/simpleRichardsnc/transportdiagnostics.hh`。
- 共享离散内核：同目录 `model.hh`；时间循环、VTK 修正、停止逻辑：`main.cc`。
- 来源、补丁：`upstream_manifest.json`、`upstream/`、`provenance/`、`transport_v1.patch`，均位于 `model/transport_v1`。
- 编译证据：`attempt_01/`、`attempt_02/`，各自保留命令、stdout/stderr 日志、配置、源文件哈希、编译器/依赖版本和 ELF 信息。
- 最终静态检查：`static_checks.json`；最终本地 ELF 身份：`local_artifact.json`。

## 1. 来源边界

20 个 upstream 源码/配置文件从 frodo 非饱和模块只读取得。远端父仓库观察到 HEAD 为 `854d4fdd330392d899f5c239b296829e90bf7a39`，与 Paper 2 原冻结 handoff 提交不是同一身份；本次明确新增独立源码快照，**未修改 `paper1_sources.lock.json`，也未将这个 HEAD 自动升级为数据基线**。

原始文件的 SHA-256 在获取后保存，构建结束后再次读取共享源文件，20 个文件全部一致。Paper 2 开发仅发生于本地 `model/transport_v1/source/` 与远端 `/temp/zhixin/paper2_transport_20261102/source/`。饱和冻结模型、非饱和共享源码与历史结果未修改。换行归一化后的开发增量用补丁记录，upstream 快照保留原始字节。

## 2. 实现的物理量与控制体

使用实际模型的 `phi Sw x` 和 `q x - D_eff grad x`，不加入变化的摩尔密度。盐库存/通量按单位水平面积归一化，分别为 m 和 m/s，不能直接标记为 kg 或 mol。

三组上部控制体对应目标 z/H=0.9、0.75、0.5。实际切面选最近顶点层中点并写入文件；另外输出不含底部 Dirichlet 节点的自由控制域。SCV 储存与 SCVF 通量直接调用求解残差的共享计算函数，保持原上风 mobility、浓度、面平均扩散系数和面面积。用实际接受 dt 的步末通量进行后向 Euler 累计，输出 advective/diffusive/total 和每步/累计残差。

全域库存额外包含底部控制体，底部交换由 Dirichlet 节点控制体平衡反推。此全域收支与自由域平衡等价，不构成另一项独立边界通量测量。所有解释在输出元数据中保留。

修正可视化速度缺少 krw 的问题，并在 t=0 写 VTK 前计算；它仍是单元中心重构量，不能替代守恒面通量。原 sigmatop 数据列保留，新的表面文件增加面积均值、局部与全域 max(x)、全域 min(x)。

默认表面或内部最大浓度达到 0.1 后终止，保存最后低于阈值和首次达到阈值的步区间；跨阈接受步本身保留并标记，不宣称精确事件根或盐结晶动力学。stop record 区分 horizon、表面/内部阈值、异常/不完整和启动失败。对 SIGKILL 等不可捕获中断，只能保留最后写入的检查点。

## 3. 编译结果

| 项目 | 结果 |
|---|---|
| 主机/目录 | frodo，`/temp/zhixin/paper2_transport_20261102/` |
| 编译器 | GCC 13.3.0，`/usr/bin/c++` |
| CMake | 3.28.3 |
| 配置 | Release，C++17 能力，`-O3 -DNDEBUG`；`-Wall -Wextra -Wpedantic` |
| 核心依赖 | DuMuX 3.9，DUNE 2.10；现有可选 dune-foamgrid 为 2.9.1，精确 Git IDs 见 result.json |
| 依赖发现 | 对现有 DUNE/DuMuX build-cmake 前缀使用 CMake package discovery，没有修改共享依赖 |
| 配置 / 编译退出码 | 两次尝试均为 0 / 0 |
| 第一次 | 主诊断功能编译成功 |
| 第二次（最终） | 增加诊断对象构造前的 startup failure 记录后重新编译成功 |
| 最终程序 | `build-gcc-release/appl/simpleRichardsnc/paper2_transport` |
| ELF 大小 | 1,615,760 bytes |
| 最终 SHA-256 | `d78ce5122cf69e4dee5a4d12bcc726779959f304b343246009cb45bb37e894cc` |

已将最终 Linux ELF 拉回 `model/transport_v1/artifacts/paper2_transport` 并重新校验哈希；本地和远端程序均未执行。第一次产物的本地副本与其构建记录也保留。二进制目录忽略于 Git，源码和审计证据保留。

完整命令在 `attempt_02/result.json/commands`，实际编译参数在 `compile_commands`。运行的构建命令只涉及 CMake 配置与 `--target paper2_transport --parallel 2`；未运行 `ctest`、应用程序或案例驱动。CMake 的编译器/依赖探测不属于模型模拟。

## 4. 静态检查与警告

`scripts/check_transport_static.py` 的 14 项检查全部通过，包括：原始源码身份、共享源码未变、当前本地源码与编译源码逐文件一致、残差和诊断共用内核、实际接受 dt 调用顺序、单 MPI rank 和覆盖保护、初态速度更新、阈值区间记录、底部 reaction 标签、Python 语法和补丁 `git apply --check`。

这些是结构与来源检查，不是运行测试。构建日志包含 **223 条 warning 记录**，不是零警告构建；主要位于既有 DUNE/DuMuX 依赖和上游代码，包含 VTK 模板的 signedness / maybe-uninitialized 提示。本次新增 `transportdiagnostics.hh` 没有编译警告。项目内的 unused-parameter 警告保留，未用关闭警告来掩盖。未修改共享依赖来消除这些警告，也未声称已验证相应运行路径。

## 5. 未执行、未证明的部分

- 没有运行任何模拟，连单步和零步求解器启动也没有。
- 没有生成真实的 transport_budget/surface_events 时序，没有测量收支残差量级。
- 没有独立 1D 程序、1D/均匀 2D 一致性或 MPI 汇总；当前显式限定单 rank 二维矩形 BOX。
- 没有网格/时间步加密，没有新 t_detect、t_eff 或 A/B 输运结论。
- 没有跨运行的累计差分或有效输运阈值冻结。准备脚本只支持未来的均匀/扰动 2D 案例资料准备，本轮未调用它建立真实案例。

下一阶段在获准开展数值验证后，先验证背景解、单位/符号、初始库存、残差面通量与库存闭合，再检查接受 dt 和阈值事件区间；通过后才运行 Pilot A/B。cxhpc 的 116 条历史文件及旧场来源冲突仍是 Phase 0 的独立待办，未被这次编译解决。
