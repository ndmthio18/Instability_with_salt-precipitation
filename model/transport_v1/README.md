# Paper 2 isolated transport model v1

本目录是非饱和 Richards–salt 模型的隔离开发版本。实现、静态检查及编译后，已获准在 frodo 完成 8 项小规模背景/事件配置，见 [验证报告](../../analysis/transport_validation_20261102/VALIDATION_REPORT_ZH.md)。归一化守恒量采用 `phi*Sw*x`，不宣称真实盐质量守恒或 Pilot 可行性已验证。完整 Pilot A/B 尚未启动。

## 源码与来源

- `upstream/`：2026-11-02 从 frodo 非饱和模块获取的 20 个源码/配置文件的原始字节快照，只读基线。
- `upstream_manifest.json`：远端父仓库 HEAD、读取路径、获取时间及每个文件 SHA-256。远端工作文件的逐文件哈希才是本隔离版本的精确来源；父 HEAD 不能单独证明这些工作文件已提交。
- `provenance/`：Paper 1 已冻结的非饱和 `CHANGES.md` 和原补丁，说明 darus-4610 来源及此前改动。
- `source/`：可编译的 Paper 2 独立模块。DUNE 内部模块名称沿用 kiemle2023b，但源目录、构建目录和可执行目标 `paper2_transport` 均与 Paper 1 分离。
- `transport_v1.patch`：从本目录 upstream 快照到 source 的完整增量，包含新增诊断头文件；静态检查脚本重建并执行 `git apply --check`。
- `artifacts/paper2_transport`：已校验的 Linux ELF 本地副本，忽略于 Git；不能作为 Windows 程序运行。

不能把本目录的当前 upstream 快照称为原始 DaRUS 未修改模型；它已经包含 Paper 1 的非饱和闭合参数、Eigenprofile 和 onset instrumentation。原饱和模型、Paper 1 源文件和原结果目录未修改。

## 实现内容

`source/appl/simpleRichardsnc/model.hh` 将储存和面通量提取为同一个共享内核：求解残差与诊断调用相同的 `normalizedStorage` 和 `transportFlux`。保持原来的浓度/压力梯度、面密度/扩散系数算术平均、上风相对渗透率、上风浓度和面面积。分量按 advective、diffusive、total 输出。公式保持一致，但未运行数值回归，不能宣称与旧二进制逐位相同。

`transportdiagnostics.hh` 实现：

1. `z/H=0.9, 0.75, 0.5` 三个目标深度对应的上部控制体，外加 `free_domain`。按最近的相邻顶点层中点选择切面，输出实际位置；不是直接用几何插值截面。
2. 按 BOX 的 SCV 体积累加 `phi Sw x` 和 `phi Sw`；按 SCVF 的 inside/outside 归属对共享面通量定向。每个局部内部面只遍历一次，完整遵循残差对内外两个节点一加一减的约定。
3. 使用 Newton 返回后的实际接受 dt 和步末面通量进行后向 Euler 累计；没有按稀疏 VTK 时间序列做梯形积分。
4. 输出盐/水的累计与逐步收支残差，初始库存差留给配对后处理显式处理。每条运行记录自身初始库存，尚未计算跨运行的 Delta C 或 t_eff。
5. 顶部面积加权均值、局部最大值、全域最大/最小值，以及与原检测器口径一致的顶部节点均值/标准差。
6. 保存表面和全域浓度首次达阈的接受步区间，默认达到 0.1 后停止。跨阈那一步仍保存，标记已观测越界；不是把越界时刻当作精确根，也没有沉淀动力学。
7. 每步刷新 `stop_record.json`；正常终点、表面阈值、内部阈值、异常/不完整和启动失败分开记录。SIGKILL/断电只能留下最后成功写入的 `in_progress_or_interrupted`，不能自动推断具体失败原因。

VTK 的 `velocity` 计算现在包含单元中心重构的 `krw/mu`，且在初始场输出之前计算。字段名称保持 `velocity`；**它仍不是残差的守恒面通量，不能替代 transport_budget。**

## 单位、底部边界和能力限制

库存和通量均除以水平宽度 W，相当于固定单位出平面厚度下的单位水平面积量。归一化盐库存单位 m，通量 m/s；没有乘变化的摩尔密度，也没有标成 kg/m2 或 mol/m2。使用固定参考浓度换算物理单位需另行声明。

底部 Dirichlet 节点不是自由演化方程。`free_domain` 排除这些节点的全部控制体，最下控制面位于第一对顶点层中点。全域库存包含底部控制体；底部输入由其储存变化与指向自由域的面通量反推。全域 residual 与自由域收支等价，**不是另一项独立边界通量验证**。元数据已明确这一限制，初始步底部 reaction 标记未定义。

当前支持 **单 MPI rank、二维矩形 YaspGrid、BOX、顶部与侧面零盐通量、底部 Dirichlet、指数 krw**。不支持 MPI ghost 汇总、任意网格、任意盐边界或真正的一维程序。横向均匀二维参考已通过短背景测试；真正的 1D 模型与 1D/2D 等价性仍未实现/证明。

不检查/冻结 t_eff 的持续时间阈值，不生成有效输运标签；这些属于 Pilot 后续工作。已做一次横向加密和时间步减半的短背景敏感性检查，未做垂向收敛。预算门槛使用预设绝对下限；没有证明所有近零输运信号均大于数值误差。

## 输出接口

| 文件 | 内容 |
|---|---|
| `transport_budget.csv` | 每接受步×4 个区域；目标/实际深度、盐/水库存、盐通量三分量、向下水通量、其他边界通量、累计量、逐步及累计残差、全域库存、底部 reaction、越界标志 |
| `surface_events.csv` | 每步顶部面积均值、max(x)、全域 max/min、原节点统计、面积与达阈标志 |
| `transport_metadata.json` | 单位、控制体定义、实际深度、积分方式、底部 reaction 的解释、阈值与停止开关 |
| `stop_record.json` | 最后接受时刻、步数、状态及表面/全域阈值包围区间；逐步原子替换 |
| `run_manifest.json` | 通用准备脚本记录输入/profile/二进制/构建及无量纲量；本次小验证由专用脚本记录配置、源码/ELF/输入/输出哈希和退出原因，E=0 等合成测试不解释为生产无量纲案例 |

`cumulative_net_water_out_m` 包含控制面向下水流和顶部/侧面水流；不能解释为只有向下水流。盐顶/侧边界当前为零，因此 `cumulative_cut_salt_m` 为控制面累计盐通量。`inventory_initial_salt_m` 是本次运行的初始库存，不自动等同于以后配对使用的共同 I_ref。

## 编译与静态检查

本地根目录命令：

```powershell
python scripts/build_transport_remote.py
python scripts/collect_transport_build.py
python scripts/check_transport_static.py
```

第一条只在 `/temp/zhixin/paper2_transport_20261102/` 同步本隔离源码并启动配置/编译，返回 PID；第二条在编译结束后取回日志；第三条检查来源哈希、源码/编译一致性、关键调用顺序、补丁可应用性和 Python 语法。没有执行 `ctest` 或 `paper2_transport`。构建尝试的日志分别保存，不覆盖上次证据。静态检查不是运行期守恒验证。

准备脚本 `scripts/prepare_transport_run.py` 只准备新目录，拒绝覆盖、检查 ELF 与构建结果哈希、要求显式流体常数、复制 profile 并记录实际输入。它没有任何执行求解器的代码路径。小规模验证使用独立的 `scripts/run_transport_validation.py`，冻结 8 项配置后顺序执行；`collect_transport_validation.py` 与 `analyze_transport_validation.py` 只收集和分析。专用验证允许 E=0、初始达阈、非法时间等明确登记的测试，不放宽通用案例准备器的检查。

完整编译命令、依赖 Git 版本、警告及二进制 SHA-256 见 `../../analysis/transport_build_20261102/BUILD_REPORT_ZH.md` 和其机器可读结果。
