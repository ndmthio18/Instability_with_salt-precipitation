# 输运诊断约定 v0.1 — 待实现与验证

> 2026-11-02 命名核对：LSA `t_crit` 是全局临界模态首次 crossing；`t_seed(0.10)` 是增长筛选（含 fallback）选中的 observable 模态自身 crossing；`t_detect` 是 DuMuX 检测时间。`t_eff` 独立且尚未定义。原始字段按[统一映射](../onset_semantics_20261102/CONTRACT_ZH.md)阅读，不能把 `t_seed` 解释为达到指定幅值的时间。

日期：2026-11-02。适用于 Paper 2 的非饱和、指数形式 krw 分支。此约定冻结概念和数据接口；t_eff 的主深度、阈值、持续性窗口与数值容差尚未通过 Pilot 校准。没有新模拟结果。

## 1. 实际守恒方程与量纲

依据本次 `source_snapshot/model.hh`，模型使用

\[
u_s=\phi S_w x_s,\qquad
\mathbf j_s=\mathbf q x_s-D_\mathrm{eff}\nabla x_s,
\qquad D_\mathrm{eff}=\phi^{1.5}S_w^{3.5}D.
\]

`u_s` 无量纲，`j_s` 为 m/s。它们是模型采用的归一化组分守恒量与通量。**不能把 x 当作 kg/m3，也不能把下述积分直接标成 mol 或 kg。** 主收支使用实际求解器的这套量。

设统一水平面积 A（2D 必须声明同一外延厚度 b，例如 b=1 m；1D 换算到相同 A），定义单位水平面积库存与通量：

\[
I_U(t)=A^{-1}\int_{V_U}\phi S_w x_s\,dV\quad[\mathrm m],\qquad
F(z_*,t)=A^{-1}\int_{A_*}\mathbf j_s\cdot\mathbf n_\downarrow\,dA\quad[\mathrm{m/s}].
\]

如果报告恒定参考摩尔密度 c_ref 乘以 I_U，可称“固定参考浓度换算的等效摩尔库存”，单位 mol/m2，并给出 c_ref 的来源；它不是变密度状态下真实摩尔库存的独立守恒验证。物理库存 `integral(phi Sw rho(x)/M_mix(x) x)` 可以另作诊断，但它不由当前离散守恒式保证。若要把真实质量严格守恒作为科学目标，须先重新审查模型近似，不能在后处理中偷偷换方程。

## 2. 配对模型与初始条件

- M0：同一 Richards–salt 方程、同一物性/闭合与竖向网格的 1D 参考。
- M0-2D：同一二维网格、横向均匀初值的背景对照。侧边界沿用实际零通量边界；不能擅自改为周期。
- M1：同一物理点的扰动 2D。建议候选主方案继承 theta=0.10 选定模态、对应垂向本征函数、t=0 初始化、绝对摩尔分数幅度 1e-6。这是待冻结方案，不混用 critical、seeded 或另行优化的最快增长模态。

必须保存扰动 profile 文件及哈希、归一化方式、实际 k 或 nx、初始库存差。波数名与波长须按当前边界和初始化公式映射，不能仅凭参数名称推断。初始场离散平均也必须检查，不能以连续余弦均值为零替代。

任何修改在 Paper 2 隔离版本完成；保留上游提交/归档、源码补丁、编译器、DUNE/DuMuX 版本、可执行文件完整 SHA-256。原运行目录只读。

## 3. 控制体、离散面与守恒检查

坐标 z=0 在底部，z=H 在顶部。初始探索 z*/H=0.9、0.75、0.5，即距表面 0.1H、0.25H、0.5H。BOX/CVFE 的控制体边界不一定与这些几何平面重合：优先选最近的闭合子控制体集合，保存实际 z*/H；否则明确几何切割方法及重构误差。跨案例/网格比较不得隐藏控制面偏移。

从 `computeFlux` 使用的同一状态、上风规则、有效扩散系数和面面积提取 advective、diffusive、total 三项。每个内部面只计一次，符号向下为正。库存使用实际 sub-control-volume 储存项和体积，不以可视化节点均值替代。

上部控制体的收支残差定义为

\[
R_U(t)=I_U(t)-I_U(0)+\int_0^t F(z_*,\tau)d\tau-B_U(t),
\]

其中 B_U 为除下控制面外的累计净盐输入（含边界和源项）。当前顶部与侧面无盐通量时 B_U=0。全域检查必须加入底部 Dirichlet 的交换/反作用通量；不能把 Dirichlet 方程替换后的残差行简单当作零通量。

配对差异满足

\[
\Delta I_U(t)-\Delta I_U(0)=-\Delta C(z_*,t)+\Delta B_U(t)+\Delta R_U(t),
\qquad \Delta C=\int_0^t(F_{2D}-F_{1D})d\tau.
\]

按求解器时间离散，在每个接受的时间步累计通量，输出 actual dt。使用后向 Euler 时应优先检验 dt_n F_n 对离散方程的闭合；稀疏 VTK 上梯形积分只作独立近似，不称离散精确收支。非线性容差、初值差和两套背景解差异分别报告。

0.5% 仅是原方案的探索性总量残差目标；是否足够必须相对于 Delta C 的信号判断。建议后续检验误差包络明显低于最小待解释信号，例如小于信号的 10%，但该数值尚未验证，不是本次通过标准。

## 4. 输运指标与解释

\[
\Delta F=F_{2D}-F_{1D},\quad
\mathcal T=\Delta C/I_{\rm ref},\quad
I_{\rm ref}=I_U^{1D}(0).
\]

I_ref 随案例和深度计算，但定义一致，且 1D/2D 共用同一分母。必须同时给出 Delta C 的绝对值、I_ref 和归一化结果，标记过小分母。不得取正部积分替代有符号积分；若另算只计向下贡献的指标，须另命名。

Delta F>0 只表示相对参考更向下，F_2D 本身仍可能向上。报告二者及 adv/diff 分解。模型连续平流项可以分解为 mean(q_z x)=mean(q_z)mean(x)+mean(q'_z x')，但实际离散上风面分解与节点协方差不完全相同，应单列差额。Delta F 不等于扰动协方差。

VTK velocity 在本次审阅源码中缺少 krw，不能直接输入上述诊断。若后续重构 q，需单独验证其与残差面通量的差异。

## 5. 事件时钟与适用窗口

- t_crit：中性稳定性首次穿越；t_seed：指定增长资格规则选中模态的 crossing；不是实际施加扰动时间。
- t_detect：既定检测器的 onset，保留原类别和时间区间。
- t_sat_LSA：线性基态的阈值参考，不等于非线性 1D 或 2D 的真实阈值时钟。
- t_sat_1D、t_sat_2D：各自表面局部最大 x 首次达到 0.1；保存包围步区间和插值约定。
- t_peak：有效窗口内第一个经过确认的表面平均浓度局部峰值。若最大值位于截断终点，标 `no_interior_peak`，不能宣称进入 flux-growth 阶段。

主配对窗口为 `min(horizon, t_sat_1D, t_sat_2D)`；数值失败更早则记质量截断，不当作正常删失。为模型有效性同时监测全域 max(x)；若内部先过阈，也停止物理有效解释，保留表面事件的独立状态。

t_eff 先从绝对时间累计曲线定义，不强制大于 t_seed。探索 1%、3%、5% 和三个深度，并保留阈值首次跨越区间、持续时间与后续回落。持续窗口尚待 Pilot，不能事后为获得分区而挑阈值。没有事件时不得填无穷大或赋成最终时刻用于普通回归。

主事件状态：`effective_before_sat`、`sat_before_effective`、`horizon_censored`、`numerical_failure_or_incomplete`、`unresolved_event_order`。另一参考模型先达到阈值导致配对诊断终止时，使用 `comparison_window_closed_by_reference`；不能据此断言二维自身一定没有有效输运。时间区间重叠时保留无法确定先后的状态。

## 6. 最小新增输出接口

| 文件 | 每条记录的最低内容 |
|---|---|
| run_manifest.json | case/source IDs、7 个无量纲量、实际物理参数、所有闭合分支、源码/二进制/输入/profile 哈希、网格、边界、dt 和停止政策 |
| transport_budget.csv | 接受步号、time_s、dt_s、控制体/实际深度、I_U、I_ref、F_adv/F_diff/F_total、累计 C、其余边界输入、R_U、总域水和盐收支 |
| surface_events.csv | time_s、面积加权表面 mean(x)、max(x)、全域 max(x)、原 t_detect 所需节点统计、阈值包围区间 |
| stop_record.json | 请求终点、最后接受步、终止原因、局部阈值标志、失败与重试来源；无论成功失败都保留 |
| fields.pvd + VTU | 保存 x、Sw、p、rho 与经过定义核验的速度；时序频率与轻量预算输出独立 |

先验证代码输出与方程一致，再校准网格/时间步和指标；现阶段不生成或提交模拟任务。
