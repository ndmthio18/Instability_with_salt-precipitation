# P2-H 双移动边界守恒扩展报告

## Target

本工作检验三个问题：

1. 给定盐壳底部向土壤的扩散盐通量后，大 \(Da\) 下顶部生长是否对该通量弱敏感；
2. 底部晶体损失应表示为孔隙率增加还是几何底面上移；
3. 能否在不重复计算溶解的前提下，实现顶部生长、底部退缩和盐壳整体上移。

交付对象是独立解析验证、守恒双前沿诊断原型、A--G 数值测试、收敛表和可复现图件。生产 DuMuX 土壤耦合不是本阶段的验收对象。

## Status

**COHERENT AFTER REFRAMING / EXTRA ASSUMPTION**

解析边值问题和规定速率的双移动边界守恒模型是闭合的；底部界面溶解的动力学、水源和产物分配尚无独立物理闭合。因此，本实现是 **PASS diagnostic, NOT physically predictive**。这与任务的 acceptance rule 一致，但不支持“已建立预测性土壤--盐壳双前沿模型”的结论。

## Invariant Object

贯穿推导和实现的主对象是每单位投影面积的守恒库存：

\[
E_i=\epsilon_iV_i,\qquad
Q_i=\rho_l C_iE_i,\qquad
M_{cr}=\rho_{cr}\sum_i(V_i-E_i).
\]

任何几何运动、体反应或界面反应都必须在溶解盐、晶体盐和水账本中出现一次且仅出现一次。顶部高度、底部高度和厚度是几何量，不能代替晶体库存。

## Assumptions

- 一维、单位投影面积、液相饱和的盐壳；
- \(\rho_l\) 和 \(\rho_{cr}\) 在诊断算例中为常数；
- \(a_v,\epsilon_s,D_{s,eff},k_r\) 在解析切片中均匀且为常数；
- 解析问题采用准稳态 \(Pe\ll Da\) 约化；
- \(B_{soil}>0\) 表示从盐壳向下进入土壤的扩散盐通量；
- 顶部外法向向上，\(z\) 向上；
- M1 中底部界面反应和体反应互斥。原型通过底部排除区避免重叠；
- 底部退缩速率 \(v_b\) 和溶解产物向盐壳孔隙分配比例 \(\alpha\) 是规定输入，不是由土壤状态预测；
- 在常密度饱和孔隙中加入 \(\alpha D_{base}\) 的溶解盐时，以等质量水排出作为诊断闭合。该闭合不是活度/密度模型。

## Notation and sign convention

\[
z_b(t)=\text{bottom},\quad z_t(t)=\text{top},\quad h=z_t-z_b,
\]

\[
v_b=\dot z_b,\qquad v_t=\dot z_t,\qquad \dot h=v_t-v_b.
\]

底部向上退缩时 \(v_b>0\)，顶部向上生长时 \(v_t>0\)。\(G_{top}\ge0\) 是顶部界面固化率，\(D_{base}\ge0\) 是底部界面晶体溶解率，\(R_{bulk}>0\) 表示体内沉淀、\(R_{bulk}<0\) 表示体内溶解。

## Derivation Strategy

1. 从给定常系数反应--扩散方程独立求解下端 Neumann、上端 Robin 边值问题；
2. 将其与 Prat 下端 Dirichlet 问题逐边界比较；
3. 以总晶体库存为不变量，先写质量恒等式，再引入几何退缩；
4. 由 \(\xi\) 坐标映射推导仿射网格速度和有限体积相对通量；
5. 用规定前沿速率建立诊断原型，先验证几何守恒和库存守恒，再讨论土壤耦合所缺的物理闭合。

## Derivation Map

1. 解析顶部速度依赖于顶部过饱和度；顶部过饱和度由下端通量和上端盐平衡共同决定。
2. \(v_t(u)/v_t(0)\) 是同一边值问题内的精确比值，不依赖 \(Da^*\)。
3. Prat 结果只有在改变下端边界条件或 \(Da\to\infty\) 时恢复。
4. 双前沿固体账本来自 \(M_{cr}=\rho_{cr}(h-\sum E_i)\)，不是独立添加的经验源项。
5. \(D_{base}\) 只有在从体反应中排除同一反应区后才能加入。
6. 数值验证分成解析有限差分验证、M0 算子恢复、前沿运动守恒和网格/时间步敏感性。

## Main Derivation

### 1. 量纲审查

现有 ale.hh 使用

\[
R=\rho_l a_vk_r(C-C_{eq}).
\]

因为 \([\rho_l]=\mathrm{kg\,m^{-3}}\)、\([a_v]=\mathrm{m^{-1}}\)，要使 \([R]=\mathrm{kg\,m^{-3}\,s^{-1}}\)，必须有

\[
\boxed{[k_r]=\mathrm{m\,s^{-1}}}.
\]

若把 \(k_r\) 当作 s\(^{-1}\)，源项会错误地得到 kg m\(^{-4}\) s\(^{-1}\)。本分支已在 ale.hh 和 P2_H_FORMULATION.md 中明确该单位。对应的

\[
Da=\frac{(1-C_{eq})a_vk_rh^2}{\epsilon_sD_{s,eff}},
\]

\[
Da^*=\frac{\beta^2(1-\epsilon_s)^2(1-C_{eq})k_r}
{\epsilon_sD_{s,eff}a_v}
\]

均为无量纲。

### 2. 给定底部扩散通量的解析解

令

\[
\xi=\frac{z-z_b}{h},\qquad
\theta=C-C_{eq},\qquad \lambda=\sqrt{Da}.
\]

准稳态方程为

\[
\theta''-\lambda^2\theta=0,
\]

其中撇号表示对 \(\xi\) 求导。由于 \(B_{soil}>0\) 定义为向下扩散通量，底部边界为

\[
\rho_l\epsilon_sD_{s,eff}\frac{1}{h}\theta'(0)=B_{soil}.
\]

顶部盐平衡和界面动力学分别为

\[
\rho_l\epsilon_sD_{s,eff}\frac{1}{h}\theta'(1)
=C_{eq}j-\rho_{cr}(1-\epsilon_s)(1-C_{eq})v_t,
\]

\[
v_t=\frac{\beta\rho_lk_r}{\rho_{cr}}\theta(1).
\]

通解写为

\[
\theta=A\cosh(\lambda\xi)+B\sinh(\lambda\xi).
\]

消去 \(A,B\) 后，顶部梯度满足精确恒等式

\[
\theta'(1)=\lambda\tanh\lambda\,\theta(1)+\theta'(0)\operatorname{sech}\lambda.
\]

定义

\[
A_d=\frac{\rho_l\epsilon_sD_{s,eff}\lambda}{h},\quad
K=\beta\rho_l(1-\epsilon_s)(1-C_{eq})k_r,
\]

则 \(K/A_d=\sqrt{Da^*}=X\)。代入两端边界得到

\[
\theta(1)=\frac{C_{eq}j-B_{soil}\operatorname{sech}\lambda}
{A_d\tanh\lambda+K}.
\]

令 \(u=B_{soil}/(C_{eq}j)\)，最终得到

\[
\boxed{
v_t=\frac{C_{eq}j}{\rho_{cr}(1-\epsilon_s)(1-C_{eq})}
\frac{X[1-u\operatorname{sech}(\sqrt{Da})]}
{X+\tanh(\sqrt{Da})}
}.
\]

因此，任务中给出的新公式在上述符号和边界定义下成立。

### 3. 与 Prat 有限-\(Da\) 问题的关系

本问题的 \(u=0\) 是 \(\theta'(0)=0\)；Prat 是 \(\theta(0)=0\)。Prat 解给出

\[
v_{Prat}=v_{ref}
\frac{X\tanh\lambda}{1+X\tanh\lambda}
=v_{ref}\frac{X}{X+\coth\lambda},
\]

而新问题在 \(u=0\) 时是

\[
v_{N,0}=v_{ref}\frac{X}{X+\tanh\lambda}.
\]

二者在有限 \(Da\) 下一般不相等。若把底部通量设为依赖解的特殊值

\[
u_{Prat}=\frac{1}{\cosh\lambda+X\sinh\lambda},
\]

通量边界问题可精确复现 Prat；这不是预先给定常数 \(u\) 的问题。另一个共同极限是 \(Da\to\infty\)，此时 tanh 和 coth 都趋于 1。

### 4. 大-\(Da\) 弱耦合

同一通量边值问题内有精确比值

\[
\boxed{\frac{v_t(u)}{v_t(0)}=1-u\operatorname{sech}(\sqrt{Da})}.
\]

这不是渐近式。大 \(Da\) 时才可进一步写为

\[
1-2u e^{-\sqrt{Da}}+O(e^{-3\sqrt{Da}}).
\]

在最强测试 \(u=1.25\) 下：

| \(Da\) | 顶部速度降低 |
|---:|---:|
| 0.1 | 119.0%，速度反向 |
| 1 | 81.01% |
| 4 | 33.23% |
| 10 | 10.56% |
| 20 | 2.855% |
| 45 | 0.3052% |
| 100 | 0.01135% |

若把“弱耦合”操作性定义为测试范围内 \(u\le1.25\) 的速度变化低于 1%，则需要

\[
Da>[\operatorname{arcosh}(1.25/0.01)]^2=30.49.
\]

所以“\(Da>10\) 即可忽略底部扩散”并不由本结果支持。\(Da\approx20\) 是弱但仍可量化的耦合；\(Da\ge45\) 才在本测试范围内降到 0.31% 以下。

### 5. 3.5 μm 实际参数切片

代码实际使用

\[
C_{eq}=0.26494087,\quad \epsilon_s=0.43,\quad
D_{s,eff}=\sqrt{\epsilon_s}(1.3\times10^{-9})
=8.52467\times10^{-10}\ \mathrm{m^2s^{-1}},
\]

\[
a_v=36(1-\epsilon_s)/(25\ \mu\mathrm m)=820800\ \mathrm{m^{-1}},
\]

\[
k_r=10^{-3}\ \mathrm{m\,s^{-1}},\quad h=3.5\ \mu\mathrm m,
\quad \beta=2.
\]

由此

\[
\boxed{Da=20.16275,\qquad Da^*=3.17504,\qquad X=1.78186},
\]

\[
\ell_r=h/\sqrt{Da}=0.77946\ \mu\mathrm m.
\]

以归档算例的 \(j=1.7815\times10^{-5}\ \mathrm{kg\,m^{-2}\,s^{-1}}\) 代入，\(u=0\) 的新边界问题给 \(v_t=3.3332\times10^{-9}\) m/s；Prat 给 \(3.3326\times10^{-9}\) m/s。两者在该 \(Da\) 下很接近，是因为 tanh/coth 已接近 1，不是边界条件相同。对 \(u=1.25\)，速度降低 2.804%，到 \(3.2397\times10^{-9}\) m/s。

### 6. 引入底部退缩前的固体与溶解盐平衡

对剩余盐壳体积，固体账本必须写成

\[
\boxed{
\frac{dM_{cr}}{dt}=G_{top}-D_{base}
+\int_{z_b}^{z_t}R_{bulk}(z)\,dz
}.
\]

以向上为正的相对盐通量 \(S^{rel}\) 定义移动控制体的溶解盐平衡：

\[
\frac{dM_{diss}}{dt}=S_b^{rel}-S_t^{rel}
-\int_{z_b}^{z_t}R_{bulk}(z)\,dz.
\]

界面反应需要出现在相应跳跃条件中。若底部反应产物有比例 \(\alpha\) 进入盐壳孔隙，则 \(\alpha D_{base}\) 进入盐壳侧边界通量，其余 \((1-\alpha)D_{base}\) 进入土壤侧账本。

关键限制：现有 M0 的负 \(R_{bulk}\) 已使 \(M_{cr}\) 增加一个负贡献。若再把同一底部网格区的 \(-R_{bulk}\) 积分定义成 \(D_{base}\)，会重复扣除晶体。本原型只有在体反应排除层之外才允许独立 \(D_{base}\)。

### 7. M0 与 M1

| 模型 | 底部 | 体反应 | 物理解释 | 本阶段状态 |
|---|---|---|---|---|
| M0 原 ALE | \(z_b=0\) | 有符号、全域 | 溶解表现为孔隙率增加；无几何退缩 | 已验证基线 |
| M1 双界面 | \(z_b(t)\) | 有符号但排除界面区 | 独立界面溶解导致底面上移；体反应只改变剩余体积内部孔隙/固体 | 规定速率诊断原型 |

底部溶解究竟应表现为孔隙率增加还是底面上移，不是数学离散自动决定的：

- 弥散、体积分布的欠饱和反应应进入 \(R_{bulk}<0\)，主要增加孔隙率；
- 集中在清晰底部界面、且能定义界面材料转化的反应才可进入 \(D_{base}\)，产生 \(v_b>0\)；
- 二者可以共存，但反应区必须互斥或采用守恒重分类。

当前数据没有给出可区分这两类机制的底部界面动力学或成像证据，因此不能选择唯一闭合。

### 8. 双移动边界 ALE

映射

\[
z=z_b+\xi h
\]

给出仿射网格速度

\[
\boxed{w(\xi)=v_b+\xi(v_t-v_b)}.
\]

定义 \(H=-w\epsilon_{upwind}\)，则孔隙体积守恒为

\[
\frac{dE_i}{dt}=H_{i-1/2}-H_{i+1/2}
-\frac{R_iV_i}{\rho_{cr}}.
\]

边界值为

\[
H_b=-v_b\epsilon_b,\qquad H_t=-v_t\epsilon_t.
\]

这保证

\[
\frac{d}{dt}\sum_iE_i=-v_b\epsilon_b+v_t\epsilon_t
-\frac{1}{\rho_{cr}}\int R_{bulk}dz.
\]

由于 \(\dot h=v_t-v_b\)，直接从 \(M_{cr}=\rho_{cr}(h-\sum E_i)\) 得到

\[
\frac{dM_{cr}}{dt}
=\rho_{cr}(1-\epsilon_t)v_t
-\rho_{cr}(1-\epsilon_b)v_b+\int R_{bulk}dz.
\]

因此候选几何关系

\[
v_t=\frac{G_{top}}{\rho_{cr}(1-\epsilon_t)},\qquad
v_b=\frac{D_{base}}{\rho_{cr}(1-\epsilon_b)}
\]

只有在界面新生/移除材料的固相体积分数确为 \(1-\epsilon_t\) 和 \(1-\epsilon_b\) 时成立。

溶解盐控制体为

\[
\frac{dQ_i}{dt}=S_{i-1/2}^{rel}-S_{i+1/2}^{rel}-R_iV_i,
\]

其中

\[
S^{rel}=J_s-w\rho_l\epsilon C.
\]

同一个 \(w\) 同时进入孔隙体积和盐的相对通量，这是几何守恒的必要条件。只替换内部网格速度而不重写两端边界通量是不守恒的。

## Numerical verification

### 1. 解析解的独立数值验证

对任务指定的 7 个 \(Da\)、6 个 \(u\)、4 个 \(X\)，共 168 个组合，使用 500 区间二阶有限差分直接求解下端 Neumann、上端 Robin 问题。顶部速度因子和全剖面相对解析解的最坏绝对误差为 \(1.8504\times10^{-4}\)。精确比值 \(1-u\operatorname{sech}\sqrt{Da}\) 在全部组合中同时由闭式因子和数值解复现。

### 2. A--G 测试

| 测试 | 结果 | 关键量 |
|---|---|---|
| A M0 恢复 | PASS（算子级） | \(v_b=D_{base}=0\)，固定底部孔隙算子误差小于 \(10^{-18}\) |
| B 纯顶部生长 | PASS | 1000 s 内 \(h:100\to102\) μm |
| C 纯底部退缩 | PASS | \(z_b:0\to1\) μm，\(h:100\to99\) μm |
| D 等前沿速度 | PASS | \(v_t=v_b\)，\(h\) 恒定；因 \(\epsilon_t\ne\epsilon_b\)，晶体库存仍改变 \(1.5155\times10^{-4}\) kg m\(^{-2}\) |
| E 退缩快于生长 | PASS | \(z_t\uparrow,z_b\uparrow,h:100\to99.5\) μm |
| F 有符号体反应 + 界面溶解 | PASS | 底部排除区无体反应；界面和体账本独立 |
| G 网格/时间步 | PASS（诊断收敛） | \(N=40,80,160\)，\(\Delta t=0.5,0.25\) s |

A--F 最大绝对盐守恒残差为 \(1.90\times10^{-14}\) kg m\(^{-2}\)，最大绝对水守恒残差为 \(3.34\times10^{-15}\) kg m\(^{-2}\)。

G 中，以 \(N=160,\Delta t=0.25\) s 为参考，固体库存绝对差从 \(N=40,\Delta t=0.5\) s 的 \(4.39\times10^{-9}\) 降至 \(N=80,\Delta t=0.5\) s 的 \(1.05\times10^{-9}\) kg m\(^{-2}\)；同网格将 \(\Delta t\) 从 0.5 降到 0.25 s 产生更小的独立时间步差。盐/水守恒残差均保持在约 \(10^{-15}\) kg m\(^{-2}\) 的浮点水平。

### 3. 基线回归证据边界

本机没有 C++ 编译器或 SUNDIALS，因此没有声称重新运行生产 C++ BDF/DuMuX 回归。归档的 BUILD_VERIFICATION.txt 报告：

- flux_quadrature_test：PASS，最坏相对局部盐平衡残差 \(3.31591\times10^{-15}\)；
- progressive_test：PASS。

当前 M0 的 ale.hh、bdf.hh、flux_quadrature_test.cc、progressive.hh 和 progressive_test.cc 的 SHA-256 与该归档构建输入完全一致。该证据支持“起始基线与已验证归档相同”，但不等价于本轮本机重编译。

## Figures and reproducibility

代码：

- model/p2_h_dual_fronts/analytical.py
- model/p2_h_dual_fronts/dual_front.py
- model/p2_h_dual_fronts/run_verification.py
- model/p2_h_dual_fronts/test_dual_front.py

结果目录：analysis/p2_h_dual_moving_fronts_20261112/。

主要图件：

- normalized_top_growth_vs_u.png
- growth_reduction_vs_Da.png
- concentration_profiles.png
- signed_internal_reaction_profiles.png
- finite_Da_vs_asymptotic_growth.png
- penetration_length_ratio.png
- dual_front_geometry_and_inventory.png
- grid_time_convergence.png

Each PNG has a vector PDF companion. No target journal or final print width was
specified, so these are provisional general scientific figures rather than a
claim of publisher compliance.

复现命令：

    python -m unittest model.p2_h_dual_fronts.test_dual_front -v
    python model\p2_h_dual_fronts\run_verification.py

## Answers to the three scientific questions

### 问题一：大 \(Da\) 时顶部生长是否真的对 \(u\) 不敏感？

**条件性成立，不是普遍成立。** 精确敏感度为 \(u\operatorname{sech}\sqrt{Da}\)。对 \(u\le1.25\)，\(Da=20\) 仍有最高 2.86% 影响；\(Da=45\) 降到 0.305%；\(Da=100\) 为 0.011%。解析解和独立有限差分解一致。小 \(Da\)、大 \(u\) 甚至会使顶部速度变负，表明继续生长假设失效。

### 问题二：底部溶解表现为孔隙率增加还是底面上移？

**两者是不同物理闭合，现有证据不能唯一决定。** 空间分布式溶解属于 \(R_{bulk}<0\)，表现为剩余盐壳孔隙率增加；集中界面溶解属于 \(D_{base}>0\)，表现为 \(z_b\) 上移。必须用反应区尺度、界面材料转化和实验/微观证据区分，不能把底部若干体单元的负反应再次作为 \(D_{base}\)。

### 问题三：能否守恒地实现顶部生长、底部退缩甚至整体上移？

**数值上可以，且 A--G 诊断测试已经证明库存守恒；物理预测上尚不可以。** 仿射网格速度和两端相对通量允许 \(v_t>0,v_b>0\)，包括 \(v_t=v_b\) 的定厚平移和 \(v_t<v_b\) 的上移变薄。缺失的是 \(D_{base}\)、\(\alpha\) 和移动界面土壤跳跃条件的独立闭合及全土壤--盐壳联合验证。

## Remarks and Interpretation

- 大 \(Da\) 使底部通量对顶部速度的直接影响指数衰减，但不使底部固体损失消失。因此顶部速度几乎不变与总晶体库存显著减少可以同时发生。
- 盐壳定厚平移不保证晶体库存不变；当两端材料孔隙率不同，即使 \(v_t=v_b\)，\(G_{top}\) 与 \(D_{base}\) 也可不同。
- 顶部高度 \(z_t\) 单独增加不能证明盐壳变厚。必须同时报告 \(z_b,h,M_{cr},\bar\epsilon\)。

## Boundaries and Non-Claims

- 未把 \(D_{base}=uS_{ref}\) 当作普遍定律；代码中没有此闭合。
- 未把负体反应自动转成底部退缩。
- 未把规定速率原型接入生产 DuMuX 土壤求解器。
- 未证明移动底部的真实水力连续性、盐活度、界面曲率或孔隙拓扑闭合。
- 未重新运行本地 C++/SUNDIALS 基线；只核对了与已通过归档完全相同的源文件哈希。
- 当前外部储库账本代表土壤侧/供给侧的诊断守恒项，不是已验证的土壤本构关系。

## Open Risks and next steps

1. 从独立物理机制建立 \(D_{base}(C_b,S_w,p,\ldots)\) 或清晰的规定速率实验协议；
2. 决定底部溶解产物的向上/向下分配 \(\alpha\)，并用移动界面跳跃条件替代诊断常数；
3. 在生产 C++ 中加入默认关闭的双前沿状态 \(z_b,z_t\)，并重新推导顶部和底部边界通量；
4. 在 HPC 工具链上重跑 M0 回归、局部通量积分和 M1 新测试；
5. 最后才进行土壤--盐壳联合守恒测试，并检查移动底面是否需要重映射或释放土壤控制体。
