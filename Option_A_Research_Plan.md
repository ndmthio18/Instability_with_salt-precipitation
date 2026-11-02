# Option A 研究实施方案

**建议工作题目：** *From Linear Instability to Effective Salt Redistribution in Evaporating Unsaturated Porous Media*

**研究背景：** 本方案基于 Chen、Ciriello 和 Helmig 合作的 Paper 1《Evaporation-Driven Density Instability in Unsaturated Porous Media: Observable Onset and Competition with Salt Precipitation》、2026 年 11 月 1 日报告《Transfer Learning Results and Proposed Extensions of Evaporation-Driven Instability Modelling》的 Option A，以及随后关于研究创新性、DuMux 盐输运和 LSA-informed prediction 的讨论。

**时间定义更新（2026-11-02）：** 两种 LSA onset、DuMuX `t_detect` 与独立的 `t_eff` 严格区分，见 [统一约定](analysis/onset_semantics_20261102/CONTRACT_ZH.md)。下文输运阈值和时间线属于待验证的候选设计，当前未定义 `t_eff`，没有训练标签。

**研究策略：** 小规模可行性验证 → 定量物理机制研究 → LSA 预测能力验证 → 代理模型开发 → 稳健性和论文整合。

> **核心原则：** 第一阶段不是证明 Option A 一定成立，而是用少量计算尽早判断，它是否包含独立、可靠、值得发表的新科学问题。DuMux 是非线性盐输运的主要计算依据；LSA 提供 onset、增长率、模态和候选的线性输运指标；机器学习仅在物理目标明确之后用于建立可推广的预测关系。

---

## 一、研究总体目标

**核心科学问题：**

> Under what conditions does evaporation-driven density instability produce consequential downward salt redistribution before the near-surface salt concentration reaches the solubility limit?

即：在什么样的土壤水力性质、蒸发条件和地下水位条件下，密度不稳定性能够在表层盐浓度达到溶解度极限之前，产生足够显著的向下盐分输运？

### 三个研究问题（Research Questions）

- **RQ1 — From instability onset to effective transport：** 线性不稳定性出现之后，需要经历怎样的非线性发展过程，才会产生显著的向下盐输运？这种延迟是否具有明确的物理规律？
- **RQ2 — Controls on transport effectiveness：** 水饱和度、渗透率、蒸发强度、地下水位和不稳定模态，如何共同控制盐输运效率？
- **RQ3 — Predictability from linear stability information：** 能否利用 LSA 提供的 onset、增长率和模态信息，预测 DuMux 中的非线性盐输运，而不需要为每一种条件都运行昂贵的二维模拟？

**论文定位：** RQ1 和 RQ2 是论文主体（physical understanding）；RQ3 是与 Valentina 合作的预测方法贡献（predictive framework）。

### 分阶段路线与起始规模

| 阶段 | 核心任务 | 建议计算规模 | 主要产出 |
|---|---|---|---|
| Phase 0 | 现有数据、代码、方法与文献审计 | 优先利用已有 436 次运行 | 数据和模型审计报告 |
| Phase 1 | Feasibility Pilot：有效盐输运是否可定义、可观测？ | 约 12–20 组物理条件，先做 2 组 Pilot-0 | 可行性报告、GO/NO-GO |
| Phase 2 | Systematic Physical Investigation | 初步规划 40–80 组物理条件 | 机制分解、定量规律、区域图 |
| Phase 3 | LSA-Informed Prediction | 利用 Phase 2 数据，按需追加 | 预测模型、学习曲线 |
| Phase 4 | Validation & Publication | 少量独立验证与必要补算 | 论文完整证据链 |

上述样本量是起始预算，非硬性规模。是否扩展必须依据前一阶段的物理信号、数值可靠性和计算成本。

---

## 二、Phase 0 — 研究基础与现有结果审计

### 0.1 已有基础

Paper 1 已经具有：

- 43,253 个有效 LSA 案例；9,354 个 instability-first 案例。
- 临界起始时间 $t_{\mathrm{crit}}$、observable onset（增长筛选模态自身的中性穿越时间） $t_{\mathrm{seed}}(0.10)$ 与模态信息。
- 采用特定初始扰动和检测协议的 DuMux onset 对比验证。
- 固定配置：$H=0.2\,\mathrm m$，横向宽度 $W=3H=0.6\,\mathrm m$。
- 表面盐浓度达到溶解度阈值 $x_{\mathrm{sat}}=0.1$（摩尔分数）的竞争事件定义。

2026 年 11 月 1 日的报告另汇总了 436 次 DuMux 运行、253 个物理案例和 293 个用于残差分析的配对观测。**这些数字代表不同统计单位，不应混用。**

现有模拟主要面向 onset detection，不一定保存了足够长的非线性输运演化。

### 0.2 数据审计清单

- [ ] **浓度场：** 是否有完整的二维盐浓度时序数据（VTK 等）？
- [ ] **水流与饱和度：** 是否保存了 Darcy velocity、water saturation？
- [ ] **守恒盐通量：** 是否可提取不同控制面上的 advective、diffusive 和 total salt flux？
- [ ] **时间覆盖：** 每次模拟在 onset 后持续了多久？是否出现充分发展的 fingers？
- [ ] **LSA 元数据：** 能否匹配 $t_{\mathrm{crit}}$、$t_{\mathrm{seed}}$、模式、增长率和初始扰动？
- [ ] **停止条件：** 是 90 天 horizon、表面溶解度阈值，还是数值失败造成停止？

优先提取 DuMux **离散控制面守恒通量**，避免仅通过速度和浓度场插值造成额外误差。若只能重构通量，需要另行校验重构误差。

### 0.3 与已有工作的创新性边界

已有研究（尤其 Bringedal et al., 2025）已经研究非饱和多孔介质中的 *linear growth*, *flux growth*, *finger merging* 等阶段，并使用表层平均浓度峰值指示对流对表层盐积累的影响。

因此，Option A 不应仅重复“表层盐浓度何时开始下降”。必须争取新增：

1. **盐分输运的空间和累计量化：** 有多少盐跨越不同深度、从近表层转移到内部？
2. **与溶解度阈值的竞争：** 额外向下输运是否在表层到达溶解度阈值之前达到有意义的程度？
3. **机制与可预测性：** 是否存在可解释、可验证的控制参数、无量纲规律和 LSA-informed predictor？

**Phase 0 交付物：** 简短审计报告，包括可用数据、缺失数据、可复用模块、必要新增输出、已知文献与潜在新贡献的对照。

---

## 三、Phase 1 — Feasibility Pilot

### 阶段目的

检验：**能否在少量 DuMux 案例中，可靠识别线性 onset 和具有实际意义的盐输运效果之间的区别？** 暂时不做机器学习。

### 1.1 成对建立 1D–2D DuMux 数值实验

**M0 — 1D reference simulation**

- 使用完整非线性 Richards–salt 模型，不引入横向不稳定性。
- 输入与二维模型一致：土壤水力参数、蒸发需求、底部水位/压力、初始盐浓度、组分与物性关系。
- 输出背景盐浓度、盐库存、盐通量、表面达到溶解度阈值的时间。

**M1 — 2D nonlinear simulation**

- 保留 Paper 1 的二维计算域与物理模型。
- 按照原论文协议，在初始时刻（不是在 $t_{\mathrm{seed}}$ 时）施加幅值 $10^{-6}$ 的 LSA 单模态浓度扰动。
- 输出非线性指状结构、盐库存、盐通量、表面浓度变化。

**额外数值一致性对照：** 选少量算例，运行横向均匀初始条件的 2D 模型，并与 1D 结果比较，以排除维度与离散误差造成的假信号。

### 1.2 代表性案例选择（起始 16 组）

| 类别 | 数量 | 选择依据 | 目标 |
|---|---:|---|---|
| P1：明显不稳定 | 4 | $t_{\mathrm{seed}}$ 远早于基态溶解度阈值 | 是否存在充分发展的盐输运 |
| P2：竞争边界 | 4 | $t_{\mathrm{seed}}$ 接近 $t_{\mathrm{precip}}$ | 对流是否赶得上表面盐积累 |
| P3：模态差异 | 4 | $k_{\mathrm{seed}}/k_{\mathrm{crit}}$ 差异明显 | 波数比与输运效率关系 |
| P4：水力性质对照 | 4 | 不同饱和度、渗透率和土壤类别 | 跨条件的物理规律 |

这是选样目标，不应机械地把波数比与土壤类别等混杂因素混为因果关系。必要时在同一个物理案例中采用不同模态初始化进行对照。

### 1.3 三类主要诊断量

#### A. 起始和阶段时间

已有：

$$
t_{\mathrm{crit}},\qquad t_{\mathrm{seed}},\qquad t_{\mathrm{detect}}.
$$

新增表层平均浓度峰值时间：

$$
t_{\mathrm{peak}}=\operatorname*{arg\,max}_t\overline{x}^{2D}_{\mathrm{top}}(t).
$$

此量用于与前人 flux-growth 研究比较，不应单独作为创新点。

#### B. 额外盐通量与累计输运

采用 Paper 1 的竖向坐标（底部 $z=0$，表面 $z=H$），初始观测深度：

$$
\frac{z^*}{H}=0.9,\qquad 0.75,\qquad 0.5.
$$

统一设**向下为正**，将二维与一维的总盐通量按照一致截面积进行积分换算：

$$
\Delta Q_s(z^*,t)=Q_s^{2D}(z^*,t)-Q_s^{1D}(z^*,t).
$$

累计额外盐输运：

$$
\Delta M_s(z^*,t)=\int_0^t\Delta Q_s(z^*,t')\,dt'.
$$

无量纲指标：

$$
\boxed{\mathcal T(z^*,t)=\frac{\Delta M_s(z^*,t)}{M_{\mathrm{ref}}}.}
$$

初步可将 $M_{\mathrm{ref}}$ 定义为指定深度以上区域的初始盐库存，并在全部案例中一致使用。实际数值处理中必须使库存和盐通量使用相同的守恒量与单位；不能直接将盐摩尔分数视为盐质量。

#### C. 表层盐积累的实际响应

指定深度以上区域的盐库存差异：

$$
\Delta M_{\mathrm{upper}}(t)=M_{\mathrm{upper}}^{2D}(t)-M_{\mathrm{upper}}^{1D}(t).
$$

同时保存二维表面**平均**和**局部最大**盐浓度，并记录：

$$
t_{\mathrm{sat}}^{1D},\qquad t_{\mathrm{sat}}^{2D}.
$$

$t_{\mathrm{sat}}^{1D}$ 表示一维参考表面首次达到溶解度阈值的时间；$t_{\mathrm{sat}}^{2D}$ 表示二维表面任一点首次达到溶解度阈值的时间。

两者均不是完整盐结晶动力学预测。原模型没有沉淀过程，若达到溶解度阈值，应记录适用边界，不应将长期过饱和的无沉淀模拟当作可靠的后续结晶预测。

### 1.4 Transport-effective onset 的探索性定义

先比较候选阈值：

$$
\mathcal T(z^*,t)=1\%,\quad 3\%,\quad 5\%.
$$

这些阈值只用于初期敏感性检验，**尚无已验证的物理最优值**。正式定义应满足：

- 超过离散数值噪声和守恒误差水平。
- 对观察深度和阈值的合理变化保持主要规律稳健。
- 能反映持续的实际输运，而不只是瞬时跨阈值。
- 与近表层盐库存和溶解度阈值竞争形成清晰联系。

**不要预先强制 $t_{\mathrm{eff}}\ge t_{\mathrm{seed}}$**。先分析从模拟开始时的绝对累计输运演化，再研究相对延迟；对从未跨越阈值的案例保留右删失信息。

### 1.5 Pilot 必须给出的四组图

1. **F1 — 完整演化时间轴：** $t_{\mathrm{seed}}$、$t_{\mathrm{detect}}$、$t_{\mathrm{peak}}$、$t_{\mathrm{eff}}$、表层浓度和指状流动快照。
2. **F2 — 1D/2D 盐库存与通量对照：** 不同深度以上的盐库存、额外累计输运与守恒闭合。
3. **F3 — Transport delay comparison：** $t_{\mathrm{eff}}-t_{\mathrm{seed}}$ 在不同物理条件下的变化，包含未达到阈值的删失案例。
4. **F4 — Event competition：** Paper 1 的 instability-first 分类与新的 transport-effective-before-saturation 判断之间的区别。

### 1.6 GO / NO-GO 标准

- [ ] **G1 数值可靠性：** 1D 与横向均匀 2D 的背景解一致，离散盐收支误差显著小于额外输运信号。0.5% 可作为待校准的相对守恒残差初始目标，但不能替代针对信号量级的误差判断。
- [ ] **G2 物理信号明确：** 若干不同条件下存在明显高于数值噪声的额外向下盐输运，并能在溶解度阈值前观察其发展。
- [ ] **G3 指标稳健：** 主要结论对观测深度、阈值、持续性与合理的网格/时间步加密不发生根本变化。
- [ ] **G4 有独立科学增量：** 除既有表层浓度峰值判据外，新指标揭示累计输运、空间穿透或事件竞争的实质性新信息。

**GO：** 四项均有证据支撑，进入 Phase 2。

**条件性 NO-GO：** 若数值可靠性不足，先修复模型和后处理；若输运信号微弱或缺乏相较已有文献的创新增量，先调整指标或问题，不盲目扩大算例。

---

## 四、Phase 2 — 系统的非线性盐输运机制研究

### 2.1 分解输运机制

对二维速度和浓度做水平平均与扰动分解：

$$
q_z=\overline{q}_z+q_z',\qquad x=\overline{x}+x'.
$$

对流盐通量的水平平均：

$$
\boxed{\overline{q_zx}=\overline{q}_z\,\overline{x}+\overline{q'_zx'}.}
$$

- $\overline{q}_z\,\overline{x}$：平均流动携带的盐分。
- $\overline{q'_zx'}$：流速扰动与浓度扰动之间的相关输运。

同时分别量化扩散项及二维平均浓度剖面的变化。**$Q_s^{2D}-Q_s^{1D}$ 是总的净差异，并不严格等于 $\overline{q'_zx'}$**。

核心机制研究链条：

> DuMux 计算真实扰动相关输运 → LSA 提出线性阶段输运 proxy → 比较一致和偏离的时间、条件与机制。

### 2.2 系统性实验矩阵

继续使用 Paper 1 的七个无量纲输入：

$$
\boldsymbol\xi=(R_E,\beta,\gamma,\widehat P_B,n,S_{wr},L).
$$

实验分为：

**A. 真实参数组合：** 从 LSA 数据库中选择代表性的联合土壤参数、蒸发和地下水位组合，避免把相关水力参数当作独立变量随意采样。

**B. 受控机制实验：** 在同一参考配置中单独控制少数参数或扰动模式，用于区分因果机制。偏离真实土壤联合分布的组合应明确标记为机制实验。

初步计划约 40–80 组物理条件，根据 Phase 1 结果和可用计算资源调整。

### 2.3 三组竞争机制

1. **Instability growth vs. salt accumulation：** 扰动是否足够快地增长，还是表面已经逼近溶解度极限？
2. **Convective transport vs. diffusive transport：** 指状对流是否显著改变总盐通量？饱和度及相对渗透率如何影响输运？
3. **Salt redistribution vs. surface concentration increase：** 非线性盐输运能否在局部达到溶解度阈值之前明显降低近表层盐库存或改变其分布？

### 2.4 候选 Transport-effectiveness regime diagram

| 候选区域 | 失稳表现 | 输运表现 |
|---|---|---|
| Regime I | 无可检测有效失稳 | 背景蒸发与扩散主导 |
| Regime II | 不稳定性出现并增长 | 额外向下输运仍较弱 |
| Regime III | 显著非线性对流 | 盐分有效向下重新分布 |
| Regime IV | 表面先达到溶解度极限 | 无沉淀模型的后续适用性受限 |

这是**待验证的候选分区**，不可强行假设一定存在清晰边界。

探索是否能用 $R_E/\beta$、局部饱和度、增长率特征、波数和输运时间尺度建立可解释的判据。区分局部条件规律与跨土壤统计规律。

**Phase 2 交付标准：**

- 经过守恒和阈值稳健性验证的主要输运指标。
- 对不同土壤条件下输运差异的机制解释。
- 至少一个经独立案例检验的定量关系、条件性规律或区域图。
- 明确对照已有 flux-growth 研究，证明新增的是实际累计输运及其竞争后果，而不是已有阶段判据的重新命名。

---

## 五、Phase 3 — LSA-informed 预测框架

### 3.1 新的预测任务

此前从 $t_{\mathrm{seed}}(0.10)$ 迁移学习 DuMux detected onset 的收益有限，因为两者在匹配协议下差距较小。

本阶段不是再修正一个已经很准确的 onset 时间，而是预测非线性输运结果：

$$
\left(t_{\mathrm{eff}},\mathcal T,\Delta M_s\right).
$$

### 3.2 四种模型对照

| 模型 | 输入 | 对照目的 |
|---|---|---|
| M1 — Direct surrogate | 物理参数 $\xi$ | 最基础、可解释的预测基线 |
| M2 — LSA-informed | $\xi$ + onset、模态、增长率等 | 测量 LSA 物理信息的增量价值 |
| M3 — Flux-proxy assisted | $\xi$ + LSA 扰动相关输运 proxy | 检验线性通量近似是否有用 |
| M4 — Transfer learning | LSA 预训练表示 + 少量 DuMux 输运数据 | 检验稀缺高保真数据下的迁移收益 |

建议优先使用 Paper 1 已验证的多项式代理方法作为 direct baseline。只有当直接代理与 LSA 特征增强方案不足以达到目标时，再重点探索 transfer learning。

**删失处理：** 对达到溶解度阈值前尚未产生显著有效输运的案例，不应删除或误记为一个确定的输运时间。可考虑事件分类 + 条件回归，或适当的生存分析方法。

### 3.3 学习曲线与泛化检验

在数据允许时，采用例如 20、40、60、80 组 DuMux 训练样本比较模型学习曲线；所有模型使用同一批独立物理案例作为测试集。

重点指标：

- 预测误差（对数误差、相对误差、必要的 $R^2$）。
- 不同土壤类别与不同不稳定机制的泛化能力。
- 未发生显著输运案例的分类与删失建模质量。
- 达到规定误差所需的 DuMux 样本量。
- 完整计算成本：LSA 特征计算、DuMux 训练数据、模型训练与预测。

**Phase 3 交付标准：** 量化 LSA 是否对非线性盐输运的预测提供实质性增益。若没有，论文的物理机制贡献仍可独立成立。

---

## 六、Phase 4 — 稳健性、实验对照与论文整合

### 4.1 数值稳健性

- [ ] 网格和时间步加密后的 $t_{\mathrm{eff}}$、$\Delta M_s$ 收敛。
- [ ] 累计盐通量与对应上部盐库存变化满足离散守恒。
- [ ] 调整初始扰动幅值、模式或多模态扰动后，主要物理判断仍可解释。
- [ ] 横向计算域宽度变化下，模式选择与输运响应的变化得到识别。
- [ ] 观察深度、输运阈值和持续性准则下的主要结论稳健。

应区分**对扰动的物理稳健性**与**LSA 对既定扰动协议的预测精度**：它们是不同问题。

### 4.2 独立实验支持

可考虑 Chaudhry et al. (2026) 的 $^{23}\mathrm{Na}$-MRI 盐分重新分布实验及相应数值模型，验证盐分重新分布的部分机制或浓度剖面。

- 论文 DOI：https://doi.org/10.1029/2025WR041207
- MRI 数据：https://doi.org/10.18419/DARUS-5094
- 相关 DuMux 模型：https://doi.org/10.18419/DARUS-4744

**限制：** 对应实验主要是饱和条件，不能作为非饱和研究的直接全面验证；需要明确适用范围。

### 4.3 候选论文结构

| 章节 | 主要内容 | 科学价值 |
|---|---|---|
| Introduction | 为什么线性 onset 不能充分描述实际盐输运 | 提出新问题 |
| Methods | LSA、1D/2D DuMux、输运指标、样本设计 | 一致且可复现的方法 |
| Transport evolution | 从线性扰动到显著向下输运 | 非线性机制 |
| Transport regimes | 土壤条件、模态、输运与溶解度阈值竞争 | 定量规律 |
| Predictive framework | LSA-informed 与 direct surrogate 的比较 | 预测贡献 |
| Discussion | 模型边界、实验联系、盐渍化意义 | 应用价值 |

**期刊定位：** 若只有新 onset 定义和少量相关性分析，创新性可能不足以支撑 AWR；若建立可靠的输运机制、跨土壤规律和有价值的 LSA-informed 预测，AWR 是合理目标；若进一步建立与表层积盐和水文应用的强联系并得到充分验证，亦可考虑 WRR。

---

## 七、立即执行：Pilot-0（优先做两组案例）

为了尽快验证方向，建议先从最小可执行任务开始：

### Pilot-0 / Case A — 明显发生非线性指状输运

1. 从现有运行中寻找一个 onset 明显早于参考溶解度阈值，且保存完整浓度场的砂土案例。
2. 运行匹配的 1D DuMux 参考解，并检查横向均匀 2D 解与 1D 的一致性。
3. 必要时延长 2D 运行（仅在物理模型适用的阈值范围内），提取完整盐通量。
4. 画出 1D/2D 表层浓度、额外输运与累计盐库存时间曲线，标出 LSA onset 和表面浓度峰值。

### Pilot-0 / Case B — 接近表面溶解度阈值竞争边界

1. 选择一个 $t_{\mathrm{seed}}$ 接近 $t_{\mathrm{precip}}$ 的案例。
2. 同样建立 1D/2D 对照，比较是否有足够时间形成有效输运。
3. 对比 Case A 和 B 的模式、增长率、饱和度、盐库存与输运延迟。

### Pilot-0 最小交付物

- 一份包含模拟参数、模型版本、输出变量和停止条件的数据清单。
- 两个案例的 1D/2D 浓度和盐库存对照。
- 不同深度的累计额外盐通量及守恒检查。
- 一张包含 $t_{\mathrm{seed}}$、$t_{\mathrm{detect}}$、$t_{\mathrm{peak}}$、$t_{\mathrm{eff}}$、$t_{\mathrm{sat}}$ 的时间轴图。
- 对是否继续 Phase 1 的明确 GO / CONDITIONAL NO-GO 建议。

---

## 八、建议协作分工

| 研究任务 | 建议主导 |
|---|---|
| DuMux 模型与数值实验、盐输运机制分析 | Zhixin |
| LSA 特征提取、代理模型、迁移学习和不确定性量化 | Valentina |
| 物理机制讨论、既有模型经验、研究方向把关 | Rainer |
| 科学问题、论文创新点与结果解释 | 共同讨论 |

建议在 Pilot-0 和 Phase 1 的首批物理结果产生后，与 Valentina 和 Rainer 讨论，再决定是否扩大计算。

---

## 九、关键科学结论与决策逻辑

计划形成的科学链条：

$$
\underbrace{\text{Linear instability}}_{\text{Paper 1}}
\rightarrow
\underbrace{\text{Nonlinear fingering}}_{\text{DuMux}}
\rightarrow
\underbrace{\text{Effective salt transport}}_{\text{Option A}}
\rightarrow
\underbrace{\text{Near-surface salt response}}_{\text{Physical consequence}}.
$$

最后利用 LSA-informed surrogate 将可靠的输运规律推广至更多土壤条件。

**第一优先级：** 确认 DuMux 是否给出稳定、可量化、独立于数值误差、且相对已有研究具有新物理价值的 *effective transport signal*。

**继续扩展的前提：** Pilot-0/Phase 1 能够证明有效输运延迟、累计盐库存变化，或溶解度阈值竞争的有意义差异，并且这些结果经守恒和稳健性验证。

**如果未通过：** 优先调整诊断定义、物理模型或研究问题，不要直接扩大算例，也不要把开发机器学习模型作为替代。

---

## 十、主要参考资料

1. Chen, Z., Ciriello, V., & Helmig, R. *Evaporation-Driven Density Instability in Unsaturated Porous Media: Observable Onset and Competition with Salt Precipitation.* Paper 1。原方案引用名为 `onstime(4).pdf`；本轮定义以用户指定的合作稿 `AWR_onstime/onstime.pdf` Section 2.2–2.4 为准，文件哈希和字段规范见 [命名约定](analysis/onset_semantics_20261102/CONTRACT_ZH.md)。
2. *Transfer Learning Results and Proposed Extensions of Evaporation-Driven Instability Modelling.* 1 November 2026，用户提供的 `report2(1).pdf`，尤其 Section 7.1（Option A）。
3. Bringedal, C., Schollenberger, T., Pieters, G. J. M., van Duijn, C. J., & Helmig, R. (2022). *Evaporation-driven density instabilities in saturated porous media.* Transport in Porous Media. https://doi.org/10.1007/s11242-022-01772-w
4. Bringedal, C., Kiemle, S., van Duijn, C. J., & Helmig, R. (2025). *Impact of saturation on evaporation-driven density instabilities in porous media: Mathematical and numerical analysis.* Transport in Porous Media. https://doi.org/10.1007/s11242-025-02207-y
5. Chaudhry, M. A., Kiemle, S., Pohlmeier, A., Helmig, R., & Huisman, J. A. (2026). *Non-invasive imaging of solute redistribution below evaporating surfaces using ²³Na-MRI.* Water Resources Research. https://doi.org/10.1029/2025WR041207
