# 原机制与实现审查

## 原始机制

- `solvePrecipitation_()`：以表面固体通量为未知量求互补根；原边界总质量通量为 `Jevap+Jprec`，NaCl 通量为 `Jprec`。
- `precrust_()`：每个接受步把 `Jprec dt` 累加到单位面积表面固体库存。
- `top_only`：达到阈值后按顶部 BOX DOF 当前液体库存比例一次性取液。
- `adaptive_multilayer`：达到阈值后在 5 mm 范围内从上到下直接修改多个 BOX DOF，直至取得 ALE 初始孔液。
- `topDofs_()`：选择物理上边界的 BOX 顶部顶点；一水平单元时存在两个同高顶点，各自携带 BOX SCV 体积。
- 顶边界两个 SCVF 使用同一单位面积 `heldLiquid_` 与 `heldSalt_`；SCVF 面积求和等于物理柱宽。

## 原网格依赖来源

`top_only` 的可用质量与顶部 BOX 体积成正比；`adaptive_multilayer` 的混合浓度、donor 数量和实际深度由网格切分决定。二者均在阈值跨越时瞬时改变土壤状态，且激活浓度不是由此前真实界面通量历史形成。

## 新实现的隔离范围

- 新模式：`P2H.DonorMode=progressive_top_boundary_filling`。
- 仅允许与 `P2H.PrecipitationMode=implicit_surface` 联用。
- 增加固体、孔液盐、孔液水、所需/已填孔体积独立账本。
- 增加 `Jprec`、`Jfill,liquid`、`Jfill,salt` 独立输出。
- 蒸汽账本改为仅累加 `Jevap`，修正加入充液后旧差值表达式会产生的重复计水。
- 新模式激活时不调用 `top_only` 或 `adaptiveDonors_()`，不修改任何土壤 DOF。
- 原两种模式与历史路径保留，作为回归对照。

## 尚未验证

代码实现本身不证明网格收敛，也不证明真实毛细供液。必须先通过独立单元测试和 HPCC DuMuX pilot，再决定是否运行完整矩阵。激活后的 ALE 底部欠饱和、内部溶解与 `epsilon_max` 失败仍是独立问题。
