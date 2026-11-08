# P2-H 渐进前驱盐层公式与网格依赖预审计

## 1. 实施前结论

禁止使用 `available_top_cell_water/dt`。该式随顶部 BOX 控制体体积和时间步直接改变，因而不存在所要求的网格、时间步极限。旧 `top_only` 与 `adaptive_multilayer` 也把激活时的离散 DOF 库存直接映射为盐层孔液，初始化浓度会随网格、 donor 深度与 BOX 体积改变。

本阶段只实现一个**零前驱层水力阻力、饱和出生的诊断性运动学闭合**。它不是经验证的盐层毛细供液本构。其目的，是检验：若新生固体所需孔隙液与固体同步、且全部液体只通过 DuMuX 顶边界进入，初始化是否有离散极限。

## 2. 每单位表面积守恒关系

表面互补条件给出已收敛的固体析出通量 `Jprec`（kg m-2 s-1）。接受步内：

```text
dms       = Jprec dt
dh        = dms/[rho_crystal(1-epsilon_birth)]
dVpore    = epsilon_birth dh
dml       = rho(C_Gamma) dVpore
dmsalt,p  = C_Gamma dml
dmwater,p = (1-C_Gamma) dml
```

`C_Gamma` 和 `rho(C_Gamma)` 取自同一次已收敛的 DuMuX 顶边界状态。隐式表面析出互补条件在 `Jprec>0` 时给出 `C_Gamma=Ceq`；代码不在求解后覆盖浓度，而是检查实际离散顶边界状态与用于构造通量的状态一致。

土壤顶边界向外的两个质量通量为：

```text
J_total = Jevap + Jprec + Jfill,liquid
J_NaCl  = Jprec + Jfill,salt
```

因此土壤水损失为 `Jtotal-JNaCl = Jevap+Jfill,water`。蒸汽库存只累加 `Jevap`，前驱层分别累加 `Jfill,salt` 与 `Jfill,water`；同一份水不能同时进入蒸汽和前驱层。

## 3. 网格极限

运动学需求只依赖物理表面通量、材料密度、出生孔隙率和物理时间，不包含顶部 BOX 体积、donor 层数或 `1/dt` 可用量。因此若 DuMuX 的表面析出通量和界面状态收敛，该初始化具有连续极限。

这并不保证结果自动网格无关。剩余误差来自：

1. DuMuX 顶边界析出互补条件自身的空间离散误差；
2. Richards–盐输运对额外顶边界孔液通量的离散响应；
3. 时间分裂与非线性根求解容差；
4. 激活阈值的时间步越界误差。

这些误差必须由 `Nsoil=80/160/320` 和 `dt=0.5/0.25/0.125 s` 的实际收敛表判定，不能由公式形式代替。

## 4. 水力供给闭合的边界

当前模型没有前驱盐层毛细压、渗透率或界面动态润湿关系，因而不能从压差独立预测部分充液速率。本实现把“孔隙与固体同步饱和出生”作为诊断性 Neumann 需求，并让 Richards 方程隐式响应。只有在非线性求解收敛、界面密度一致、顶部状态通过声明的有效性筛查时，该步才接受；否则整步回滚，不从深层 DOF 补液，也不采用任意系数或等待时间。

因此，即使数值收敛，本阶段最多支持“该诊断性零阻力近似未被计算反驳”，不能支持“真实盐层供液机理已验证”。若供液筛查失败，当前物理模型未闭合；最小缺失关系是前驱层饱和度/毛细压力与界面液通量之间的动态水力关系。

## 5. ALE 激活

仅在真实固体等效高度达到 `3.5e-6 m`、前驱层所需和已填孔体积一致、顶部状态有效、且下列转移残差满足门槛时激活：

```text
Msolid,ALE = Msolid,precursor
Msalt,ALE  = MporeSalt,precursor
Mwater,ALE = MporeWater,precursor
Vpore,ALE  = VporeFilled,precursor
Cinit,ALE  = MporeSalt/(MporeSalt+MporeWater)
```

激活不再修改任何土壤 DOF。事务状态只在整个 DuMuX/ALE 步通过后提交，拒步不会重复转移。
