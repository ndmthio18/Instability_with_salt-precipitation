# 修正增长律与 ALE 的同口径比较

## Target

说明修正增长律与 ALE 应比较哪个高度，并推导核心固相分数
\(\phi_{\mathrm{core}}\) 的来源。

## Status

**COHERENT AFTER REFRAMING / EXTRA ASSUMPTION。**

把高度限定为“晶体库存等效净高度”后，两种模型总量一致；固定底面 ALE
几何顶高不是同一个物理量。当前结果是守恒一致性检查，不是独立动力学验证。

## Invariant Object

单位面积盐壳中的晶体质量库存 \(M_{\mathrm{cr}}\)。

## Assumptions

- 界面液相取 \(C=C_{\mathrm{eq}}\)。
- 壳内浓度扰动小，且 \(Pe_\ell\ll1\)。
- 盐壳内反应—扩散场准稳态。
- 厚壳极限 \(h/\ell\gg1\)，顶部和底部活性层不重叠。
- 顶部活性区形成稳定推进结构。
- 新生顶面材料孔隙率为 \(\varepsilon_t\)。
- 顶部活性区最终把进入该区的盐全部转化为晶体。
- \(j(t)\) 和 \(B_{\mathrm{soil}}(t)\) 是既有土壤计算的单向规定输入。

## Notation

- \(j\)：向上的水质量通量。
- \(B_{\mathrm{soil}}>0\)：土壤侧向下的盐扩散通量。
- \(P\)：进入盐壳并导致净固体增加的盐通量。
- \(S_{\mathrm{ref}}=C_{\mathrm{eq}}j/(1-C_{\mathrm{eq}})\)：平衡界面的参考盐通量。
- \(X=\sqrt{Da^*}\)：顶部表面反应与内部反应—扩散的分配参数。
- \(\varepsilon_t\)：新生顶面材料孔隙率。
- \(\phi_{\mathrm{core}}\)：厚壳稳定推进解预测的核区固相分数。

## Derivation Strategy

先由界面水盐守恒得到 \(P\)，再用顶部表面析出决定推进速度 \(v_t\)，
最后用全部盐转化后的晶体库存闭合求出 \(\phi_{\mathrm{core}}\)。

## Derivation Map

1. 界面水盐守恒给出净供盐 \(P\)。
2. 厚壳解把 \(S_{\mathrm{ref}}\) 分为表面析出和层内析出。
3. 表面析出与新生材料堆积密度决定 \(v_t\)。
4. 顶部区的总晶体产量决定后方核区固相分数。
5. ALE 使用总晶体库存定义同口径等效净高度。

## Main Derivation

### Step 1. 土壤反扩散后的净供盐

界面总液体质量通量为 \(m\)：

\[
P=mC_{\mathrm{eq}}-B_{\mathrm{soil}},\qquad
j=m(1-C_{\mathrm{eq}})+B_{\mathrm{soil}}.
\]

消去 \(m\)：

\[
\boxed{
P(t)=\frac{C_{\mathrm{eq}}j(t)-B_{\mathrm{soil}}(t)}
{1-C_{\mathrm{eq}}}}.
\tag{1}
\]

### Step 2. 顶部区的盐通量分配

厚壳极限下，

\[
\boxed{G_{\mathrm{top}}=\frac{X}{1+X}S_{\mathrm{ref}}},
\qquad
\boxed{R_{\mathrm{up}}=\frac{1}{1+X}S_{\mathrm{ref}}},
\qquad X=\sqrt{Da^*}.
\tag{2}
\]

因此 \(G_{\mathrm{top}}+R_{\mathrm{up}}=S_{\mathrm{ref}}\)。

### Step 3. 顶面推进速度

新生材料的固相分数为 \(1-\varepsilon_t\)。表面析出控制顶面推进：

\[
\rho_{\mathrm{cr}}(1-\varepsilon_t)v_t
=G_{\mathrm{top}}
=\frac{X}{1+X}S_{\mathrm{ref}}.
\tag{3}
\]

所以

\[
v_t=\frac{X}{1+X}
\frac{S_{\mathrm{ref}}}
{\rho_{\mathrm{cr}}(1-\varepsilon_t)}.
\tag{4}
\]

### Step 4. 核区固相分数

稳定推进时，表面析出和顶部层内析出最终都留在后方核区。因此

\[
\rho_{\mathrm{cr}}\phi_{\mathrm{core}}v_t
=S_{\mathrm{ref}}.
\tag{5}
\]

把式 (4) 代入式 (5)：

\[
\boxed{
\phi_{\mathrm{core}}
=\frac{S_{\mathrm{ref}}}{\rho_{\mathrm{cr}}v_t}
=(1-\varepsilon_t)\frac{1+X}{X}
=(1-\varepsilon_t)
\frac{1+\sqrt{Da^*}}{\sqrt{Da^*}}}.
\tag{6}
\]

这里最后两行之间是乘法，不是相加。该因子表达：顶面只直接接收
\(X/(1+X)\) 的盐，但顶部活性区整体接收全部 \(S_{\mathrm{ref}}\)。

本算例

\[
\varepsilon_t=0.43,\qquad
Da^*=3.17503557,\qquad
\sqrt{Da^*}=1.78186295,
\]

因此

\[
\boxed{\phi_{\mathrm{core}}=0.88988992},\qquad
\boxed{\varepsilon_{\mathrm{core}}=1-\phi_{\mathrm{core}}
=0.11011008}.
\tag{7}
\]

物理可行性要求 \(\phi_{\mathrm{core}}\le1\)，即

\[
\sqrt{Da^*}\ge
\frac{1-\varepsilon_t}{\varepsilon_t}.
\tag{8}
\]

本算例右端为 1.32558，小于 1.78186，因此没有违反体积分数上界。

### Step 5. 修正模型净高度

\[
\boxed{
\frac{dh_{\mathrm{net}}}{dt}
=\frac{P(t)}{\rho_{\mathrm{cr}}\phi_{\mathrm{core}}}}
\tag{9}
\]

以及

\[
h_{\mathrm{net}}^{\mathrm{corr}}(t)
=h_0+\frac{1}{\rho_{\mathrm{cr}}\phi_{\mathrm{core}}}
\int_0^tP(\tau)\,d\tau.
\tag{10}
\]

### Step 6. ALE 的同口径净高度

允许局部固相完全溶解时，以非负固相体积

\[
S_i=V_i(1-\varepsilon_i)\ge0
\]

为守恒变量；\(S_i=0\) 对应 \(\varepsilon_i=1\)，该单元不贡献晶体库存。
与式 (10) 同口径的 ALE 高度定义为

\[
\boxed{
h_{\mathrm{net}}^{\mathrm{ALE}}(t)
=h_0+
\frac{M_{\mathrm{cr}}^{\mathrm{ALE}}(t)
-M_{\mathrm{cr}}^{\mathrm{ALE}}(0)}
{\rho_{\mathrm{cr}}\phi_{\mathrm{core}}}}.
\tag{11}
\]

## Remarks and Interpretation

- \(\phi_{\mathrm{core}}\) 不是从 ALE 拟合出来的，而是厚壳稳定推进解析闭合的预测。
- 本次 ALE 末态最小孔隙率为 0.16275，对应最大局部固相分数约 0.83725，
  并未直接复现理论值 0.88989。
- 因此把 ALE 晶体库存也除以同一个 \(\phi_{\mathrm{core}}\)，只是在定义共同的
  等效高度；这不能验证 \(\phi_{\mathrm{core}}\)。
- 不使用 \(\phi_{\mathrm{core}}\) 时，可直接比较晶体质量增量：
  修正式累计供盐为 \(0.606399\,\mathrm{kg\,m^{-2}}\)，ALE 晶体增量为
  \(0.606438\,\mathrm{kg\,m^{-2}}\)。

## Boundaries and Non-Claims

- 固定底面 ALE 的几何顶高不等于式 (11) 的净库存高度。
- 若要验证 \(\phi_{\mathrm{core}}\)，应从 ALE 的稳定核区平台或实验截面独立测量
  固相分数，再与式 (6) 比较。
- 若要预测实验室坐标中的壳顶位置，还需底部退缩、壳体下沉或接触条件闭合。
- 本次 \(B_{\mathrm{soil}}\) 为单向输入，不是同步双向耦合解。

## Open Risks

- 当前 ALE 在有限计算时段内尚未形成与理论假设完全相同的均匀核区平台。
- 反应层尺度接近或小于晶粒尺度时，连续介质 \(\phi_{\mathrm{core}}\) 的解释可能失效。
- 原土壤直接 SCVF 净盐通量与式 (1) 重构供盐累计相差约 1.92%，论文中不可混用。
