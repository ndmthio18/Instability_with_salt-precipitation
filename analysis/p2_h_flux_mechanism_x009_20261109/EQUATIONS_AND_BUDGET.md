# ALE 通量与盐质量预算

## 符号与方向

面编号从盐壳底部 `0` 到顶部 `n`，正通量方向为向上。单元 `i` 的孔液溶解盐库存（单位土壤面积）为

\[
Q_i=\rho\,\varepsilon_i C_i h\Delta\xi_i .
\]

生产算子的离散守恒式为

\[
\dot Q_i=J_{s,i}-J_{s,i+1}-R_i h\Delta\xi_i .
\]

其中 `R>0` 为内部沉淀，`R<0` 为内部溶解。

## ALE 内部面

移动网格体积通量为

\[
J_{mesh,i}=-\xi_i\dot h\,\varepsilon_{upwind,i},
\]

液相相对于移动网格的质量通量为

\[
J_{rel,i}=J_{liquid,i}+\rho J_{mesh,i}.
\]

因此内部面输盐通量是

\[
J_{s,i}=J_{adv,i}+J_{diff,i}
=J_{rel,i}C_{upwind,i}
-\rho\frac{C_i-C_{i-1}}{
h[(\xi_i-\xi_{c,i-1})/(\varepsilon_{i-1}D_{eff,i-1})
+(\xi_{c,i}-\xi_i)/(\varepsilon_iD_{eff,i})]} .
\]

这里 `D_eff=sqrt(epsilon) D_s`，所以在几何条件相同的近似下，扩散传输系数含有约 `epsilon^(3/2)` 的依赖。报告中的梯度比、精确离散电导比和 `epsilon^(3/2)` 比分别给出，避免仅凭通量衰减把原因归于反应。

## 底部面与顶部面

底部面是液相输盐：

\[
J_{s,0}=C_bJ_{liquid,0}
-\rho\varepsilon_0D_{eff,0}\frac{C_0-C_b}{h\xi_{c,0}}.
\]

顶部面不是普通液相面。生产代码中的

\[
J_{s,n}=G_{top}=\rho_{cr}(1-\varepsilon_{birth})\dot h
\]

表示新生顶部固体通量。因此 `FACE_FLUX_PROFILES.csv` 将顶部标为 `top_solid_birth`，且顶部 `J_adv`、`J_diff` 留空，不把 `G_top` 错标为液相通量。

## 接受步积分守恒

本分支把所有面通量和所有单元反应作为附加保守状态，与 `E=epsilon*V`、`Q=rho*C*E` 由同一个 CVODE BDF 历史积分。对接受步 `[t_n,t_{n+1}]`：

\[
Q_i^{n+1}-Q_i^n=
\int J_{s,i}dt-\int J_{s,i+1}dt-\int R_iV_i dt .
\]

从底部累加至任意单元 `k`：

\[
Q_{0:k}^{n+1}-Q_{0:k}^n=
\int J_{s,0}dt-\int J_{s,k+1}dt-\sum_{i=0}^{k}\int R_iV_i dt .
\]

这与“用步末瞬时通量乘以 `dt`”不同。`LOCAL_STORAGE_BALANCE.csv` 保存所有步的整壳残差、每步最差局部残差，以及选定时刻的逐单元和逐前缀累积残差。

整壳孔液、固体和总盐满足

\[
\dot Q_{pore}=J_{s,0}-G_{top}-R_{net},\qquad
\dot M_{solid}=G_{top}+R_{net},
\]

\[
\frac{d}{dt}(Q_{pore}+M_{solid})=J_{s,0}.
\]

## 机制诊断尺度

局部反应–扩散尺度定义为

\[
\ell(z,t)=\sqrt{\frac{\varepsilon D_{eff}}
{(1-C_{eq})a_vk_r}}.
\]

`h/ell` 仅作为当前非均匀数值解的诊断量，不是自动成立的解析解或新的物理闭合。

参考比值 `G_top/[C_eq j_evap/(1-C_eq)]` 的分母仅是饱和盐水蒸发参考通量，不被解释为普适生长上限。
