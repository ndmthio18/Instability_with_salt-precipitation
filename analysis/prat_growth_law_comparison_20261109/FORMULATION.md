# M1/M2/M3 诊断模型定义

本轮只替换盐层顶部增长效率，不改变 DuMuX 土壤方程、气膜模型、表面析盐激活、库存账本、回滚或有效域筛选。

共同参考速度为

\[
v_{\rm ref}=\frac{C_{\rm eq}j_{\rm evap}}
{\rho_{\rm cr}(1-\epsilon)(1-C_{\rm eq})}.
\]

三个候选增长律分别为

\[
v_1=v_{\rm ref}\frac{\mathrm{Da}_{\rm top}}{1+\mathrm{Da}_{\rm top}},
\]

\[
v_2=v_{\rm ref}\frac{\sqrt{\mathrm{Da}^{*}}\tanh(\sqrt{\mathrm{Da}})}
{1+\sqrt{\mathrm{Da}^{*}}\tanh(\sqrt{\mathrm{Da}})},
\]

\[
v_3=v_{\rm ref}\frac{\sqrt{\mathrm{Da}^{*}}}{1+\sqrt{\mathrm{Da}^{*}}}.
\]

其中

\[
\mathrm{Da}_{\rm top}=\frac{\beta(1-\epsilon)(1-C_{\rm eq})hk_r}{\epsilon D_{\rm eff}},
\quad
\mathrm{Da}=\frac{h^2(1-C_{\rm eq})a_vk_r}{\epsilon D_{\rm eff}},
\]

\[
\mathrm{Da}^{*}=\frac{\beta^2(1-\epsilon)^2(1-C_{\rm eq})k_r}
{\epsilon D_{\rm eff}a_v},
\quad
a_v=\frac{36(1-\epsilon)}{d_{\rm sub}},
\quad D_{\rm eff}=\sqrt{\epsilon}D_s.
\]

每个接受步通过二分法联立

\[
C_{\rm top}=C_{\rm eq}+\frac{\rho_{\rm cr}v}{\beta\rho_l k_r},
\qquad j_{\rm evap}=\mathrm{GasFilm}(C_{\rm top}),
\]

而不是给 M2/M3 复用 M1 的蒸发通量。

## 质量账本

三种增长律都暂时使用相同的新生体积账本：

\[
\Delta M_{\rm solid}=\rho_{\rm cr}(1-\epsilon)\Delta h,
\]
\[
\Delta M_{\rm pore,salt}=\rho_l\epsilon C_{\rm eq}\Delta h,
\qquad
\Delta M_{\rm pore,water}=\rho_l\epsilon(1-C_{\rm eq})\Delta h.
\]

因此 M2/M3 是“增长公式替换诊断”，不是 Prat 内部反应—扩散质量闭合。通过全局质量守恒只能验证账本实现，不能验证此物理替代。

## 有效域

- 土壤界面与盐层顶端质量分数偏离 `Ceq` 不得超过相对 2%。
- 土壤界面 `Sw >= 0.99`。
- 既有 Peclet 数限制保持不变。
- M3 另逐步报告 `Da >= 10`；低于该阈值时不得把 Eq. (29) 当作受支持的大 Da 近似。
- 任一越界试算均拒绝并回滚；不裁剪浓度、速度或盐库存。
