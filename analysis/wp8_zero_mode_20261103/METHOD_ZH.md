# WP8 方法与解释边界

本轮只重用 Case 36877 已有数值解。canonical nx=28 是原完整域物理波数协议；R24/R48 是一个物理波长内的横向网格单元数，不是将物理波数改成24/48。W=0.6/28 m，H=0.2 m，垂向200个非均匀单元；初始M1幅值保持1e-6，M0幅值0。T=136423.0637723484 s；4T约6.316日。

## 平均、扰动及三个不同参照

采用BOX顶点双控制体宽度，端点半权：xbar=sum(wx*x)/W；每条解的扰动为x'=x−xbar。A=M1−M0的水平平均；B=M0−xb；C=M1−xb。A仅表示初始播种敏感性，不能在没有长期xb时当作对流相对于无对流背景的因果效应。B/C使用原一格挤出背景，仍是未通过严格横向均匀性门槛的候选参照。

零均值扰动使用与无通量侧边相容的Neumann余弦基 cos(m*pi*x/W)，m=1..N；canonical m=2，k=293.215314335 rad/m。采用双体积权重的离散DCT-I投影及Nyquist端点归一化，E_m=∫z <x'_m²>_x dz。Parseval已逐帧核查。垂向RMS幅值=sqrt(E_m/H)；表面有符号余弦系数另存，两者不能混用。

k_peak对应E_m最大项的实际k；k_E=sum(k_m E_m)/sum(E_m)。总非零能量<=1e-26 (mol/mol)^2 m时不报告波数；这是浮点噪声显示约定，不是物理检测阈值。均值差的能量E0=∫(Delta xbar)^2 dz单独存储，绝不塞入非零能量分母。只有以同一背景为参照时，E0+sum E_m才是该参照下的浓度差平方范数；A下也可作与水平平均M0的数学分解，但不是相对完整二维M0场的差能量。

## 守恒量与符号

qbar=phi*<Sw*x>=phi*(Swbar*xbar+cov(Sw,x))，协方差由同一帧同一权重直接求出，所有计算使用乘积的真实平均。库存I=∫区域phi*Sw*x dV/W，单位m，是归一化守恒量，不是kg盐。累计C=∫J_down dt，单位m；J=J_adv+J_diff，单位m/s，向下离开上部区域为正。现有顶部盐通量为零。采用求解器接受步的右端通量积分，沿用SCVF预算定义。

Delta I_c=(I_case(t)−I_ref(t))−(I_case(0)−I_ref(0))，Delta C=C_case−C_ref，Delta R=R_case−R_ref，核查Delta I_c+Delta C=Delta R。初始与当时库存两个归一化分母都明确输出；不以小差值作分母报告夸大的百分误差。机器恒等式闭合不等于连续模型误差小，也不证明SCVF核的物理真实性。

仅用既有same_time(32 ULP)契约匹配实际接受时刻，不做累计通量插值，不拼接外推背景。浮点接近而超过契约的终点也丢弃。

## 文献动机

[Tilton & Riaz (2014), JFM 745, 251–278, DOI 10.1017/jfm.2014.72](https://doi.org/10.1017/jfm.2014.72)讨论瞬态扩散边界层中单一横向模式的非线性相互作用产生零波数响应。本轮借用“将均值响应与非零模式分开”的分析思路；其重力不稳定扩散问题不能直接证明本项目蒸发、非饱和、密度耦合模型的定量非线性反馈。没有把该文结果用作本项目验证数据。

## 复现

已存在jobs.json时禁止再次launch。原场重算需要原始归档；已完成的sal只读产品可运行`python scripts/wp8_campaign.py collect`按SHA恢复忽略目录缓存。随后运行`analyze_wp8.py`、`wp8_observables.py`、`plot_wp8.py`、`report_wp8.py`、`test_wp8.py`及`verify_wp8.py`。源码与缓存SHA保存在manifest/audit；不提交大型原场和NPZ。
