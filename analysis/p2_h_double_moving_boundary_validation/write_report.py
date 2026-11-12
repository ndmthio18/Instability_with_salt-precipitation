"""Write the task-27 evidence-bounded report from computed result tables."""
from pathlib import Path
import json
import pandas as pd
from analyze import CR, PHI

R=Path(__file__).resolve().parent

def main():
    s=pd.read_csv(R/'results/summary.csv')
    h=pd.read_csv(R/'results/history.csv.gz')
    cv=pd.read_csv(R/'results/convergence.csv')
    sensitivity=pd.read_csv(R/'results/sensitivity.csv')
    meta=json.loads((R/'results/summary.json').read_text())
    base=s[(s.kappa==1)&(s.gap==1e-9)&(s.rtol==1e-8)&s.id.str.contains('_n')].sort_values('nc').groupby('u').tail(1).sort_values('u')
    def label(u):return 'x008' if u<0 else f'u={u:g}'
    table=['| 工况（最细网格） | 终点 h | 顶部 μm | 底部 μm | 实际厚度 μm | ALE 净速度 nm/s | 修正律 nm/s | 净速度相对差 |',
           '|---|---:|---:|---:|---:|---:|---:|---:|']
    core=['| 工况 | 顶部速度相对差 | 实际平均固相分数 | 核平均固相分数 | 核空间标准差 | 最后20%时间核均值相对漂移 |',
          '|---|---:|---:|---:|---:|---:|']
    for a in base.itertuples():
        table.append(f'| {label(a.u)} ({a.nc}) | {a.end_s/3600:.5f} | {a.final_zt*1e6:.3f} | {a.final_zb*1e6:.3f} | {a.final_hc*1e6:.3f} | {a.final_vnet*1e9:.6f} | {a.final_vnet_corr*1e9:.6f} | {a.final_net_error_relative:.2%} |')
        core.append(f'| {label(a.u)} | {a.final_top_error:.2%} | {a.final_phi_mean:.6f} | {a.final_core_phi:.6f} | {a.final_core_std:.6f} | {a.late_core_relative_drift:.3%} |')
    a=base[base.u==.6].iloc[0];d=h[h.id==a.id].iloc[-1]
    grid=['| 网格对 | 顶部速度差 / vref | 底部速度差 / vref | 净速度差 / vref | 净速度相对差（避开零值） | 厚度相对差 | 晶体库存相对差 | 溶解盐库存相对差 |',
          '|---|---:|---:|---:|---:|---:|---:|---:|']
    for row in cv.itertuples():
        grid.append(f'| {row.id} → {row.reference_id} | {row.vt_scaled_max:.2%} | {row.vb_scaled_max:.2%} | {row.vnet_scaled_max:.2%} | {row.vnet_relative_max_away_zero:.2%} | {row.hc_scaled_max:.2%} | {row.solid_scaled_max:.3%} | {row.dissolved_scaled_max:.3%} |')
    sens=['| 对照 | 共同比较时长 h | 净速度差 / vref（最大） | 厚度相对差（最大） |','|---|---:|---:|---:|']
    for row in sensitivity.itertuples():
        sens.append(f'| {row.id} | {row.common_end_s/3600:.6f} | {row.vnet_scaled_max:.4%} | {row.hc_scaled_max:.4%} |')
    text=f'''# Task 27：双移动边界 ALE 的几何净增厚验证

## 结论与适用范围

**科学状态：守恒和基本物理验证通过；长时几何计算已获得，但严格时间容差下的u=0.6数值鲁棒性未解决。因此修正净增长律的最终定量验证为NO-GO。**

已完成双区域、双移动边界的独立守恒计算。共登记 {len(s)} 个正式/加密/容差算例，其中 {(s.status=="end").sum()} 个达到其输入支持终点；失败记录保留。所有基准与空间加密算例完成，恒定供给为50 h，x008为44.70875 h，未外推土壤输入。x008时间原点为土壤首次析盐19048.5 s，其终点对应原土柱模拟50 h。
顶面、底面分别由局部动力学和移动界面守恒计算，主观测量为 **hc=zt-zb**。
解析增长律和理论 phi_core 从未用于推进界面，全部晶体库存换算仅用于后处理对照。

这一研究得到的是**附加底部界面动力学闭合下的条件性验证**。原有内部面积律 av∝phi 使局部晶体量在有界欠饱和度下指数衰减，不能唯一给出有限时间耗尽前沿。因此选用任务允许的 Option 2：暴露面积为 kappa_b*phi_b、界面采用原有本征 kr。kappa_b=1 是投影面积假设，0.5和2用于敏感性，未对增长律拟合。

**H1 与 H2 必须分开解释。** 本次 ALE 保留动态孔隙率、动态扩散/反应面积、溶解盐储存及液层输运；解析式使用常系数反应层和近稳态近似。两者同用基准材料参数，不等于二者的所有近似完全相同。观察到的偏差衡量解析式在此演化问题中的近似能力，不能据此声称其在严格假设下的数学推导错误。

## 独立几何结果

{chr(10).join(table)}

速度相对差以实际ALE末态净速度为分母；不是以两个大速度之差接近零时的百分数来掩盖绝对误差。完整的参考速度归一化误差也见 results/history.csv.gz。
上述位置是静止骨架假设下的实验室坐标；核未作刚性下沉处理。初始底面为1 nm、初始壳厚3.5 μm。

在 u=0.6 的最细网格，50 h 几何净速度为 {d.vnet*1e9:.6f} nm/s，解析预测为 {d.vnet_corr*1e9:.6f} nm/s；实际壳厚 {d.hc*1e6:.3f} μm，积分解析预测 {d.h_corr*1e6:.3f} μm。
这不是“以同一个 phi_core 换算两遍库存”的结果。库存等效厚度为 {d.h_inventory*1e6:.3f} μm，定义为 h0+[Mcr(t)-Mcr(0)]/(rho_cr*phi_core)，以相同初始壳厚对齐，作为单独曲线展示。

## H1：顶部增长律

Da 使用实际活动壳厚 hc，sqrt(Da)=hc/ell；ell={meta['ell_m']:.9g} m，Da*={meta['Dastar']:.9g}。
主预测使用题目规定的 Bsoil；另以计算得到的壳侧扩散通量 Fplus 代入，只作条件性边界诊断。
液层与界面反应存在时 Fplus一般不等于Bsoil，不能把该诊断称为独立预测。
厚壳末期底部扰动受到屏蔽，两个版本接近，剩余偏差还反映活性层系数随孔隙率演化，不能全归于土壤反扩散。

{chr(10).join(core)}

H1表现为有界近似，而非本次长时问题中的精确等式。短时固定底部恢复验证与长时解析近似误差是两种不同检验；前者通过不代表后者自动通过。

## H2：为何库存关系不能自动等同于几何净增厚

理论 phi_core={PHI:.9f}。表中独立计算的核平均值和全壳平均值均应与它比较。
核统计取 [zb+3ell, zt-3ell]，按单元重叠体积加权；没有这个区间的早期时刻记为缺失，而非填入理论值。
晚期核可以出现较平坦区域，但这并不证明其数值等于理论 phi_core，也不保证边界层库存可以忽略。

精确关系是：

    Mcr = rho_cr * phi_mean * hc
    hc_dot = Mcr_dot/(rho_cr*phi_mean) - hc*phi_mean_dot/phi_mean
    Mcr_dot = Gtop - Db + integral(R dz) = Jsoil - Mdiss_dot

所以要把 P/(rho_cr*phi_core)解释成实际几何净增厚，需要实证支持平均固相分数近似恒定且取该理论值，以及 Jsoil-Mdiss_dot≈P；并不要求Jsoil和Mdiss_dot分别等于P和零。

u=0.6 最细网格末态的独立分解为：

- 实际晶体质量增长率 / P = {d.solid_gain_rate/d.P_corr:.8f}。
- 实际全壳平均固相分数 = {d.phi_mean:.8f}，理论值 = {PHI:.8f}。
- 几何净速度减去解析值 = {(d.vnet-d.vnet_corr)*1e9:.8f} nm/s。
- 平均固相分数不同的贡献 = {d.mean_density_term*1e9:.8f} nm/s。
- 结构继续演化的贡献 = {d.geometry_structure_term*1e9:.8f} nm/s。
- 溶解盐储存贡献 = {d.storage_term*1e9:.8f} nm/s；实际边界盐供给与P的差贡献 = {d.supply_term*1e9:.8f} nm/s。

这些项按精确恒等式求得，其总和另行核验。供给与储存贡献会显著抵消，不能仅选其中一项解释偏差。
个别工况净速度较接近，并不能独立证实 phi_core 闭合；不同贡献也可能抵消。**H2不能无条件升级为几何增长定律。**

## 空间与时间分辨率

三档基准网格为64/128/256个壳单元，液层为其1/4。因为u=0.6的公式偏差与残余网格差仍需分离，且u=0.9的公式偏差较小，追加两者的512/128网格，未改变方程、参数或时间容差。
下表比较相邻网格在共同区间 t≥1 s 的最大差；vref=Sref/(rho_cr*phi_core)，位置差除以细网格hc，库存差除以细网格库存。
这不是Richardson外推误差，也不是“最细格就等于真值”。完整表包含顶/底位置、局部界面残差及各区域守恒残差。

{chr(10).join(grid)}

细化差应与公式偏差一起引用。特别是净速度由两个速度相减，不能只凭总晶体质量收敛就宣布几何结果完全收敛。

时间容差1e-8→1e-9并将最大步长30→15 s；其中u0.6_strict在内部t=2.5201 s发生CVODE重复误差检验失败，最后输出2 s。追加同一可执行程序、rtol=1e-10、最大步长15 s的u0.6_strict_tighter仍在内部t=7.26963 s失败，最后输出5 s。两条失败均保留，未再放宽阈值。**u=0.6的50 h时间容差收敛尚未建立，不能据短前缀的小差异声称通过；以下偏差及机制分解属于条件性计算证据，尚不足以完成定量验证。**液层、初始种子厚度和底部面积闭合的对照如下。
所有对照以同u的128/32基准为参照，kappa变化是物理闭合敏感性，不是数值误差。

{chr(10).join(sens)}

## 瞬时通量的独立收敛限制

早期壳侧扩散通量有明显数值敏感尖峰。x008在1 s的Fplus：128格为+1.58146e-5、256格为-1.05635e-4 kg/m²/s，而128格严格容差为+2.35416e-6 kg/m²/s。符号和幅度均不稳定，不能解释成已经证实的物理反向扩散事件，也不能仅凭守恒残差很小就认为该通量可信。

在x008的128格时间容差对照中，1 h以后Fplus最大差为1.215e-10 kg/m²/s；128→256格的对应最大差为7.543e-9 kg/m²/s。1 h是明确展示的诊断窗口，不是人为修改求解方程或删除早期输出。所有时刻仍保留在数据表中。Figure 8单独展示0–1 h未收敛通量及1 h以后通量，避免初始尖峰压扁后期曲线；results/flux_resolution.csv给出1 s、60 s和1 h起算的对照。

因此，早期尖峰的具体成因尚需数值分析；不能将所有早期波动都归因于准稳态近似失效。底部界面残差是代数匹配检验，瞬时通量精度需要另外的网格和容差证据。

## 守恒与历史失败

全体正式轨迹的全局/分区质量残差最大为 {meta['budget_max']:.6g} kg/m²，低于登记门槛2e-9 kg/m²；界面代数残差另列，未把代数恒等式当作独立积分守恒证明。
物理前置验证、失败版本和重映射测试见 CONSERVATION_AND_REMAP_TESTS.md；原始失败日志与源代码快照保留。
所有生产轨迹均未发生exact-zero前缀重映射事件（次数见summary.csv）；真正的底部运动由连续锐界面动力学产生。
守恒转移已用制造的精确零固相数据验证，不能把这个测试说成真实轨迹中已经观测到完整耗尽事件。

## 对任务十个问题的直接回答

1. **能否恢复原固定底部ALE？** 能，最大顶部速度相对差2.40e-7，底面位移仅约2.8e-33 m舍入量级。
2. **能否独立确定底部速度？** 在明确增加暴露面积和局部界面动力学后能；仅靠原守恒方程不能唯一确定。
3. **能否守恒处理零固相转移？** 精确连通前缀转移测试通过；不删除正晶体量，不把孤立内部液体单元当作整个底部已移动。当前生产没有触发该转移。
4. **能否超过原epsilon=0.999终止时间？** 新闭合的u>0基准均算至50 h。原128格u=0.3/0.6/0.9分别在2640/1320/900 s停于该筛查，但两组底部物理闭合不同，不能宣称只是修复了同一模型的数值终止。
5. **顶部是否持续生长而底部同时退缩？** 是，u>0和x008的计算轨迹直接展示；u=0仅有浮点量级的底部移动。
6. **顶部修正律是否仍准确？** 见H1误差表；是受近似条件约束的预测，长时存在动态系数及网格误差，不是精确恢复关系。
7. **净增长律是否预测真实几何净增厚？** 需要附加结构条件，本次不能给出无条件肯定。库存增长较准不足以保证几何高度同样准确。
8. **何时、为何偏离？** 初始条件违背准稳态且早期通量尖峰数值上未收敛，不能把全部早期波动当成物理信号；长时固相分布、活性层系数及供给/储存偏离常系数闭合，精确分解和网格对照用于量化，而不是仅凭曲线猜测原因。
9. **理论核固相分数是否被独立支持？** 计算可形成近平台，但其数值与理论0.88989有系统差别；局部核存在和取理论数值是两项不同命题。
10. **真实土壤支持盐壳还缺什么？** 底部暴露面积与界面动力学的实验约束、骨架应力/位移及支撑机制、连续供液的压力/毛细/渗透率闭合，以及需要时的液区成核和底部再沉淀。这里没有求解这些机制。

## 可复现材料

- baseline c7da956041b7e42f85278603a068c3915394ec20；独立分支 research/p2-h-double-moving-boundary-validation。
- HPCC生产数组16572089、追加512格数组16572116；全部编译和求解使用cnmix计算节点。
- parameters.json、cases.csv、cases_refine.csv、forcing.csv、execution_record.json、provenance/下的源码/原始结果哈希及Slurm账目。
- results/history.csv.gz和profiles.csv.gz保留全部正式轨迹与剖面；summary、convergence、sensitivity为汇总CSV/JSON。
- figures/figure1至figure8提供PNG、PDF和可编辑文字SVG，Times New Roman字体；所有曲线来自实际算得数据。
- `python analyze.py`、`python audit_fluxes.py`、`python write_report.py`、`python audit_delivery.py`；原始数据本地raw/与远程目录保留，原有验证模型和结果未修改。
'''
    (R/'VALIDATION_REPORT.md').write_text(text,encoding='utf-8')

if __name__=='__main__':main()
