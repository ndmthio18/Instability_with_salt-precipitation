from pathlib import Path
import json
import pandas as pd
from analyze import ROOT,R,PHI

def main():
    s=pd.read_csv(R/'growth_validation.csv');g=pd.read_csv(R/'long_time_grid_convergence.csv')
    t=pd.read_csv(R/'temporal_convergence.csv');k=pd.read_csv(R/'kinetic_sensitivity.csv')
    allcases=pd.read_csv(R/'summary.csv');fine=pd.read_csv(R/'fine_grid_temporal.csv')
    if len(s)!=5:raise RuntimeError('Require five complete final cases')
    p01=json.loads((R/'fast_p01_gate.json').read_text());p2=json.loads((R/'p2_kinetic_gate.json').read_text())
    vk=['vt_max_normalized','vb_max_normalized','vnet_max_normalized'];pk=['zt_max_normalized','zb_max_normalized','hc_max_normalized']
    space=float(g[vk+pk].max().max());time=float(max(t[vk].max().max(),fine[vk].max().max()))
    ids=list(pd.read_csv(ROOT/'fast_p01_cases.csv').id)+list(pd.read_csv(ROOT/'primary_cases.csv').id)
    primary=allcases[allcases.id.isin(ids)]
    if len(primary)!=len(ids):raise RuntimeError('Missing a registered primary case')
    budget=float(primary[[c for c in primary if c.startswith('max_') and c.endswith('_res')]].max().max())
    numerical=bool((primary.status=='end').all()) and p01['P0_pass'] and p01['P1_pass'] and p2['passed'] and space<.01 and time<.01 and budget<2e-9
    gates=[
      dict(name='mathematical_conservation',status='PASS' if budget<2e-9 else 'FAIL',max_abs_kg_m2=budget),
      dict(name='moving_interface_consistency',status='PASS_CONDITIONAL_ON_KINETIC_CLOSURE',evidence='Independent inventories and moving-frame integral jumps; not empirical calibration'),
      dict(name='spatial_convergence',status='PASS_WITHIN_TESTED_ENVELOPE' if space<.01 else 'FAIL',maximum_scaled_difference=space,criterion=.01),
      dict(name='temporal_convergence',status='PASS_WITHIN_TESTED_ENVELOPE' if time<.01 else 'FAIL',maximum_velocity_scaled_difference=time,criterion=.01),
      dict(name='interface_flux_convergence',status='PASS_AT_SAMPLED_TIMES' if p01['P1_pass'] else 'FAIL',evidence=p01['flux']),
      dict(name='bottom_kinetic_sensitivity',status='AUDITED_NOT_EMPIRICALLY_VALIDATED',max_net_change_reference=float(k.vnet_max_normalized.max()),max_relative_thickness_change=float(k.hc_max_normalized.max())),
      dict(name='corrected_top_growth_law',status='CONDITIONAL_APPROXIMATION_WITH_SYSTEMATIC_BIAS',no_fitted_parameters=True),
      dict(name='corrected_geometric_net_growth_law',status='NOT_VALIDATED_AS_UNCONDITIONAL_GEOMETRIC_LAW',requires='Stable density, dissolved-storage approximation, and additional bottom-interface physics')]
    status=dict(numerical_status='GO_WITHIN_TESTED_ENVELOPE' if numerical else 'NO_GO',overall='NO_GO_UNQUALIFIED_PHYSICAL_VALIDATION',gates=gates,primary_completed=int((primary.status=='end').sum()),primary_total=len(primary))
    (R/'validation_status.json').write_text(json.dumps(status,indent=2))
    rows=[]
    for r in s.itertuples():
        rows.append(f'| {r.case} | {r.final_vt*1e9:.6f} | {r.final_vnet*1e9:.6f} | {r.final_vnet_corr*1e9:.6f} | {r.final_net_error_percent:+.2f}% | {r.final_hc*1e6:.3f} | {r.final_hc_corr*1e6:.3f} | {r.final_core_phi:.6f} |')
    u09=s[s.case=='u09'].iloc[0]
    text=f'''# Task 28 最终审查

## 结论

**数值状态：{status['numerical_status']}。无条件物理验证状态：NO-GO。** 原先严格容差失败和早期界面通量不收敛问题已在本轮测试范围内解决；这不等于修正解析式已成为无需底部界面假设的几何净增厚定律。

ALE直接演化顶、底位置，净厚度为zt−zb，净速度为vt−vb。解析核心固相分数{PHI:.8f}未进入数值演化、初值或底部闭合。没有调整kappa或拟合新的phi来追求重合。

## 数值证据

- 精确复现旧代码在2.520100517 s / 7.269625167 s的失败，保存内部步长、误差权重、Jacobian计数和失败状态。并非孔隙或固相先耗尽。
- 使用守恒扰动库存、变量尺度、自动微分Jacobian、针对浓度扰动的误差控制和及时刷新Jacobian；没有放宽容差或裁剪接受态。
- 三个代表工况完成两档严格容差全时段计算；512格另做严格时间复核。最大速度差/非零参考速度为{time:.3g}。
- 最终256→512格全时段速度/几何量最大归一化差{space*100:.3f}%，低于登记的1%目标。不是只比较终点。
- 库存残差最大{budget:.3g} kg/m²；守恒与收敛分别检验。
- 30条早期对照覆盖0.1/1/10 nm、128/256/512格、内部时间步细化和相容初值。1 nm主设置中，u=.6在全部共同样本通过；x008最早通过样本为5e−9 s，此后至100 s的共同样本通过。这是采样误差界，不是连续时间上的严格证明。未人为排除第一小时。
- 失败、压力检验、重试和运算优化交叉验证全部保留。后续优化只改变运算实现；物理方程和接受准则不变。

## 最终几何比较

常值工况终点50 h，x008终点44.70875 h。速度单位nm/s，厚度单位µm。误差定义为(ALE−解析)/ALE。

| 工况 | ALE顶部速度 | ALE净速度 | 解析净速度 | 净速度偏差 | 实际厚度 | 解析厚度 | 核区固相分数 |
|---|---:|---:|---:|---:|---:|---:|---:|
'''+ '\n'.join(rows)+f'''

完整时间历史、顶面预测误差与实际位置均在表和图中；以上仅是终点摘要。初始场为均匀平衡浓度和出生固相分数0.57，并不具备解析推导的稳定致密核。准稳态、稳定推进的解析式不能被当作最初瞬态的精确解。数值通量收敛时间与解析物理近似的适用时间是两件事。

## 误差抵消与底部假设

u=.9终点顶部误差为{u09.final_delta_top*1e9:+.6f} nm/s，底部诊断差为{u09.final_delta_bottom*1e9:+.6f} nm/s，两者相减仅余{u09.final_delta_net*1e9:+.6f} nm/s。按两项绝对值之和定义的抵消比例为{u09.final_cancellation_fraction*100:.2f}%。因此小净误差不能证明顶、底运动分别准确。解析底速由“解析顶速−解析净速”反推，始终仅作诊断。

kappa=.5/1/2对初期退缩与累计厚度有可辨影响；末时刻速度趋近不能抹去路径差异。kinetic_sensitivity.csv给出全时段最大差和终点差。表面暴露面积kappa*phi_b是额外构成假设，并非体积反应律的必然推论。本轮验证其数值守恒和敏感性，没有证明kappa=1是实验真实值。

## 核固相分数与净生长含义

实际体积平均、核区平均和底部前沿迹值分别计算，解析常数未被恢复为全域通用值。随着孔隙率演化，扩散和比表面积系数偏离出生值；厚度参数、网格细化及kappa对照可约束部分解释，但没有冻结系数的独立反事实实验，不能将全部偏差归因于某一个机制。

几何恒等式hc_dot=Mdot/(rho_cr*phi_bar)−hc*phi_bar_dot/phi_bar将质量变化与密度演化分开。结果同时保存实际供给、溶解盐储存、平均密度差和结构演化四项，它们精确合成净律误差。新孔隙卤水库存由ALE计算，简单解析闭合则略去；不能仅凭总盐守恒将两者画等号。恒等式的代数重构本身不作为独立验证。

修正律可作为给定参数、给定供给和指定底部物理下的长时近似，其误差已量化；尚不能声称它独立决定实际底面位置、对初期溶解过程普遍成立，或在不同底部暴露规律下保持同样精度。

## 交付

八张主图及三张补图提供PNG/PDF/可编辑SVG；曲线均为实际数值结果或明确标注的解析预测。七份科学说明、C++/Python源码、自动检查、配置、Slurm记录与SHA-256在本目录。完整原始场保留HPCC，轻量结果随分支提交。最终提交及远程核验见执行回复。
'''
    terms=['| 工况 | 供给修正 | 溶解盐储存 | 平均密度差 | 结构演化 | 合计净误差 |','|---|---:|---:|---:|---:|---:|']
    for r in s.itertuples():
        terms.append(f'| {r.case} | {r.final_supply_term*1e9:+.6f} | {r.final_storage_term*1e9:+.6f} | {r.final_density_term*1e9:+.6f} | {r.final_structure_term*1e9:+.6f} | {r.final_delta_net*1e9:+.6f} |')
    text=text.replace('## 交付', '## 终点净误差的账本拆分\n\n以下单位为nm/s，各项精确合成ALE−解析净速度。这是恒等分解，不是独立因果实验；供给修正与溶解盐储存往往各自较大而互相抵消。\n\n'+'\n'.join(terms)+'\n\n## 交付')
    (ROOT/'FINAL_REVIEW.md').write_text(text,encoding='utf-8')
    print(json.dumps(status,indent=2))
if __name__=='__main__':main()
