"""Render the measured historical-analysis findings and their limits as Chinese Markdown."""
import csv,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'analysis/dense_history_20261102'
def rows(n):
    with (OUT/n).open(encoding='utf-8',newline='') as f:return list(csv.DictReader(f))
runs=rows('run_assessment.csv');inv=rows('inventory_summary.csv');timing=rows('time_reconciliation.csv')
summary=json.loads((OUT/'summary.json').read_text());audits=json.loads((OUT/'evidence/summary.json').read_text())
transfer=json.loads((OUT/'transfer_check.json').read_text())
iv={(r['id'],r['region']):r for r in inv}
consistent=[r for r in runs if r['status'].startswith('field_input_consistent')]
assert len(consistent)==14,'Quarantined runs require an individually scoped report'
all_monotone=all(int(r['intervals_decreasing_over_1e_minus8_initial'])==0 for r in inv)
text=['# 14 条历史密集时序：盐库存与再分配分析',
      '', '日期：2026-11-02。范围：已有 14 条运行、10 个物理案例、1,414 个 VTU 场；仅后处理，未启动任何新模拟或 Pilot A/B。',
      '', '时间定义遵循 [统一命名约定](../onset_semantics_20261102/CONTRACT_ZH.md)。`t_crit` 是全局临界模态首次 crossing；`t_seed(0.10)` 是经过增长筛选（含 fallback）的 observable 模态自身 crossing，既不是 DuMuX 检测时间，也不是达到指定扰动幅值的时间。`t_detect` 是数值检测时间，`t_eff` 仍未定义。', '', '## 结论', '',
      ('**这些记录可以分析，且已完成全部时序分析。在首次观测越界前，14 条运行的三个上部控制体均没有超过描述性筛选下限的相邻帧库存下降。已有结果的主要表现为表层盐分累积，尚不能据此证明不稳定性产生了显著的额外向下再分配。**' if all_monotone else '**部分控制体出现库存回落，须结合下表、来源与窗口逐例检查；没有背景对照，仍不能直接归因为不稳定性。**'),
      '', '这里的“没有下降”只指保存帧之间，筛选下限为初始库存的 1e-8，不是经过网格/时间步验证的物理显著性阈值。净库存增加不排除“相对于无扰动背景的累积减缓”，也不排除同时存在向下扩散和向上平流。现有数据缺少匹配背景，不能识别这部分额外输运。',
      '', f"14 条运行的输入、profile、PVD、每步诊断哈希与冻结审计一致；所有场的密度和饱和度与明确输入相符。历史可执行文件身份仍未建立。{summary['observed_threshold_runs']} 条运行在保存场中达到 x=0.1，越界帧及其后数据已从主汇总中剔除。",
      '', '## 数据可用性和复核', '',
      f"- 读取原始 VTU 总计 {sum(a['bytes_read'] for a in audits)/1e9:.3f} GB，逐帧 SHA-256 留存在 `evidence/runs/*/field_manifest.json`；没有回传整套原始场。",
      f"- 密度公式逐节点最大误差 {max(a['max_density_error_kg_m3'] for a in audits):.4e} kg/m³；Sw(p) 最大误差 {max(a['max_sw_error'] for a in audits):.4e}。这两项支持输入—字段一致性，不证明未保存的历史二进制全部实现路径。",
      f"- BOX 节点体积与单元四角面积积分的初末库存最大相对差 {max(a['max_quadrature_relative_error'] for a in audits):.4e}。解析采样场自检也通过，详见 `postprocessing_selfcheck.json`。这验证后处理积分，不等于旧模拟的离散守恒残差已验证。",
      f"- 顶部节点均值/标准差与每步记录吻合，最大绝对差分别为 {max(a['max_sigmatop_mean_error'] for a in audits):.4e} / {max(a['max_sigmatop_sigma_error'] for a in audits):.4e} mol/mol。",
      f"- 在所有竖向层、横向展示采样点上，用实际 profile 线性插值、输入模态和幅值重建初始 x，最大误差 {max(float(r['initial_profile_sample_max_x_error']) for r in runs):.4e} mol/mol，详见 `initial_protocol_check.csv`。",
      f"- PVD 仅以约 6 位有效数字保存时间；最大舍入差 {max(float(r['absolute_difference_s']) for r in timing):.6g} s。全部 {len(timing)} 帧通过唯一最近时刻、六位舍入回算和顶部统计一致性核对，采用 sigmatop 的高精度时间。原始 `structural_review_required` 检查记录保留，映射另存 `time_reconciliation.csv`。",
      f"- 本地收集核验 {transfer['verified_products']} 个分析产品哈希；没有发现完全相同的整条时序副本（重复初态不算重复轨迹）。14 条记录仍只对应 10 个物理案例，不能作为 14 个独立土壤样本。",
      '', '## 各运行的前阈库存变化', '',
      '“前阈终点”是首次在保存场中观测到全域 x>=0.1 之前的最后一帧；没有观测越界时采用记录终点。这不能排除保存间隔内更早、随后回落的局部峰值。下表都是每条轨迹相对自身初始库存的变化，不是相对于无扰动背景的变化；不同窗口不能直接作物理强弱排名。',
      '', '| 运行 | 模态 / 横向单元 | 前阈终点 (day) | 上部 0.1H 库存变化 | 上部 0.25H | 上部 0.5H | 前阈末帧表面 max(x) |',
      '|---|---:|---:|---:|---:|---:|---:|']
for r in runs:
    values=[100*float(iv[r['id'],region]['relative_inventory_change']) for region in ['upper_0.9','upper_0.75','upper_0.5']]
    text.append(f"| {r['id']} | {r['mode_nx']} / {r['cells_x']} | {float(r['last_pre_observed_threshold_day']):.5f} | {values[0]:+.4f}% | {values[1]:+.4f}% | {values[2]:+.4f}% | {float(r['top_x_max_last_pre_threshold']):.7f} |")
text+=['', '表中深度是目标标签，实际离散切面 z/H 分别为 0.8995154060、0.7522817589、0.5005319057。所有库存用 `phi Sw x`，单位 m，除以同一水平宽度。',
       '', '顶部和侧面零盐交换、暂忽略未测残差时，库存反推累计向下净交换为 `C_proxy=I(0)-I(t)`。本表正的库存变化对应负的 C_proxy，即这些控制体在记录期间净累积盐，而非净失盐。它不提供平流/扩散分解，也不能替代成对背景的 Delta C。',
       '', '## 扰动增长是否已经转化为强再分配', '']
best=next(r for r in runs if r['id']=='r06_c36877_n60')
text += [f"case 36877 的模态 60 是这批中最值得追踪的辅助记录：记录终点为 {float(best['final_time_day']):.5f} day，顶部标准差相对记录中的最小值放大 {float(best['sigma_amplification_full_record']):.4f} 倍。但末帧 sigma={float(best['top_sigma_last_pre_threshold']):.4e} mol/mol，仅为表面均值的 {100*float(best['top_sigma_over_mean_last_pre_threshold']):.5f}%。绝对扰动仍很小，放大倍数不能直接解释成成熟指状对流或强输运。",
         '',f"其三个上部区域库存增加分别为 {100*float(iv[best['id'],'upper_0.9']['relative_inventory_change']):.4f}%、{100*float(iv[best['id'],'upper_0.75']['relative_inventory_change']):.4f}%、{100*float(iv[best['id'],'upper_0.5']['relative_inventory_change']):.4f}%。原 DuMuX t_detect 检测区间与增长曲线已画在同一时间轴；库存没有因出现 std 极小值而自动转为下降。",
         '',f"全体前阈保存帧的最大横向层标准差/该层均值为 {100*summary['max_relative_layer_sigma_pre_threshold']:.5f}%，且可能由初始扰动而非后期增长决定。需要结合各运行末值和整条时间线判断，不能把全时段最大值当成非线性增长幅度。",
         '', '模态 60 采用 1200×200 网格；同一案例模态 30/5 的记录采用不同模式或时间设置，原主协议模态是 28。这些轨迹不能拼接，也不能只凭它们之间的差别宣称某个物理参数的因果作用。',
         '', f"该记录仅保存到 {float(best['final_time_day']):.3f} day；此案例的 LSA 基态阈值参考约为 {float(best['lsa_sat_reference_s'])/86400:.3f} day。这个对比说明旧记录未覆盖全部候选后续窗口，但 LSA 时钟不是已验证的二维实际阈值，不能据此直接设定延长终点。原运行以 DuMuX t_detect 检测为目标、扰动很小，不能默认已经进入强非线性输运阶段。",
         '', '## 阈值与已有时间步信息', '']
text += ['', '下表将每条运行的三个参考时间与实际初始化模态分开；所有时间单位 s。弱区间、immediate 和记录内未检测到的情况不赋予已分辨的 t_detect 点值。', '',
         '| 运行 | 全局 n_c | t_crit | observable n_x | t_seed(0.10) | 实际模态 | 实际模态 crossing | t_detect（resolved） | 检测类别 |',
         '|---|---:|---:|---:|---:|---:|---:|---:|---|']
for r in runs:
    fmt=lambda k: f"{float(r[k]):.2f}" if r[k] else '—'
    text.append(f"| {r['id']} | {r['lsa_global_critical_mode_nx']} | {fmt('lsa_t_crit_s')} | {r['lsa_observable_mode_nx_theta010']} | {fmt('lsa_t_seed_theta010_s')} | {r['initialized_mode_nx']} | {fmt('lsa_t_initialized_mode_crossing_s')} | {fmt('dumux_t_detect_s')} | {r['dumux_detection_class']} |")
text += ['', 'Case 36877 的 canonical observable 模态为 28，t_seed(0.10)=45474.35 s；本表的模态 30/60 运行不实现该 canonical 初始化。t_eff 不从任何一列复制或推断。', '']
for r in runs:
    if r['threshold_high_s']:
        text.append(f"- {r['id']} 首次观测全域阈值的相邻保存帧区间为 [{float(r['threshold_low_s']):.6f}, {float(r['threshold_high_s']):.6f}] s；{r['post_observed_threshold_frames']} 帧从主汇总排除。它是稀疏场的首次观测区间，不是精确首次局部达阈时间。")
comp=rows('shared_window_comparison.csv') if (OUT/'shared_window_comparison.csv').exists() else []
if comp:
    text += ['', '同物理参数、同模态、同网格且 profile 字节一致的既有运行，允许在重叠保存时刻比较不同时间设置；具体差异见 `shared_window_comparison.csv`。']
    for key in sorted(set((r['a'],r['b']) for r in comp)):
        c=[r for r in comp if (r['a'],r['b'])==key]
        text.append(f"- {key[0]} 与 {key[1]}：{c[0]['common_saved_times']} 个公共保存时刻，库存差/初始库存最大 {max(float(r['max_inventory_difference_over_initial']) for r in c):.4e}。差异包含不同时间步的影响，仅适用于这一重叠窗口，不能移用为模态 60 或后续非线性阶段的误差上限。")
text += ['', '## 对下一步的具体建议', '',
         '1. 继续使用旧数据做筛选，不重跑 400 多条运行。这次已经证明它们能够支持库存演化和弱扰动阶段诊断。',
         '2. 项目另一路工作已启动 case 36877 的 600×200、136423 s 无外加扰动试算（`../unperturbed_36877_20261102/RESULT_ZH.md`，本次收尾读取时仍在计算）。优先利用其完成后的背景结果，与原 600×200 主协议模态 28 在共同窗口核对库存差；本次不重复提交背景计算。该试算不属于本次 14 条历史记录，也不能直接充当 1200×200、2.605 day 的模态 60 轨迹的严格配对背景。',
         '3. 现有历史运行缺少独立面通量预算与可执行文件身份；如需发布级机制结论，少量候选的扰动腿也应在可追溯隔离模型中验证，而不能仅给旧场补一个背景就宣称完整闭合。',
         '4. 只有在配对差异超过时间步、网格和预算误差后，再决定是否延长到更强非线性阶段以及定义 t_eff。目前不适合直接训练“有效输运”代理标签。',
         '', '本轮没有实施这些新模拟。结论是“已有时段尚不足以建立显著额外输运”，不是“这些物理条件永远不会产生有效输运”，也不是对 Option A 的最终否定。',
         '', '## 文件与复现', '',
         '- `run_assessment.csv`：14 条运行身份、协议、窗口、阈值和场核验状态；新增标准化字段，旧字段作为来源别名保留。旧 `t_seed_s` 是案例的 canonical observable crossing，不一定是该行实际初始化模态的 crossing，后者使用 `lsa_t_initialized_mode_crossing_s`。',
         '- `inventory_timeseries.csv` / `inventory_summary.csv`：三个区域全部库存时序与汇总。',
         '- `event_intervals.csv` / `time_reconciliation.csv`：阈值记录和时间精度修复证据。',
         '- `figures/inventory_overview.png`：全部运行的库存轨迹；`dense_history_diagnostics.pdf`：逐运行时间线及代表场/剖面。',
         '- `evidence/runs/`：逐帧指标、profile、原小型文件和哈希；`representative_fields.npz` 为仅供展示的横向抽样场，忽略于 Git，库存统计使用所有节点。',
         '- 原始场保留在 `plan.json` 指定的 frodo 路径。r01 分析失败仅为 numpy 布尔值 JSON 序列化错误，日志与原分析产品保存；r02 修复不改变物理设置或计算公式。',
         '', '在 `paper/2` 下执行以下命令可重建本地汇总与图，不启动模拟：',
         '', '```powershell', 'python scripts/analyze_dense_history.py', 'python scripts/write_dense_history_report.py', '```',
         '', '详细公式、符号和证据边界见 `METHOD_ZH.md`。']
(OUT/'REPORT_ZH.md').write_text('\n'.join(text)+'\n',encoding='utf-8')
print(OUT/'REPORT_ZH.md')
