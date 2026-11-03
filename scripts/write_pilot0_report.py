"""Write the bounded baseline research report and an unlaunched WP1 case design."""
from pathlib import Path
import csv,json,hashlib,math,shutil
from prepare_transport_run import parse
ROOT=Path(__file__).resolve().parents[1];O=ROOT/'analysis/pilot0_integrated_20261102'
a=json.loads((O/'baseline_analysis.json').read_text());d=a['detection'];R=a['regions']
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(name,text):(O/name).write_text(text,encoding='utf-8',newline='\n')
def dump(p,obj):p.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')

# Source-verified second-case design. This does not create remote jobs.
C=O/'case38852_design';src=C/'source';frozen=ROOT.parent/'1/analysis/wp1_paired_v2_raw/runs/case-38852--both--nx-16/case_manifest.json'
assert sha(frozen)=='f54d6cf75a93204203906cd2de54c162fbfc735c2dd0b3e50d8c5f14f1560070'
manifest=json.loads((src/'case_manifest.json').read_text());assert manifest==json.loads(frozen.read_text())
source_record=dict(status='verified_by_semantic_manifest_and_physical_mapping',initial_exact_manifest_hash_check='failed_due_to_serialization_difference; preserved source bytes',frozen_manifest_sha256=sha(frozen),remote_manifest_sha256=sha(src/'case_manifest.json'),manifest_json_semantically_equal=True,sha256={p.name:sha(p) for p in src.iterdir() if p.is_file()})
dump(C/'source_hashes.json',source_record)
v=parse(next(src.glob('*.input')).read_text());real=manifest['realization']
mapping={'Permeability':'permeability_m2','Porosity':'porosity','Swr':'residual_saturation','PcN':'vg_n','PcAlpha':'alpha_pa_inv','SweL':'mualem_l','WaterLevel':'water_level_m','EvaporationRate':'evaporation_m_s','SaltDiffusivity':'salt_diffusivity_m2_s','ReferenceViscosity':'reference_viscosity_pa_s','SaltMoleFraction':'initial_salt_mole_fraction','GravitationalAcceleration':'gravity_m_s2'}
for key,val in mapping.items():assert math.isclose(float(v['Problem.'+key]),real[val],rel_tol=1e-13,abs_tol=1e-15)
x=float(v['Problem.SaltMoleFraction']);mw=float(v['Problem.MolarMassH2O']);ms=float(v['Problem.MolarMassNaCl']);X=x*ms/((1-x)*mw+x*ms)
rho=float(v['Problem.PureWaterDensity'])+1000*X*(float(v['Problem.DensityLinearCoefficient'])+float(v['Problem.DensityQuadraticCoefficient'])*X)
assert abs(rho-manifest['constitutive_audit']['reference_density_kg_m3'])<1e-10
assert v['Problem.KrwMualemExponent']=='true' and float(v['Problem.WaveNumber'])==16
assert math.isclose(float(v['TimeLoop.TEnd']),manifest['dumux_numerics']['physical']['t_end_s'],rel_tol=1e-13)
plans={}
for role in ('M0','M1'):
    dest=C/('planned_'+role);dest.mkdir(exist_ok=True);p=dict(v)
    p.update({'Problem.Name':'pilot0_design_38852_'+role,'TimeLoop.DtInitial':str(float(v['TimeLoop.DtInitial'])/4),'TimeLoop.MaxTimeStepSize':str(float(v['TimeLoop.MaxTimeStepSize'])/4),'TimeLoop.OutputInterval':'10','Transport.StopAtSolubility':'true','Transport.SolubilityMoleFraction':'0.1'})
    if role=='M0':
        p.update({'Problem.PerturbationType':'Constant','Problem.Amplitude':'0','Problem.StandardDeviationSaltMoleFraction':'0'});p.pop('Problem.EigenprofilePath')
    else:p['Problem.EigenprofilePath']='eigenprofile.dat';shutil.copyfile(next(src.glob('*.dat')),dest/'eigenprofile.dat')
    sections={}
    for key,value in p.items():sec,key=key.split('.',1);sections.setdefault(sec,{})[key]=value
    text='\n\n'.join('['+sec+']\n'+'\n'.join(k+' = '+val for k,val in vals.items()) for sec,vals in sections.items())+'\n'
    (dest/'params.input').write_text(text,encoding='utf-8',newline='\n')
    plan=dict(case_id=38852,role=role,status='design_only_not_launch_ready',solver_executed=False,parameters=p,input_sha256=sha(dest/'params.input'),source_hashes=source_record['sha256'],reference_density_kg_m3=rho,intended_elf_sha256='d78ce5122cf69e4dee5a4d12bcc726779959f304b343246009cb45bb37e894cc',actual_cut_geometry='reuse 600x200 original grading; verify before execution',required_gates=['case36877 shows resolved transport and reliable diagnostics','extend case-specific matched-pair validator to38852','short resource/cost preflight before long run'],wall_limit_s=43200)
    if role=='M1':plan['profile_sha256']=sha(dest/'eigenprofile.dat')
    dump(dest/'design.json',plan);plans[role]=plan
allowed={'Problem.Name','Problem.PerturbationType','Problem.Amplitude','Problem.StandardDeviationSaltMoleFraction','Problem.EigenprofilePath'}
assert all(plans['M0']['parameters'].get(k)==plans['M1']['parameters'].get(k) for k in set(plans['M0']['parameters'])|set(plans['M1']['parameters']) if k not in allowed)
dump(C/'design_checks.json',dict(passed=True,physical_mapping_to_WP1=True,reference_density_verified=True,pair_differences_limited_to_initialization=True,solver_executed=False,nx=16,dt_s=float(plans['M1']['parameters']['TimeLoop.MaxTimeStepSize']),sample_steps_expected=1600,profile_sha256=plans['M1']['profile_sha256'],scope='Input design checks only; no initial-field or transport validation of case38852'))

# Preserve the candidate extension explicitly as a costed, unlaunched design.
T=136423.0637723484
extension=dict(status='costed_design_not_executed',reason='Weak localized growth and residual-scale deep transport do not justify an automatic long rerun',t_end_candidate_s=2*T,from_t0_required=True,amplitude=1e-6,nx=28,grid=[600,200],dt_max_s=227.3717729539139,stop_x=.1,cost_scenarios_hours=a['cost_scenarios'],growth_scenarios=a['growth'],no_guaranteed_nonlinear_transition=True,guard='Do not execute 4T M0: baseline-scaled15.49h exceeds12h; do not bypass by unvalidated restart')
dump(O/'extension_design.json',extension)

lines=[]
for key,vv in R.items():
 last=vv['final'];lines.append(f"| {last['actual_z_over_H']:.6f} | {vv['max_abs_delta_C_m']:.3e} | {vv['max_abs_delta_C_over_Iref']:.3e} | {vv['max_budget_envelope_m']:.3e} | {last['delta_C_m']:.3e} |")
write('PILOT_RESULTS_ZH.md',f'''# Pilot-0 基准配对结果：检测到弱扰动增长，额外输运尚未分辨

Case36877 的同主机 M0/M1 均正常完成600步、136423.064 s（37.895 h，1.579 day）。正式配对、完整初始场、输出哈希与归档检查通过，M0 跨主机桥接的预算/表面诊断差为0。两条腿各601个物理时刻完全一致，无时间插值。M0/M1墙钟分别{a['integrity']['M0']['wall_s']/3600:.3f}/{a['integrity']['M1']['wall_s']/3600:.3f} h。

## 直接观测

新M1完整sigma序列经冻结Paper1提取器独立计算，t_detect≈{d['onset_refined_s']:.3f} s（约14.00 h），分类为resolved_onset。采样最小值为{d['onset_s']:.3f} s，0.1%平坦区间为[{d['onset_low_s']:.3f}, {d['onset_high_s']:.3f}] s；这个区间是定位不确定性描述，不是置信区间。新检测值恰与历史参考一致，但来源是本次实际输出，未复制历史标签。

LSA全局临界时间31409.158 s、observable nx28时间45474.355 s与上述检测时间保留独立含义。sigma从约7.077e-7下降到6.659e-7，末值7.907e-7，仅为最小值的{d['amplification_since_min']:.4f}倍。末步sigma/表面均值={a['surface']['final_sigma_over_mean']:.3e}，约0.00407%。共同nx28基底的场投影也支持弱模态增长。

全程局部最大浓度{a['surface']['max_domain_x']:.8f}，低于停止阈值0.1；表面及内部阈值事件均未观测，时间被本次终点右删失。末步M1−M0表面均值差{a['surface']['final_mean_difference']:.3e}，最大值差{a['surface']['final_max_difference']:.3e} mol/mol。最大值差主要反映已有横向扰动，不能当作平均盐积累的明显降低。

## 三个原定观测面

守恒量为phi Sw x；单位水平宽度的归一化库存及累计通量单位m，不是真实kg或mol盐质量。向下为正。

| 实际z/H | max abs(Delta C), m | max abs(Delta C)/I_ref | max(abs(R0)+abs(R1)), m | 末步Delta C, m |
|---|---:|---:|---:|---:|
{chr(10).join(lines)}

从第一个接受步起，M0/M1三个切面的瞬时净通量均向上，累计值最终约−1.49至−1.50e-6 m。t=0的两个切面出现约1e-21 m/s的正值，属于启动时近零值，不能读为对流迁移。配对差值约1e-18至1e-17 m，与预算闭合残差同量级，且未经过dt/网格敏感性检验，不能赋予其正负号稳定的物理含义。

并非所有时刻差值都低于同一时刻的闭合包络：三个面分别有4、0、58个样本超过瞬时abs(R0)+abs(R1)。闭合包络不是严格误差界，瞬时越过它不能单独证明可分辨输运；这里保留全时序而不选择有利片段。

## 空间结构与适用范围

26帧原始场已处理。M1扰动主要位于表面以下数毫米，保留初始化的细小横向周期结构，未见明显向深部发展的非线性指状结构。最浅预算面约在20.097 mm深处，末帧其下的最大横向浓度偏差仅约9.44e-15；更深面下约2e-17。因而本结果同时受发展时间和观测深度限制，不能推出所有位置、所有时刻都没有物理效应。

当前最接近Task4的A类：扰动及弱增长可辨认，但相对背景的额外盐输运尚未分辨。B类的稳定输运差异、C类的明显扰动相关对流输运均未建立，D类浓度越界未发生。未单独记录守恒面上的扰动协方差，Delta adv不能冒充该量。

## 机制判断边界

- 检测转折已定位，但t_eff未定义，不能计算“达到有效输运的延迟”。
- 原定深度的总输运以向上对流供应为主；额外adv/diff差值均未独立分辨，不能据其大小证明对流增强或扩散抑制。
- 此处未显示可可靠辨认的adv/diff补偿机制；深部扩散接近零不代表表面薄层扩散不重要。
- 上部库存相对背景积累更慢尚未成立：修正库存差与闭合残差同量级，符号不构成可靠趋势。
- 未观察到可分辨的额外输运逐渐到达更深三个切面；场分析提示发展主要局限于浅层。
- 局部最大浓度与面平均响应确实不同：前者保留约1e-6的模式峰值，后者差仅约1e-11；不能用局部峰值替代平均盐积累或阈值到达时间。

## 六项研究判断

1. 尚未观察到可可靠解释为物理效应的三深度额外盐输运；存在数值非零差值，不能将其称为物理零。
2. 该差值未显示明显高于可确定的闭合残差量级，且U_dt/U_z/U_cut/U_x均unknown。
3. 已证明本窗口中“检测到扰动转折”不等于“已分辨有效输运”；尚不能定量给出t_eff或其相对LSA的时差。
4. Option A可继续做诊断设计和受控发展窗口研究，但尚无支持中心物理主张的正向证据。
5. 本轮不机械加密或启动第二案例；延窗成本和预期收益见NEXT_STAGE_DECISION_ZH.md。
6. Phase2系统参数研究及代理训练NO-GO。保持科学结果inconclusive，而非宣称稳定或无物理效应。

四类图及阅读约定见FIGURE_CAPTIONS_ZH.md；完整表为transport_timeseries.csv，计算摘要为baseline_analysis.json。
''')
write('NUMERICAL_VALIDATION_ZH.md','''# 数值验证与误差边界

## 已完成

- 同sal、同ELF/源码、同600×200分级网格、同物理参数、本构、边界、时间步与终点。仅初始化及输出名称有允许差异。
- M0/M1各600步、601时刻；每个区域完整预算、事件和停止记录均通过正式接口，returncode=0且wrapper无失败。
- 实际t=0网格/压力/Sw相同；M1符合原nx28本征函数和1e-6幅值；独立BOX权重重构初始库存通过。分析中扣除实测初始库存差。
- M0完整桥接的库存、累计通量、残差及表面诊断与原frodo逐值相同。这是可复现性验证，不是物理有效性验证。
- 26帧VTU、PVD和小型输出已回传；归档及成员SHA与finished.json核对一致，见sal_migration/field_archive.json及baseline_analysis.json。
- 原PVD时间只有有限小数；按已核验OutputInterval=50、13帧序号对应的接受步恢复精确物理时刻，所有舍入差≤0.5 s，映射保存在field_time_mapping.csv。没有对预算做插值。
- 每条累计预算独立重新累加dt×瞬时面通量；三深度全部601对满足D_I+Delta C=R1−R0，代数差为0。该恒等式由同一组量构成，不是独立的精度或收敛证明。
- 共同nx28基底对M0/M1原始场重新投影；不把M0 Constant分支的内置投影当作nx28。M1重投影与原每步模态输出在保存帧一致。

## 未完成/unknown

| 误差项 | 状态 | 解释 |
|---|---|---|
| U_dt | unknown | 未运行本次成对dt/2；旧单腿试验不能替代 |
| U_z | unknown | 未运行新成对垂向加密 |
| U_cut | unknown | 坐标预检量化切面移动，但未量化其输运误差 |
| U_x | unknown | 未运行横向加密 |
| 时间对齐 | 无插值、接受时刻一致 | 不包含时间离散误差 |
| 初始库存 | 实测并修正 | 不把其数值小误认为严格零 |

Task5A/5B/5C均未执行求解。理由是三深度差值处于闭合残差量级、扰动仍弱；依command/3先评估发展窗口，不机械执行全部加密。这不等于敏感性检查通过。

已有vertical_cut_precheck.csv仅为坐标嵌套预检。400层最近名义切面相对基准移动约+0.1196、−0.2683、−0.5226 mm；不保证同一控制体。600层可表示原切面，但自动选取名义深度仍可能选错，需要实际映射。没有据此宣称网格收敛。

## 诊断范围

三个面保存的是总守恒通量及adv/diff分解，未保存扰动协方差通量。扩散项在这些远离表面薄层的切面接近零，不能推广为全域扩散无关。归一化phi Sw x收支闭合不等于真实质量收支或模型物理验证。
''')
write('NEXT_STAGE_DECISION_ZH.md',f'''# 下一阶段决策与成本

基准M0/M1、正式配对分析和Figures A–D completed；额外输运的物理结论inconclusive，t_eff undefined。Phase2、批量加密及第二案例执行NO-GO。没有新增求解任务。

## 延窗候选及为什么暂不启动

本次T=136423.064 s。保持原幅值、nx28、网格和dt，将双方从t=0重算至2T=272846.128 s（3.158 day）是一个明确有界候选；本轮未验证restart，不能拼接旧VTK或只延长M1。实际x≥0.1仍立即停止。

按完整基准墙钟作线性成本情景（不是保证）：

| 终点 | M0 h | M1 h | 合计CPU核时 | 判断 |
|---|---:|---:|---:|---|
| 2T | 7.75 | 5.29 | 13.04 | 单腿在12 h内，但预期物理增益有限 |
| 3T | 11.62 | 7.94 | 19.55 | M0几乎无12 h余量，不能自动提交 |
| 4T | 15.49 | 10.58 | 26.07 | M0超过预算，禁止直接启动 |

后半段/末四分之一的log(sigma)线性拟合增长率约2.42e-6/3.18e-6 s^-1，对应倍增时间约3.32/2.52 day。若保持该增长率，到2T的sigma只有约1.10–1.22e-6，较当前约1.39–1.54倍。这个外推仅用于判断成本收益，不是LSA预测或非线性演化保证。达到当前sigma的10倍对应总物理时间约8.60e5–1.09e6 s；按相同比例M0需约24–31 h，超出现有预算，而且10倍本身也不定义有效输运。

因此本轮只登记2T候选和成本，不为追求阳性结果立即追加约13核时。没有将LSA的11.33 day阈值直接当作可安全运行的终点。

## 优先的下一步

先补充面向表面数毫米薄层的观测设计：保留原三个深度作为固定主要结果，同时预先登记更浅的守恒控制面和真正的扰动相关面通量；对新增诊断做短窗口桥接后，再判断2T配对是否有足够信息增益。这是改进测量能力，不能把新深度替换原负结果或事后挑选阳性切面。当前模型未改动。

若未来出现稳定、超出误差量级的输运差值，再按顺序执行成对dt/2、保持真实切面的垂向检查；只有明显指状结构需要时才做Nx敏感性。所有未知误差仍unknown。

第二案例38852已完成WP1输入/本征函数核验与匹配草案，见case38852_design。由于36877尚未提供明确可分辨输运，按Task7仅设计、不运行。正式报告完成后停止本轮自动轮询，保留聊天和所有原始证据。
''')
write('CASE_COMPARISON_ZH.md','''# 两案例的证据与设计

| 项目 | 36877 | 38852 |
|---|---|---|
| 当前状态 | 同主机600步M0/M1完成、配对分析完成 | 仅WP1来源核验和匹配输入设计，未模拟 |
| canonical observable nx | 28 | 16 |
| LSA t_crit / t_seed, s | 31409.158 / 45474.355 | 147558.332 / 147558.332 |
| 本次DuMuX t_detect | 约50394.334 s；平坦区间44110–56616 s | 无新结果；旧WP1为weak_interval_onset，不能当作有效点标签 |
| LSA阈值参考 | 约11.33 day，不等于2D终点 | 151921.412 s，距canonical t_seed约4363 s |
| 输运证据 | 三深度差值未分辨 | 无新守恒输运数据 |
| 扩算决定 | 不机械加密，先评估薄层诊断与有界延窗 | 执行NO-GO，仅保存设计 |

38852使用frodo的paper1_wp1_paired_v2/runs/case-38852--both--nx-16输入和profile。远端manifest字节哈希与冻结本地副本不同，首次精确哈希检查因此拒绝；进一步验证JSON语义完全一致，两个哈希与原字节均留存。输入的12个物理参数逐项映射到WP1 realization，参考密度1006.024104433543 kg/m3与本构审计一致；没有使用hf1密度冲突旧场。

草案M0/M1保留600×200、分级−0.98、nx16/幅值1e-6（M0为0）、统一指数krw和流体常数。拟定首窗口151921.368 s为旧WP1窗口上限，实际表面/全域x≥0.1即停；该上限不是已验证的物理阈值。事件区间仅约4363 s，设计dt为旧379.803/4=94.951 s，约1600步，每步诊断、每10步场输出。该分辨率只是候选，不构成时间收敛。

单腿墙钟上限12 h；尚无该输运ELF下38852的成本测量。按36877每步代价仅作量级估计，1600步可达约7–10 h，干燥/临界附近的迭代与步长退缩可能更昂贵；执行前需要短运行实测、内存/磁盘核验和独立场一致性验收。当前case-specific配对接口仅接受36877，也必须先受控扩展，不能称草案为可直接生产运行。

两个案例同时改变多项土壤/强迫参数及模式，不支持独立单参数因果结论。
''')
write('FIGURE_CAPTIONS_ZH.md','''# Figures A–D：问题、证据与局限

本组为内部研究报告图，不宣称满足某个尚未指定期刊的排版要求。PNG/PDF均已导出。所有曲线使用完整时序，无挑选样本；绘图和来源哈希见figure_manifest.json，重现脚本为scripts/analyze_pilot0_baseline.py。

- Figure A：LSA两个时间、实际新检测、扰动和表面浓度在同一物理时间轴上。上图使用对数纵轴区分背景噪声和施加扰动；中图保留弱衰减/增长。绿色阴影是冻结提取器的平坦区间，不是CI。检测发生不意味着输运有效。
- Figure B：每行一个实际BOX切面，依次给出M0/M1总累计通量、Delta C与负的初始修正库存差、adv/diff差值。灰带仅为abs(R0)+abs(R1)，不是总数值误差带。不同面在差值图上使用不同纵轴范围且明确标注；不能靠微小曲线的形状赋予机制。
- Figure C：t=0、12.63 h、37.90 h的原生节点浓度。上排是完整0.6×0.2 m域；下排明确放大顶部5 mm并减去逐水平行均值，色标量级1e-6 mol/mol。各排共用色标，未做空间插值；放大后的条纹不能称为已形成宏观非线性指状结构。
- Figure D：浓度与0.1阈值、局部最大/均值差异、归一化额外通量及扰动垂向分布。原浅切面约20.1 mm深，远于主要扰动薄层。图中10^-12量级的库存比例是计算差值，不是新t_eff阈值。

已逐张检查PNG：文字和单位可读，图例及色标未遮挡关键数据，无裁切；PDF页数/尺寸和非空输出另见figure_qa.json。灰带以外的少量样本如实保留；不存在统计重复，未加置信区间。U_dt/U_z/U_cut/U_x均unknown。
''')
# Flat comparison table with explicit fresh versus historical event provenance.
with (O/'ONSET_TRANSPORT_COMPARISON.csv').open(newline='') as f:old=list(csv.DictReader(f))
for r in old:
    r['dumux_t_detect_low_s']='';r['dumux_t_detect_high_s']='';r['detection_source']='';r['event_censor_time_s']=''
    if r['case_id']=='36877':
        r['status']='completed_pair_validated_transport_inconclusive';r['solubility_event_status']='right_censored_at_requested_horizon';r['event_censor_time_s']=str(T)
        if r['run_role']=='M1':r.update(dumux_t_detect_s=str(d['onset_refined_s']),dumux_detection_class=d['onset_class'],dumux_t_detect_low_s=str(d['onset_low_s']),dumux_t_detect_high_s=str(d['onset_high_s']),detection_source='new sal M1 full sigma trace; frozen extractor SHA51e75a7f')
    else:r['status']='WP1_design_only_not_executed'
with (O/'ONSET_TRANSPORT_COMPARISON.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(old[0]));w.writeheader();w.writerows(old)
print('Wrote six deliverables, figure captions, extension cost design, and WP1 paired input design; no simulations launched.')
