"""Evidence-derived WP8 reports; never launches simulations."""
import hashlib,json
from pathlib import Path
from analyze_wp8 import ROOT,OUT,RUNS,rows

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def write(name,s):(OUT/name).write_text(s.strip()+'\n',encoding='utf-8')
def table(headers,data):return '| '+' | '.join(headers)+' |\n|'+ '|'.join(['---']*len(headers))+'|\n'+'\n'.join('| '+' | '.join(map(str,r))+' |' for r in data)+'\n'
def main():
    r=json.loads((OUT/'results.json').read_text());gr=rows(OUT/'grid_checkpoints.csv');tr=rows(OUT/'transport_final.csv');cv=rows(OUT/'covariance_checkpoints.csv');a=r['audit'];base='323fd639c81e47b201ecb0ae8b52ba8449892a59'
    write('METHOD_ZH.md',r'''# WP8 方法与解释边界

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
''')
    run_table=table(['运行','帧数','末时刻 s','横向×垂向单元'],[(k,v['frames'],f"{v['end_s']:.9f}",f"{v['nx']}×{v['nz']}") for k,v in a.items() if k!='background_gate'])
    jobs=json.loads((OUT/'jobs.json').read_text());cost=sum(json.loads((OUT/k/'audit.json').read_text())['wall_s'] for k in RUNS)
    write('DATA_AUDIT_ZH.md',f'''# 数据与有效窗口审计

审阅基线 `{base}`；当前指令 `command/12.md`。全程未重跑已有模拟、未启动任何新长模拟；sal任务 {', '.join(jobs.values())} 均为单核读场后处理，总实测读场墙钟相加 {cost:.2f} s（不是新增模拟核时）。

{run_table}
原场共20405帧，每帧重新核对旧field_hashes.json中的SHA256；五组产品回传并复核SHA。R24用本地原归档上传sal解包，R48和背景使用sal已有原场。所有旧源、归档、结果和失败门槛保持原样。Paper1固定来源核验通过。

A各4801个共有状态至545692.2550893656 s。B/C各1200个共有状态，最后136309.37788587352 s；背景末帧136423.0637723484 s与二维末帧136423.06377235046 s超过既有32 ULP约定，因此未强配。这是新的全精度时刻审计，不能用旧12位sigma时间轴掩盖差别。图F在共同窗口结束，2T–4T背景不存在。

原背景最大横向range=3.889943922530392e-14，冻结限值1e-14，**FAILED**；原profile与closure门槛通过并不能抵消此失败。R24/R48的M0也不是严格一维约束解；其晚期噪声增长必须保留。

单腿原库存重构、SCVF积分、停止记录与哈希验收重用WP6/WP7证据；新审计额外核对phi*Sw*x协方差分解、逐帧Parseval和六种配对的库存恒等式。字段缓存与原始VTU留在忽略的tmp目录和远端；小型逐步CSV、图表及来源清单提交Git。

四条二维解均到4T且未触及0.1；只可报告阈值时间右删失于各自末时刻。不同未触及阈值轨迹不能提供达到阈值的时间差。
''')
    mean_table=table(['时刻','R24表面 A','R48表面 A','绝对网格差 mol/mol'],[(str(i+1)+'T',f"{float(x['R24_mean_response_top']):.8e}",f"{float(x['R48_mean_response_top']):.8e}",f"{float(x['absolute_grid_difference']):.8e}") for i,x in enumerate(gr)])
    write('ZERO_MODE_RESPONSE_ZH.md',f'''# 水平平均响应

存在可数值测量的播种平均响应，但晚期量值不稳健。先区分A与C，不能将“可从M0区分”写成“已验证对流因果作用”。

{mean_table}
到3T，A仍仅约−1.5e-9 mol/mol，横向网格差2.29e-11；到4T，R24/R48响应幅值比约382，绝对差5.64001096e-4。4T深度L2响应分别2.30163850e-5和6.00798349e-8 mol/mol，最大响应都在表面。原始表面均值本身接近不能用来宣称这个小差量已经收敛。图A、E分别显示时间演化与浅层剖面。

B在全部共有窗口的最大深度/时间差：R24 3.8771e-15，R48 3.9309e-15 mol/mol，属于浮点/求解噪声尺度。C在共同末时刻：R24 −8.62916683e-12、R48 −8.60204269e-12 mol/mol；这是候选背景上的微小早期响应，背景严格门槛失败仍有效。图F独立展示B/C。

通过T–3T的A同号且横向网格接近，只能说这些已算轨迹在早期数值一致；没有把缺失的长期背景或时间/垂向/域宽误差写成0。4T的R48平均响应并非严格为零，也不是“无下降”，而是比R24小得多。不能据此推导结晶推迟。

selected_mean_profiles.csv保留A/B/C选定时刻的实际深度；mean_difference_energy.csv记录均值差能量。图G的扰动均相对各自水平平均，使用共同色标，避免把横向异常与零模态混淆。
''')
    spec_table=table(['时刻','R24 canonical RMS','R48 canonical RMS','R24 kE','R48 kE'],[(str(i+1)+'T',f"{float(x['R24_canonical_rms']):.7e}",f"{float(x['R48_canonical_rms']):.7e}",f"{float(x['R24_k_E_rad_m']):.6f}",f"{float(x['R48_k_E_rad_m']):.6f}") for i,x in enumerate(gr)])
    write('SPECTRAL_EVOLUTION_ZH.md',f'''# 非零谱与物理波数演化

{spec_table}
RMS单位mol/mol，kE单位rad/m。canonical是m=2，其物理k固定293.215314335。T、2T、3T的canonical RMS网格相对差约0.44%、0.65%、0.91%；3T的kE绝对差0.0098296 rad/m。4T canonical RMS差约8.56%，完整非零能量分别1.07909804e-8与4.95813592e-11 (mol/mol)^2 m，相差约218倍。

4T R24峰值m=10，kpeak=1466.07657 rad/m，能量份额93.46%；R48峰值仍m=2，m=2/6/10分别约35.70%/34.09%/26.46%。R48的kE=817.12285虽不同于kpeak，却正确反映多峰谱。不能仅用“两个网格都出现m=10”宣称短波机制已经收敛。

图B显示kE/kpeak；图C给出全时谱与实际物理k，共同色标和纵坐标，R24未表示的高k区域用阴影标出；图G给出4T浅层二维异常。modal_amplitudes.csv覆盖canonical及m=4/6/10/14/20的全时幅值，modal_checkpoints.csv含每个可表示模式的T–4T完整能量；完整逐帧谱可按SHA恢复缓存。

高次谐波、零模态和增长同时出现，与非线性耦合相容；数值离散及噪声也能激发或污染这些现象。现有谱只能支持“数值轨迹中有晚期非线性响应候选”，不能单独证明物理对流反馈，更不能提供t_eff。
''')
    transport_table=table(['网格','实际深度mm','Delta I_c m','Delta C m','Delta I/I_M0(t) %','Delta C/I_M0(0) %','Delta R m'],[(x['comparison'][:3],f"{(1-float(x['actual_z_over_H']))*200:.6f}",f"{float(x['delta_I_corrected_m']):.6e}",f"{float(x['delta_C_m']):.6e}",f"{100*float(x['delta_I_over_current']):.6g}",f"{100*float(x['delta_C_over_initial']):.6g}",f"{float(x['delta_R_m']):.3e}") for x in tr if x['comparison'].endswith('_A')])
    cov_table=table(['运行','max |cov(Sw,x)|','∫phi cov dz m','max |phi cov/qbar|'],[(x['run'],f"{float(x['max_covariance']):.6e}",f"{float(x['integrated_phi_covariance_m']):.6e}",f"{float(x['max_fraction_of_salt_mean']):.6e}") for x in cv if x['frame']=='4800'])
    write('CONSERVATIVE_TRANSPORT_ZH.md',f'''# 初始库存修正后的守恒盐输运

以下为4T的A=M1−M0；深度使用实际控制体切面，不把1/2/5/10/20/50/100 mm名义值当作精确坐标。各网格使用同一组实际坐标。

{transport_table}
最浅面R24累计差1.87086e-8 m，R48仅5.79200e-11 m；相对初始库存为4.81869%和0.0149182%，相对当时库存对应减少0.457798%和0.00141729%。分母不同会产生不同百分比，不能择大值宣传效果。5 mm及以下变化迅速减小；10–100 mm不足以在当前误差知识下解释为有效输运。图D显示两类归一化量。

六种配对全部共有时刻、七面恒等式max|Delta I_c+Delta C−Delta R|={r['max_identity_error_m']:.9e} m。两浅面输运差明显高于对应数值闭合残差，但这只能确认已算离散轨迹的盐再分配，不能抵消R24/R48差异。源预算本身的盐/水残差仍保存在旧结果和新CSV的ref/case R列；这个恒等式不是一个独立的连续解准确度验证。

以约1 mm面为例，4T R24的M0/M1净J分别−6.47254e-12/−2.45576e-12 m/s；R48分别−6.47983e-12/−6.47141e-12 m/s，均仍向上。正Delta C意味着M1比M0较少向上积累，**不是净向下通量已经为正**。adv/diff/total逐步符号及库存、累计量原值全部在conservative_comparisons.csv，未用幅值取代符号。

## 协方差审计（4T）

{cov_table}
协方差占绝对库存密度比例虽小，仍可能影响微小差量；本轮始终直接平均phi*Sw*x，没有将其设为0。逐帧乘积分解误差小于1e-16。B/C只计算至共有窗口末时刻，结论带有原背景FAILED门槛；没有生成2T–4T的虚构背景输运。
''')
    classification=[['canonical增长','T–3T','Numerically consistent','RMS差0.44–0.91%；仅两档横向网格'],['canonical增长','4T','Grid sensitive','RMS差8.56%；表面与全深投影需各自解释'],['kE/谱','T–3T','Numerically consistent','3T kE差0.00983 rad/m'],['kE/谱','4T','Grid sensitive','kE差683.610 rad/m；总能量差约218倍'],['均值A','T–3T','Numerically consistent','3T绝对网格差2.29e-11；小响应非物理有效性'],['均值A','4T','Grid sensitive','绝对差5.64001e-4；比值约382'],['M0−候选背景 B','≤T共有窗','Below numerical resolution','最大3.93e-15，且背景严格均匀性失败'],['区域库存/累计输运A','4T约1–5mm','Grid sensitive','守恒闭合，但网格差远大于可信收敛余量'],['区域库存/累计输运A','4T约10–100mm','Below numerical resolution','当前dt/z/域宽误差未知，不能用机器闭合当分辨率'],['表面均值/最大值','T–3T','Numerically consistent','原值接近；不意味着小差量已物理验证'],['表面均值/最大值','4T','Grid sensitive','M1均值差5.65234e-4；最大值差0.00614977'],['相对严格无对流背景的响应','>T','Unidentifiable because the reference state is unavailable','缺失长期严格背景；不能用M0替代']]
    save='observable,window,classification,evidence\n'
    import csv
    with (OUT/'observable_classification.csv').open('w',newline='',encoding='utf-8') as f:w=csv.writer(f);w.writerow(['observable','window','classification','evidence']);w.writerows(classification)
    write('GRID_SENSITIVITY_ZH.md','# R24/R48可观测量审查\n\n'+table(['量','窗口','分类','依据'],classification)+'''
这些是证据解释分类，不是事后拟合的自动PASS阈值。“Numerically consistent”仅限当前两档横向网格和指定窗口，绝不等于dt、垂向和域宽收敛。“Below numerical resolution”区分机器噪声与在现有验证范围内不能支持物理效应，两者均不填零；原值全部保留。

答案是否定的：晚期水平平均响应并没有比非零谱更稳定。绝对浓度均值的相对误差会被大背景掩盖，应检查A差量的绝对网格差。区域原库存相近同样不能证明新增库存变化稳健。

grid_checkpoints.csv给出原模式、非零能量、物理k、表面均值/最大值、平均响应的T–4T数值和绝对差。2T的既有dt/2证据不能外推为4T晚期误差界；U_dt(晚期)、U_z及U_width均unknown。仅R24/R48两级也不足以作可靠Richardson外推。
''')
    background=ROOT/RUNS['B'];protocol=dict(status='algebraic_reduced_Newton_validated_native_DuMuX_integration_pending',launch_authorized=False,base_commit=base,base_input_sha256=sha(background/'params.input'),source_main_sha256=sha(ROOT/'model/transport_wp7_localflux_v1/source/appl/simpleRichardsnc/main.cc'),geometry=dict(width_m=.6/28,height_m=.2,cells_x=1,cells_z=200,grading_z=-.98),formulation='same BOX residual F(Pu); reduced P^T F and P^T J P; one pressure and concentration per exact height',input_changes_allowed={'Problem.Name':'WP8-constrained-background-4T','TimeLoop.TEnd':545692.2550893936},dt_s=113.68588647695695,stop=dict(surface_or_domain_max_mole_fraction=.1,check_each_accepted_step=True,record_bracket=True),gates=dict(lateral_range_max=1e-14,expected_lateral_range=0.,initial_x_abs=1e-14,initial_pressure_Pa=1e-8,initial_storage_relative=1e-12,profile_T_abs=1e-10,salt_residual_m=1e-12,water_residual_m=1e-10,storage_reconstruction_salt_m=1e-14,flux_total_abs_m_s=1e-17,inventory_match_abs_m=1e-12,flux_reintegration_abs_m=1e-14),budget=dict(short_validation_steps=10,short_validation_wall_cap_s=120,long_window_stages_T=[1,2,3,4],long_run_wall_cap_s=900,cpus=1,memory_gib=2,local_and_remote_archive_cap_mib_each=250,restart='no claim of restart equivalence; prefer one continuous run with stage checkpoints',authorization='required for new long run; no submission script supplied'))
    (OUT/'background_protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
    write('BACKGROUND_VALIDATION_PLAN_ZH.md',r'''# 严格背景实现与4T验证计划

## 现有模型为何不严格均匀

WP7 background是同一二维YaspGrid/BOX方程的一格横向挤出，仍有左右两套独立压力、浓度未知量。main.cc用FVAssembler和普通NewtonSolver解完整二维系统，没有等高约束；迭代/装配浮点不对称可以留下或增长。既有横向range=3.88994e-14 > 1e-14，不能通过放松门槛把它变成已验证背景。

直接改成真正一维需要处理globalPos[1]、重力方向、边界检测、BOX截面积、输出和七面SCVF语义。为尽量保持原物理与垂向算子，优先采用**一格挤出上的显式等高约束**，而不是每一步解完后把左右平均。

## 已实现与未实现

`scripts/wp8_constraint.py`提供可执行的约束Newton参考实现：u每个高度含一个压力和一个浓度，P将其复制至二维顶点；R(u)=P^T F(Pu)，Jr=P^T J P，Dirichlet值在约束自由度上精确施加，接受步直接保存Pu，不做求解后投影。不等体积贡献在F中保留并求和，内部交换项相消；这与把两个已算二维解取平均不同。

10项数值测试包含非线性Jacobian差分、体积加权非线性储存/内部通量的守恒、边界值、严格逐位均匀及失败迭代不接受。它们验证约束代数与参考Newton实现，**没有验证DuMuX物理求解器**。当前参考代码尚未接入FVAssembler/NewtonSolver；没有编译一个新的DuMuX背景可执行程序，也没有产生新物理解。Task 5的实现原型已验证，原生集成/短窗物理桥接仍pending，不能将4T计划标成可直接生产运行。

下一工程步骤是在独立`model/transport_wp8_background_v1`中接入稀疏约束装配及Newton接口，沿用原本构/SCVF核，不改WP7源码。先逐项核对压力/浓度方程映射、底部Dirichlet消元、积分残差缩放与Newton步拒绝逻辑，再编译，执行不超过10步的桥接；未完成该步不得启动4T。真一维实现作为备选，也必须通过同一桥接，不能仅凭维数更低认定等价。

## 冻结协议

机器可读参数、原输入和源哈希见background_protocol.json。几何W=0.021428571428571432 m、H=0.2 m、横向1格、垂向200格、Grading1=-0.98，保留已有节点坐标与实际七面切割；dt=113.68588647695695 s。沿用原M0初始压力/浓度剖面、密度公式、盐扩散、蒸发、Krw指数形式和底部水位/浓度条件。仅新名称、TEnd及明确记录的约束实现变化允许不同，不改幅值或材料参数。

分级验证：

1. t=0及10步：与原背景及R48 M0相同初始剖面；浓度绝对差<=1e-14，压力<=1e-8 Pa，库存相对差<=1e-12；逐高度压力/浓度复制必须逐位相同，浓度range仍必须<=1e-14。失败保留证据，不调整阈值。
2. 获得长运行授权后，连续0–4T计算，T/2T/3T/4T记录验收检查点。T窗口对原挤出背景的水平均值最大差<=1e-10 mol/mol、同区域库存绝对差<=1e-12 m、净盐通量绝对差<=1e-17 m/s；原挤出背景的FAILED记录仍独立保留。
3. 每阶段盐残差<=1e-12 m，水残差<=1e-10 m；独立VTU盐库存重构差<=1e-14 m，独立累计盐通量再积分差<=1e-14 m。比较使用接受时刻；如终点截步不匹配则丢弃该点，不插值C。保存归一化残差及原值，接近阈值或有系统漂移时人工审查，PASS不覆盖物理异常。
4. 每个接受步检查表面和全域max(x)>=0.1，任一触发立即停止，保存越阈前后括区间；若只是达到4T，阈值时间保持censored。新背景独立停止，不能借M0/M1未触阈推断其结果。
5. 三段以上若盐库存误差仍与目标效应相当，先细化时间/垂向并重新预算，不发布小差量物理结论。后续与二维比较必须记录宽度归一化与边界物理相同。

## 成本、权限和停止

旧一格1200步求解实测8.764 s、后处理5.302 s、归档19.75 MB，线性外推4800步约56.3 s求解+处理、约79 MB原场归档，**这不是约束实现的实测性能**。先用10步桥接重新计时，最多120 s；拟议长跑1核/2 GiB、单条15 min硬上限、原场归档每端250 MiB。计划储存预算含余量，不计既有大数据。超限、数值失败、哈希差异、均匀性或守恒失败都停止并保留最后接受状态，不盲目重试。

本轮没有启动上述短窗或长跑，仅完成离线代数验证和计划；新的原生集成/桥接状态pending，新的长1D、R96、多波长或8T必须另获授权。没有计算资源不足的假定，也没有因此替用户缩减已授权的后处理范围。
''')
    write('SCIENTIFIC_DECISION_ZH.md',r'''# 科学判断与后续优先级

## 最终八问

1. **增长的横向扰动会不会产生可测平均浓度响应？** 会，在现有离散轨迹中可直接测出。到3T表面A约−1.5e-9，4T分别−5.65482e-4/−1.48073e-6 mol/mol；晚期响应尚非已验证物理量。
2. **是否能与无播种二维轨迹区分？** 可以数值区分A；浅层库存/通量也超过单腿闭合残差。但M0晚期本身存在横向增长，A不是“有对流减去严格无对流”。相对严格背景的因果贡献仍缺证据。
3. **哪些量在R24/R48间稳健？** T–3T原模式RMS、kE及小均值A较一致；两个网格都未达到0.1。4T的谱、均值响应、浅层新增库存/累计输运、表面最大值均不稳健。不能因原库存大体接近而宣称新增输运已收敛。
4. **浓度响应是否对应守恒盐再分配？** 约1–2 mm层确有初始修正库存减少和累计通量差相配，恒等式误差<=6.353e-22 m；以当时库存计最浅层R24减少0.4578%，R48仅0.001417%。净通量仍向上。空间浓度起伏很大并不意味着大范围或显著的净盐输出。
5. **4T前有可信的非线性反馈证据吗？** 有“数值解中均值响应/谐波/守恒再分配共现”的证据，与非线性耦合相容；尚没有网格稳健、相对严格背景的物理反馈证据。3T后的强分歧阻止定量解释，不能将谐波出现单独当证明。
6. **哪些结论受长期一维背景缺失限制？** >T的B/C、无播种二维自身偏离、真正背景修正的零模态、对流贡献、背景与二维到0.1时间差。A的数值差不受“无数据”限制，但其因果解释受限。现有短背景严格均匀性还失败，≤T的B/C也只是候选比较。
7. **优先一维背景还是R96？** 优先完成显式等高约束的原生集成和短窗物理桥接，再单独授权严格背景4T。背景便宜且提供当前完全缺失的参照；R96随后仍是晚期网格审查的必要候选，不能被背景替代。暂不开新长运行。并行资源不改变这些逻辑依赖。
8. **可靠t_eff将来如何定义？** 暂为unidentified/null。候选必须是相对已验证背景的、共同物理区域内的初始修正库存响应及保守通量变化，同时超过预登记物理影响阈值与数值误差上界，并在时间窗口内持续，不能用sigma或单个谐波触发。

## 可检验的下一版指标（尚未启用）

可令E_I=−Delta I_c/I_b(t)，E_C=Delta C/I_b(t)，预先指定约1、2、5 mm三个实际平面，选择主要平面而非观察结果后挑最大的一个。建议下一次协议讨论采用**当时区域库存1%**作为可解释影响候选，并要求持续至少0.1T；这是尚待物理问题论证的设计值，不是已验证阈值，也不回溯生成本批标签。

正式判定需要E_I与E_C符号和预算相容、严格背景通过、独立的dt/横向/垂向/域宽对照，并在所有容许协议中“效应减去保守误差界”仍超过阈值。没有证据支持误差独立时不采用平方和缩小误差，可先取各项绝对界之和；未知项保留unknown，使事件不能被识别。没有三档收敛证据时，不把两个网格差冒充Richardson误差界。

LSA临界时刻、播种扰动的增长/观测检测、均值响应出现、可分辨盐再分配、首个0.1局部越阈，是五种不同事件。4T没有越阈只意味着右删失；不宣称结晶延迟，不制造t_eff或surrogate训练标签。

## 本轮交付状态

Tasks 1–4、6–7：既有数据后处理、图表、分类与科学判断已完成。Task 5：约束Newton代数原型及执行/验收协议已准备并通过离线测试；DuMuX原生集成、编译及短窗物理桥接仍pending，不能声称新的严格背景模型已经验证。本轮未改Paper1和既有基准，未新增模拟。
''')
    native=OUT/'background_smoke10_r03/validation.json'
    if native.exists():
        v=json.loads(native.read_text());assert v['passed']
        protocol['status']='native_DuMuX_compiled_ten_step_bridge_passed_long_validation_pending'
        protocol['native_short_validation']=v
        protocol['native_build']='build/r03_attempt_01/result.json'
        protocol['internal_Newton']=dict(water_residual_max=2e-18,salt_residual_max=2e-20,relative_shift_max=1e-10,max_iterations=20,failure='abort; no timestep restart implemented')
        protocol['budget']['measured_native_10step_wall_s']=v['wall_s'];protocol['budget']['linear_4800step_estimate_s']=v['wall_s']*480
        (OUT/'background_protocol.json').write_text(json.dumps(protocol,indent=2)+'\n')
        p=OUT/'DATA_AUDIT_ZH.md';s=p.read_text(encoding='utf-8').replace('所有旧源、归档、结果和失败门槛保持原样。','所有旧源、归档、结果和失败门槛保持原样。')
        s=s.replace('全程未重跑已有模拟、未启动任何新长模拟；','未重跑已有模拟、未启动任何新长模拟；')
        s+='\nTask 5另在frodo完成隔离源编译和有界原生桥接：r01/r02各在第1步接受前失败并归档；r02仅增加残差诊断，没有改变算法。定位舍入噪声后r03完成10步至1136.8588647695694 s，墙钟1.16757 s。全部原场和失败证据均回传核验SHA，未运行T或4T新背景。\n';p.write_text(s,encoding='utf-8')
        p=OUT/'BACKGROUND_VALIDATION_PLAN_ZH.md';s=p.read_text(encoding='utf-8');start=s.index('## 已实现与未实现');end=s.index('## 冻结协议')
        s=s[:start]+'''## 已实现和短窗验收

`model/transport_wp8_background_v1/source`是WP7独立副本，仅main.cc接入等高约束求解和固定步长，新增`wp8_equal_height.hh`。model.hh本构、储存、SCVF核、边界公式、七面诊断与独立溶解度停止均未改动。每个高度共享压力/浓度，Newton残差为P^T F(Pu)，Jacobian为P^T J P；每次迭代直接提升Pu，不做事后平均。原Dirichlet两份方程同样求和，其倍数缩放不改变约束值。仅允许单MPI rank、一格挤出；非均匀初值拒绝，不收敛则异常退出且不接受该步。现实现使用402标量未知量的行缩放、部分主元稠密线性求解，局限于此有界背景。

Python参考实现及10项单元测试覆盖保守约束、Jacobian、边界和失败处理；原生C++在frodo完成三次隔离编译。r01首步失败；r02只添加日志，确定水方程残差在5.3e-19–1.2e-18范围停滞、相对位移已约1e-15，统一1e-20内部判据低于舍入噪声。r03按两方程尺度使用水2e-18、盐2e-20且相对位移<=1e-10；这些是内部Newton判据，不改变1e-14横向均匀性和独立盐/水收支门槛。两个失败运行都保留0接受步记录和原场，没有覆盖或宣称通过。

r03完成10步至1136.8588647695694 s，16项原生检查全部通过：所有11帧压力和浓度横向range均为0；盐残差max=2.09778e-17 m、水残差max=5.82730e-15 m；与原背景库存差1.49078e-19 m、净盐通量差1.48513e-19 m/s，浓度剖面差4.51462e-15 mol/mol。初始压力公式误差5.82e-11 Pa、初始浓度误差0；独立VTU盐库存重构差2.84603e-19 m、累计通量再积分误差0。输入仅Name/TEnd不同，所有物理参数和实际观测面保留。

这是**原生实现和短窗桥接通过**，不是T/4T背景已验证。源、构建、ELF、运行输入、11帧归档及SHA在build/r03_attempt_01、background_smoke10_r03和source_manifest.json。长运行仍须以下分阶段验收。

'''+s[end:]
        s=s.replace('先用10步桥接重新计时，最多120 s；','新原生10步实测1.16757 s，包含进程启动；保守线性放大到4800步约560.4 s，加旧后处理线性估计约21 s，仍须在T检查点更新实测成本。短窗硬上限120 s；')
        s=s.replace('本轮没有启动上述短窗或长跑，仅完成离线代数验证和计划；新的原生集成/桥接状态pending，新的长1D、R96、多波长或8T必须另获授权。','本轮已执行上述10步桥接，未启动新T/4T长背景。新的长1D、R96、多波长或8T必须另获授权。')
        p.write_text(s,encoding='utf-8')
        p=OUT/'SCIENTIFIC_DECISION_ZH.md';s=p.read_text(encoding='utf-8').replace('优先完成显式等高约束的原生集成和短窗物理桥接，再单独授权严格背景4T。','显式等高约束原生集成与10步桥接本轮已通过，下一步优先单独授权严格背景4T，按T检查点分阶段验收。')
        s=s.replace('Task 5：约束Newton代数原型及执行/验收协议已准备并通过离线测试；DuMuX原生集成、编译及短窗物理桥接仍pending，不能声称新的严格背景模型已经验证。本轮未改Paper1和既有基准，未新增模拟。','Task 5：约束Newton原生实现、隔离编译、16项10步桥接检查及执行/验收协议已完成；两个首步失败及诊断修正完整保留。4T长期验证尚未执行。本轮未改Paper1和既有基准，新增模拟仅限有界10步桥接，未启动任何新长运行。')
        p.write_text(s,encoding='utf-8')
    p=OUT/'CONSERVATIVE_TRANSPORT_ZH.md';p.write_text(p.read_text(encoding='utf-8').replace('conservative_comparisons.csv','conservative_comparisons.csv.gz'),encoding='utf-8')
    print('8 reports, classification and background protocol written')

if __name__=='__main__':main()
