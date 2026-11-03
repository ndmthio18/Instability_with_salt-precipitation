"""Command 14 deliverables from frozen evidence; prepares but never runs 4T."""
import csv,difflib,json,io,subprocess,sys
from datetime import datetime
from pathlib import Path
import numpy as np
import wp8_v2_regression as v
from analyze_wp8 import rows
c=v.c;O=c.OUT
def write(name,text): (O/name).write_text(text.strip()+'\n',encoding='utf-8')
def table(header,rr):return '| '+' | '.join(header)+' |\n|'+'|'.join('---' for _ in header)+'|\n'+'\n'.join('| '+' | '.join(map(str,r))+' |' for r in rr)+'\n'
def main():
    summary=json.loads((O/'regression_summary.json').read_text());eq=json.loads((O/'equation_diagnostics.json').read_text());assert summary['all_regression_checks_passed'] and eq['FD_scale_study_pass']
    runs={n:json.loads((O/n/'validation.json').read_text()) for n in v.CASES}; costs={n:r['cost'] for n,r in runs.items()}
    old=c.ROOT/'model/transport_wp8_background_v1/source';new=c.ROOT/'model/transport_wp8_background_v2/source';changed=[];patch=[]
    for p in new.rglob('*'):
        if not p.is_file():continue
        n=p.relative_to(new).as_posix();q=old/n
        if q.read_bytes()!=p.read_bytes():changed.append(n);patch.extend(difflib.unified_diff(q.read_text().splitlines(True),p.read_text().splitlines(True),fromfile='v1/'+n,tofile='v2/'+n))
    assert set(changed)=={'appl/simpleRichardsnc/main.cc','appl/simpleRichardsnc/wp8_equal_height.hh'}
    (new.parent/'v1_to_v2.patch').write_text(''.join(patch),encoding='utf-8')
    late=np.load(c.CACHE/'late_scaled/post/profiles.npz');legacy=np.load(c.CACHE/'late_legacy_window/post/profiles.npz');factor=np.load(c.CACHE/'late_scaled_16/post/profiles.npz')
    comparison={'late_profile_max_diff':float(np.max(abs(late['x_mean']-legacy['x_mean']))),'late_pressure_max_diff_Pa':float(np.max(abs(late['pw_mean']-legacy['pw_mean']))),'factor16_vs32_fields_exact':all(np.array_equal(late[k],factor[k]) for k in late.files)}
    for key in ['inventory_salt_m','cut_salt_total_down_m_s']:
        aa=rows(O/'late_scaled/transport_budget.csv');bb=rows(O/'late_legacy_window/transport_budget.csv');assert [(r['region'],r['time_s']) for r in aa]==[(r['region'],r['time_s']) for r in bb];comparison[key+'_max_diff']=max(abs(float(a[key])-float(b[key])) for a,b in zip(aa,bb))
    comparison['late_bridge_pass']=comparison['late_profile_max_diff']<=1e-10 and comparison['inventory_salt_m_max_diff']<=1e-12 and comparison['cut_salt_total_down_m_s_max_diff']<=1e-17
    comparison['probe_does_not_change_accepted_field']=all(np.array_equal(np.load(c.CACHE/'late_legacy_probe/post/profiles.npz')[k],np.load(c.CACHE/n/'post/profiles.npz')[k]) for n in ['late_equation_probe','late_equation_probe_fine'] for k in late.files)
    c.save(O/'late_protocol_comparison.json',comparison);assert comparison['late_bridge_pass'] and comparison['probe_does_not_change_accepted_field'] and comparison['factor16_vs32_fields_exact']
    tr=rows(c.ROOT/'analysis/wp8_strict_background_20261103/newton_trace.csv');stats=[]
    for label,lo,hi in [('1T',1,1200),('1T_to_2T',1201,2400),('2T_to_3T',2401,3600),('late_accepted',3601,4413),('failed_4414',4414,4414)]:
        rr=[r for r in tr if lo<=int(r['attempted_step'])<=hi];ends={int(r['attempted_step']):r for r in rr};it=[int(r['iteration']) for r in ends.values()];stats.append(dict(window=label,steps=len(ends),iterations_max=max(it),iterations_median=float(np.median(it)),final_water_max=max(float(r['water_residual']) for r in ends.values()),final_salt_max=max(float(r['salt_residual']) for r in ends.values())))
    c.save(O/'historical_trace_windows.json',stats)
    gate=json.loads((c.ROOT/'analysis/wp8_zero_mode_20261103/background_protocol.json').read_text())['gates'];build=json.loads((O/'build/attempt_01/result.json').read_text());base=c.ROOT/'analysis/wp7_shortwave_validation_20261103/runs/WP7-1D-background-T/params.input'
    inp=base.read_text().replace('Name = WP7-1D-background-T','Name = WP8-v2-strict-4T').replace('TEnd = 136423.0637723484',f'TEnd = {4*c.T:.17g}')+'\n[WP8]\nCriterion = scaled\nRoundoffFactor = 32\nMaxIterations = 20\nProbe = false\nAllowDiagnosticReplay = false\n'
    write('pending_full4T.input',inp)
    plan=dict(decision='GO_FOR_USER_AUTHORIZED_FULL_4T_VALIDATION_ONLY',launch_authorized=False,launched=False,scientific_4T_passed=False,reviewed_commit='87678761f4fa0d64a7205add3f1f0139a1b716e8',input_sha256=c.sha(O/'pending_full4T.input'),executable=build['executable'],source_hashes=build['source_sha256'],source_manifest_sha256=c.sha(new.parent/'source_manifest.json'),regression_summary_sha256=c.sha(O/'regression_summary.json'),equation_diagnostics_sha256=c.sha(O/'equation_diagnostics.json'),gates=gate,original_initial_state=True,restart=False,replay=False,dt_s=113.68588647695695,horizon_s=4*c.T,max_steps=4800,checkpoints_steps=[1200,2400,3600,4800],threshold=.1,threshold_check='every accepted step; record first bracket; no later accepted state',criterion={'relative_shift':1e-10,'scaled_residual':32*np.finfo(float).eps,'componentwise_linear_backward_error':1e-12,'regional_step_salt_m':1e-12/4800,'regional_step_water_m':1e-10/4800,'max_Newton_iteration':20},cost={'host':'frodo','cpus':1,'ranks':1,'wall_hard_cap_s':900,'solver_address_space_gib':2,'archive_cap_mib_per_end':250,'minimum_free_disk_gib':2,'estimate_wall_s_range':[180,600]},failure_action='retain outputs; no automatic retry, timestep cut, tolerance tuning or continuation',completion='hash every field; independently reconstruct seven cuts plus free domain; exact x/p equality; frozen T/common-prefix bridge; report censored/failed distinctly')
    c.save(O/'full4T_protocol.json',plan)
    cost=dict(runs=costs,frames=sum(x['frames'] for x in costs.values()),total_accepted_steps=sum(x['frames']-1 for x in costs.values()),solver_and_monitor_wall_s=sum(x['solver_and_live_audit_wall_s'] for x in costs.values()),total_worker_wall_s=sum(x['total_wall_s'] for x in costs.values()),archives_bytes=sum(x['archive']['bytes'] for x in costs.values()),max_solver_rss_kib=max(x['peak_child_rss_kib'] for x in costs.values()),future_4T=plan['cost'],monitor_RSS=None,network_wait_cost=None)
    cost['compile_wall_s']={}
    for n in ['build','probe_build','probe_fine_build']:
        b=json.loads((O/n/'attempt_01/result.json').read_text());cost['compile_wall_s'][n]=(datetime.fromisoformat(b['finished_utc'])-datetime.fromisoformat(b['started_utc'])).total_seconds()
    c.save(O/'cost_audit.json',cost)
    rt=table(['短测','接受步','结果','盐闭合上限 m','水闭合上限 m'],[[n,r['cost']['frames']-1,'预期拒绝/测试通过' if n=='forced_failure' else '通过',f"{r['metrics']['max_salt_closure_m']:.3e}",f"{r['metrics']['max_water_closure_m']:.3e}"] for n,r in runs.items()])
    write('NUMERICAL_FAILURE_ANALYSIS_ZH.md',f'''# WP8 晚期 Newton 停滞：方程级诊断

结论：原第4414步是固定绝对盐残差判据与双精度残差评估尺度冲突导致的停滞；证据支持舍入/消去主导，不支持物理发散或线性解算失败。此判断是对原完整轨迹和邻近状态诊断的综合推断，未得到原失败迭代的无损状态及矩阵，因此不能宣称逐位复现或严格证明全局残差下限。

## 原始失败与近邻回放

旧4T结果保持失败：4413步，501695.8170227767 s = 3.6775T；第4414步未保存为接受状态。原失败盐残差从1.415301e-14降至8.457137e-19，再降至2.406284e-20；后续迭代在约2.126e-20–4.855e-20间振荡，未稳定下降。压力/水残差已低于2e-18，最后修正约1e-15量级。小修正单独不能证明物理解收敛。

`historical_trace_windows.json`从原完整Newton轨迹按T分窗；`failed_step_iterations.csv`仍在原目录，保留所有拒绝迭代。

VTU以文本保存的4413步状态仅用于诊断回放。它与内存精确状态有舍入差异，首轮底部Dirichlet压力残差2.91e-11 Pa在下一迭代消除。该回放单步及10步旧判据均成功，因此**没有重现完全相同的失败**；10步中仍见旧阈值导致额外迭代（第10步直到第8次才接受）。不能把这些诊断场拼接成正式4T背景，也不声称重启等价。

## 线性系统、方程尺度与雅可比

402维约束系统仍为r=P^T F(Pu)，A=P^T J P；原BOX残差/边界/本构不变。行归一化、部分主元高斯消去未换算法。

- 邻近晚期迭代2：线性分量后向误差C++长双精度累加为4.364e-15；重新计算约{eq['linear_eta_component_recomputed']:.3e}。线性残差约{eq['linear_absolute_residual_recomputed']:.3e}，远小于非线性残差。Windows numpy longdouble有效尾数为{eq['numpy_longdouble_mantissa_bits']}位，不能把本地复算声称为额外精度验证。
- 原始带量纲矩阵2范数条件数{eq['raw_condition_2']:.4e}；以变量尺度及行和无量纲化后约{eq['row_column_scaled_condition_2']:.4e}。仍存在敏感性，不能称完全良态；但原始巨大数值主要受量纲和边界行影响，不能直接解读为物理奇异。
- 行归一化后最小主元约1.4007e-3、最大约1.0126、增长因子约1.0126，无灾难性主元增长。线性后向误差和非线性实际下降均不指向线性解失败。
- 非Dirichlet水方程尺度sum|A_ij|s_j为4.8645e-4–8.8843e-3，盐方程为1.7515e-6–4.2771e-4；单一带单位绝对阈值没有反映这些差别。
- 对底部邻近、中部、顶部两层的压力/浓度共8列，分别使用6档中央差分，覆盖全部402个约束残差分量。每列在稳定差分尺度下相对差低于{max(eq['best_fd_relative_by_column'].values()):.3e}。大压力扰动产生最高约1.13%差异，减小扰动后消失；极小浓度扰动也开始受舍入影响。完整表见`equation_diagnostics.json`，没有删除不理想尺度的数据。
- 依赖源码证据表明DuMuX采用数值雅可比，默认前向差分h=1e-10*(|u|+1)。物理通量有上风分支；大压力探针可能跨越分支，有限差分步长也有截断误差。此处分支归因属于解释，未逐面记录切换。测试未发现约束映射/装配错误，但并非穷举402列的雅可比证明。

## 浮点分辨率与非线性停滞

压力各自由节点仅改变一个浮点间隔，盐残差最大变化1.370e-20（向上）/1.440e-20（向下），已相当于旧2e-20阈值的68–72%。浓度同向一ULP扰动的影响约8e-26。后者不能代表随机逐节点扰动的上界。

源代码的压力梯度由约9.2e4 Pa节点值乘形函数梯度相加，再与重力项相减；瞬态库存也作相邻时刻差。观测到的压力ULP敏感性与这类消去一致。r/(|A|s)约1.233e-16，处在双精度epsilon=2.220e-16的量级。反复Newton修正不能强迫离散浮点评估在每个状态都落到同一固定阈值之下。

因此旧2e-20在部分晚期状态低于可可靠要求的残差尺度，而不是在所有状态都不可达到；邻近回放成功正好否定“普遍硬下限”的过强解释。未实施全残差高精度重算，不能精确分解压力消去、通量抵消和存储差分各自的贡献。

## 判据修订与实际回归

采用尺度化非线性残差、线性后向误差、上一步相对修正及每步区域守恒缺陷四重门槛，详见`SOURCE_REVIEW_AND_CONTRACT_ZH.md`。未直接提高盐残差常数，未减小时间步，未修改物理门槛。

{rt}

早期10步与冻结T的11个状态剖面、压力和预算完全相同；T附近10步有10个严格共有时刻（最后终点因浮点截断差异不配对），最大浓度剖面差1.006e-16，净通量差3.167e-19 m/s。晚期修订与旧判据10步剖面差{comparison['late_profile_max_diff']:.3e}、净通量差{comparison['cut_salt_total_down_m_s_max_diff']:.3e} m/s，仍通过原桥接门槛。factor16与32保存场完全相同，诊断探针未改变接受场。

修订晚期盐闭合最大7.886e-17 m，旧判据短测3.679e-17 m；不能声称每项残差都变小。两者均远低于原1e-12 m门槛，场/通量桥接也未发生可分辨的精度损失。长期累积是否仍合格必须由新的连续4T验收决定。
''')
    write('SOURCE_REVIEW_AND_CONTRACT_ZH.md',f'''# 隔离源码审查与停止准则

生产候选：`model/transport_wp8_background_v2/source`，ELF SHA256 `{build['executable']['sha256']}`。

## 逐项源码边界

与冻结v1相比仅`main.cc`和`wp8_equal_height.hh`变化，逐行补丁为`model/transport_wp8_background_v2/v1_to_v2.patch`。`model.hh`、`transportdiagnostics.hh`及其余文件字节完全相同。旧源、旧ELF、T通过结果与4T失败结果均不覆盖。

1. main只增加默认关闭的诊断回放入口，并向约束求解器传递当前dt。原初值生成函数、更新xOld位置、接受步骤、VTU/事件/七面诊断和0.1停止代码不变。生产input明确关闭回放且无ReplayProfile。
2. 按实际z节点建立P，原始两个横向节点必须完全相等；每次Newton状态都从同一u精确复制。求和残差与雅可比块对应P^T F及P^T J P，未使用求解后平均/投影掩盖横向误差。底部Dirichlet行保留原装配。
3. 线性求解的消去顺序、行缩放和主元选择保持v1；增加主元增长和后向误差记录。新版本在判断接受前额外求解一次当前线性系统，用于独立检查；不再以小修正单独接受。
4. `legacy`分支保留原2e-18/2e-20/1e-10逻辑；生产候选显式选择`scaled`。默认仍legacy以免漏写参数时静默改协议。失败到最大迭代直接抛出，未增加隐式重试/步长切减。
5. 三个编译版本都逐文件核验。两个`_probe`版本只改探针触发迭代及差分步长，生产算法相同；它们仅用于诊断，不是最终4T指定ELF。探针恢复原场及装配器，实际保存场与无探针完全相同。

编译退出码均为0，但并非无警告：除依赖未用参数警告外，新头文件第59/79行有两项misleading-indentation警告。逐语句检查确认装配调用位于循环之后、shift更新位于失败判断之后，均符合预期控制流；完整日志保留，未用关闭警告选项掩盖。

## 冻结候选准则

令s_j=max(|u_j|,1 Pa)（压力）或max(|u_j|,0.0036)（摩尔分数）；d_i=sum_j |A_ij|s_j。

- 上一次Newton更新的max|du_j|/s_j <= 1e-10。
- 非线性尺度指标eta_N=max_i |r_i|/d_i <= 32 epsilon = 7.105427357601002e-15。
- 当前线性解eta_L=max_i |r_i-(Adu)_i|/(|r_i|+sum_j|A_ij du_j|) <= 1e-12；误差累加使用C++long double。
- 对七个既有上部区域及不含Dirichlet层的free_domain，|sum_region r_salt|dt/W <= 1e-12/4800 = 2.0833333333333334e-16 m；水对应1e-10/4800 = 2.0833333333333334e-14 m。
- 最多20次更新；NaN残差/修正、零行、奇异主元、非均匀输入或以上门槛失败均不能接受。

32在测试前固定，不按某一步的失败值拟合；16的晚期敏感性得到相同轨迹。实际接受eta_N最高约1.53e-16，低于1 epsilon，未靠接近32 epsilon的宽松边界通过。eta_L允许多项消去累计的浮点误差；实测最高约1.09e-13，均小于1e-12。

eta_N是基于线性化雅可比的非线性尺度指标，不等于严格的非线性后向误差证明，也不能直接当物理解误差。区域求和允许局部抵消，因此必须同时满足逐方程eta_N、修正门槛、独立VTU库存/SCVF积分及场桥接。每步预算分配限制代数缺陷累计，但不能包办诊断浮点误差或时间/空间截断误差。

## 验证边界

9个有界试验覆盖原初值、T附近、旧失败附近、factor16/32、探针及强制失败；共{cost['frames']}帧原场SHA逐成员校验，原场压缩包在忽略的`tmp/wp8_newton_correction_20261103`，远端证据保留。

22项现有停止语义测试检查表面/内部首次阈值、事件后接受状态拒绝、删失及时间配对。新短测实际没有达到0.1；仅证明非越阈记录有效及停止实现未改变，不宣称做过本轮真实越阈实验。强制失败试验使用原物理初值，只令Newton最大迭代为0，确认零个时间步接受、初始场保留、退出码失败。

未证明诊断回放与无损重启等价；最终4T必须从原始冻结初值连续开始。未验证完整4T，不以短测替代全程物理验收。
''')
    write('FULL4T_GO_NOGO_AND_COST_ZH.md',f'''# 完整4T：条件准入、执行协议及预算

**数值准备GO；启动权限NO-GO（本轮不启动）；完整4T科学验收PENDING。** command/14要求先交付条件方案，本轮只执行短回归。`pending_full4T.input`和`full4T_protocol.json`已冻结；没有提交新的4T作业，也不启动R96、多波长、完整域或8T。

## 下一次获准执行时

1. 核对command、输入、v2源manifest和生产ELF SHA；实时检查frodo CPU、内存、磁盘，选择空闲一核，BLAS/OMP线程各1；不以旧负载快照代替检查。
2. 从原冻结初值重新开始，单格横向、200层垂向、原七观测面和所有物性不变。禁止使用本轮ReplayProfile或拼接旧3.6775T前缀。dt=113.68588647695695 s，目标545692.2550893936 s，共4800步；不改变固定dt/2协议。
3. 每个接受步检查表面或全域max(x)>=0.1，首次触阈立即停止后续接受步，记录区间；触阈停止是删失/适用性结果，不能宣称走完4T。
4. 每步保存事件、区域库存/通量/残差、约束Newton记录及原始场。T、2T、3T、4T检查点只审查已接受状态，不重置累计量或初值。最多900 s运行硬上限；失败、超时、hash异常或横向不等时保留证据并停止，不自动重试。
5. 完成后原场全部回传并逐成员SHA；独立重构phi Sw x库存、七面及free_domain收支、通量积分。与旧T及旧接受前缀仅在32ULP内严格共有时刻比较，排除被裁剪的非共有终点，不插值。前三个T检查点的预算/桥接不通过则不发布全程对照。

## 原门槛不变

{table(['量','限值'],[[k,v] for k,v in gate.items()])}

必须额外保持压力与浓度横向差精确0。最终通过必须同时满足正常退出/正确阈值退出、每步一次Newton成功、失败迭代未接受、源/ELF/input/原场SHA齐全、固定dt及原初值、全程区域预算、冻结T桥接。短测GO并不预判这些条件。

## 实测与预算

本轮9个短测共{cost['total_accepted_steps']}个接受步，{cost['frames']}帧；求解+在线检查合计{cost['solver_and_monitor_wall_s']:.3f} s，含归档{cost['total_worker_wall_s']:.3f} s（不是所有任务同时开始的日历跨度），原场压缩合计{cost['archives_bytes']/1024**2:.3f} MiB，求解子进程峰值{cost['max_solver_rss_kib']/1024:.2f} MiB。三次隔离配置/编译合计{sum(cost['compile_wall_s'].values()):.3f} s，明细和日志保留；网络等待及监测Python进程峰值内存未测，保持null。

旧4413步耗时156.7 s，新候选增加线性后向误差/接受前求解；据旧整段和本轮短测，下一次连续4T规划180–600 s（约3–10分钟），硬上限900 s/一核，即最多0.25核时。短测含0.5秒监测轮询开销，不能简单按10步墙钟线性给精确ETA。压缩原场预计约55–65 MiB，预算每端250 MiB；展开场约250 MiB，预检至少2 GiB可用磁盘；求解器地址空间硬限2 GiB。

以上预算仅单条严格背景验证，不含二维敏感性研究。所有既有失败保留，任何追加长任务仍需相应指令。
''')
    early=rows(c.ROOT/'analysis/wp8_strict_background_20261103/grid_comparison.csv');rr=[r for r in early if r['role']=='M1'][:3]
    write('SCIENTIFIC_BOUNDARY_ZH.md',f'''# 3.6775T接受前缀能支持什么

旧严格背景4414帧及其物理接受门槛仍有效；本轮仅修订数值收敛判据，未新增正式背景4T或二维物理结果。对早期均值可以形成有限而明确的数值结论：T、2T、3T处R24/R48的M1相对严格背景均为负，两个横向网格间量值差约0.31%、0.85%、1.51%。这支持在既有协议下可重复的早期零模态响应；并不构成时间、垂向与域宽误差受控的物理网格收敛证明。

{table(['时刻/T（近似）','R24平均差','R48平均差','横向网格绝对差'],[[i+1,r['R24_delta_mean_top'],r['R48_delta_mean_top'],r['strict_grid_absolute_gap']] for i,r in enumerate(rr)])}

1. 数值onset或诊断departure：依赖可观测量、阈值和连续5步要求；原1e-12/1e-10/1e-8门槛记录保留，不能选其中一个充当物理输运开始。
2. 零模态：末态M1平均表面差R24为−1.136288e-7、R48为−4.142072e-8 mol/mol；M0平均偏移较小，支持扰动相关的均值响应。
3. 守恒再分配：最浅面初始库存修正后M1相对背景减少4.393991e-12/2.145782e-12 m，与额外累计通量相容，约为当时背景库存的0.000115584%/0.0000564448%。自身库存仍增长，净通量仍向上；不能把“相对更低”写为净盐输出或自身库存下降。
4. 网格收敛：末态严格背景差的跨网格差7.220809e-8与M1−M0原差7.212121e-8相近。更严格的参照没有消除晚期离散敏感性；完整4T背景也不会自动解决二维网格误差。
5. 有效输运：尚无经过时间/垂向/域宽误差检验且达到预登记物理量级的稳健时间标签，t_eff保持null/unidentified。
6. 首次局部溶解度：M1局部峰值反而高于严格背景约1.976e-4/1.075e-4。平均降低不推出局部触阈推迟；共同窗口尚无0.1越阈，时间差不可识别，不生成结晶延迟训练标签。

不能用同一数据既后验调准门槛又宣称有效输运的外部验证。本轮解决的是背景参照的数值可运行性准备，物理主张仍受原证据边界限制。
''')
    write('README_ZH.md','''# command/14交付索引

已完成方程级诊断、隔离修订、编译和9个有界回归。数值准备GO；新的连续4T未启动，完整4T科学验收仍PENDING。

- [数值失败分析](NUMERICAL_FAILURE_ANALYSIS_ZH.md)：原失败、舍入敏感性、矩阵/差分尺度证据、回归结果。
- [源码审查与准则](SOURCE_REVIEW_AND_CONTRACT_ZH.md)：公式、变量尺度、守恒约束和证据限制。
- [完整4T条件方案与预算](FULL4T_GO_NOGO_AND_COST_ZH.md)：固定输入、GO/NO-GO、原门槛、预计成本。
- [科学解释边界](SCIENTIFIC_BOUNDARY_ZH.md)：早期零模态、再分配、网格和局部阈值的区别。
- [数值诊断图](Figure_Newton_diagnostics.png)：原失败残差、近邻回放、差分尺度和独立守恒预算；同时提供PDF。

复算：恢复archive_receipt所指的忽略缓存后，运行`validate_wp8_v2.py`、`report_wp8_v2.py`和`verify_wp8_v2_delivery.py`。这些脚本不启动模拟。`wp8_v2_regression.py`仅含有界短测名单，已有目录阻止重复prepare/launch；不提供完整4T自动启动入口。
''')
    print(json.dumps(dict(comparison=comparison,cost={k:x for k,x in cost.items() if k!='runs'},full4T_launched=False),indent=2))
if __name__=='__main__':main()
