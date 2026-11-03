"""Regenerate WP5 reports and figures from collected evidence; never launch a solver."""
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator, FormatStrFormatter
from analyze_single_wavelength import OUT, CACHE, ROOT, BASE, loop_time
from analyze_post_onset import readcsv, writecsv

def load(p): return json.loads(Path(p).read_text(encoding='utf-8'))
def col(rows,key): return np.array([float(r[key]) for r in rows])
def report(name,body): (OUT/(name+'_ZH.md')).write_text(body.strip()+'\n',encoding='utf-8')
def save(fig,name):
    fig.savefig(OUT/(name+'.png'),dpi=200,bbox_inches='tight')
    fig.savefig(OUT/(name+'.pdf'),bbox_inches='tight')
    plt.close(fig)

def main():
    short=load(OUT/'short_validation.json'); full=load(OUT/'full_validation.json')
    assert short['advance_to_full'] and full['numerical_passed']
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42})
    series=[('Full 600 x 200',OUT/'reference/M1/weighted_diagnostics.csv','#222222','o-'),('R24, full window',OUT/'runs/R24-M1-full/post/weighted_diagnostics.csv','#0072B2','-'),('R48, short window',OUT/'runs/R48-M1-short/post/weighted_diagnostics.csv','#D55E00','--')]
    fig,axs=plt.subplots(1,3,figsize=(12,3.7),layout='constrained')
    for label,path,c,sty in series:
        r=readcsv(path)
        for ax,key,title in zip(axs,['weighted_top_sigma','weighted_fourier','weighted_projection'],['a  Weighted surface standard deviation','b  Physical-wavenumber amplitude','c  Normalized eigenmode projection']):
            ax.plot(col(r,'time_s')/3600,col(r,key)*1e6,sty,color=c,label=label,markersize=3,lw=1.5)
            ax.set(xlabel='Physical time (h)',ylabel='Concentration amplitude (10$^{-6}$ mol/mol)',title=title)
    axs[0].legend(fontsize=8);fig.suptitle('Case 36877: identical physical mode; full reference has 13 saved fields',fontsize=12)
    save(fig,'Figure1_growth')
    fig,axs=plt.subplots(2,3,figsize=(12,7),layout='constrained')
    for label,where,c,sty in [('Full',CACHE/'reference/M1/post/profiles.npz','#222222','o-'),('R24',CACHE/'runs/R24-M1-full/post/profiles.npz','#0072B2','-')]:
        p=np.load(where);z=p['z_m']/.2
        axs[0,0].plot(p['x_mean'][-1],z,sty,color=c,label=label,markevery=15,markersize=3)
        axs[0,1].plot(p['x_sigma'][-1]*1e6,z,sty,color=c,label=label,markevery=15,markersize=3)
    axs[0,0].set(xlabel='Mean x (mol/mol)',ylabel='z/H',title='a  Final vertical mean');axs[0,0].legend()
    axs[0,0].xaxis.set_major_locator(MaxNLocator(4));axs[0,0].xaxis.set_major_formatter(FormatStrFormatter('%.3f'))
    axs[0,1].set(xlabel='Weighted sigma (10$^{-6}$ mol/mol)',ylabel='z/H',title='b  Final perturbation profile')
    region='upper_0.500000'
    for label,path,c,sty in [('Full',BASE/'M1/transport_budget.csv','#222222','--'),('R24',OUT/'runs/R24-M1-full/transport_budget.csv','#0072B2','-')]:
        r=[x for x in readcsv(path) if x['region']==region];t=col(r,'time_s')/3600
        axs[0,2].plot(t,col(r,'inventory_salt_m')*1e5,sty,color=c,label=label)
        axs[1,0].plot(t,col(r,'cut_salt_total_down_m_s')*1e10,sty,color=c,label=label)
    axs[0,2].set(xlabel='Physical time (h)',ylabel='Normalized salt inventory (10$^{-5}$ m)',title='c  Upper-region inventory (z/H > 0.5005)')
    axs[1,0].set(xlabel='Physical time (h)',ylabel='Downward normalized flux (10$^{-10}$ m/s)',title='d  Conservative cut flux (z/H = 0.5005)')
    r=[x for x in readcsv(OUT/'full_paired_transport.csv') if x['region']==region];t=col(r,'time_s')/3600
    axs[1,1].plot(t,col(r,'reference_delta_C_m')*1e17,'--',color='#222222',label='Full pair')
    axs[1,1].plot(t,col(r,'delta_C_m')*1e17,color='#0072B2',label='R24 pair')
    axs[1,1].plot(t,col(r,'residual_envelope_m')*1e17,':',color='#D55E00',label='R24 |R0| + |R1|')
    axs[1,1].set(xlabel='Physical time (h)',ylabel='Normalized cumulative salt (10$^{-17}$ m)',title='e  Extra transport: UNRESOLVED');axs[1,1].legend(fontsize=8)
    axs[1,1].set_ylim(-3,6)
    axs[1,1].text(.03,.97,'Floor: 1000 axis units\nTime / z errors unknown',transform=axs[1,1].transAxes,va='top',fontsize=8)
    for reg,c in zip(['upper_0.900000','upper_0.750000','upper_0.500000'],['#0072B2','#D55E00','#009E73']):
        r=[x for x in readcsv(OUT/'runs/R24-M1-full/transport_budget.csv') if x['region']==reg]
        axs[1,2].plot(col(r,'time_s')/3600,col(r,'salt_residual_m')*1e17,color=c,label=reg.replace('upper_','z/H ~ '))
    axs[1,2].set(xlabel='Physical time (h)',ylabel='Salt closure residual (10$^{-17}$ m)',title='f  Three control-volume budgets');axs[1,2].legend(fontsize=8)
    fig.suptitle('Mean transport agrees with the fixed-grid reference; extra transport remains below resolution',fontsize=12)
    save(fig,'Figure2_transport')
    fig,axs=plt.subplots(1,3,figsize=(12,3.9),layout='constrained');xx=np.arange(3)
    for j,(role,c) in enumerate([('M0','#0072B2'),('M1','#D55E00')]):
        vals=[loop_time((OUT/'runs'/f'{m}-{role}-{s}/solver.log').read_text(),10) for m,s in [('F600','bench'),('R24','short'),('R48','short')]]
        bars=axs[0].bar(xx+(j-.5)*.34,vals,.34,color=c,label=role)
        axs[0].bar_label(bars,fmt='%.2f',fontsize=8)
    axs[0].set(xticks=xx,xticklabels=['Full','R24','R48'],ylabel='First 10 accepted-step loop time (s)',title='a  Same frodo host, one physical core/run',ylim=(0,76));axs[0].legend()
    vals=[0,short['runs']['R24-M1-short']['metrics']['weighted_fourier_max_relative']*100,short['runs']['R48-M1-short']['metrics']['weighted_fourier_max_relative']*100]
    bars=axs[1].bar(xx,vals,color=['#666666','#0072B2','#D55E00']);axs[1].bar_label(bars,fmt='%.3f%%',fontsize=9)
    axs[1].set(xticks=xx,xticklabels=['Full reference','R24','R48'],ylabel='Max weighted Fourier difference (%)',title='b  Common short-window saved times',ylim=(0,1.05));axs[1].text(.02,.93,'Relative to fixed 600 x 200 grid;\nnot continuum error',transform=axs[1].transAxes,va='top',fontsize=9)
    vals=[2*601*201,2*25*201,2*49*201];bars=axs[2].bar(xx,np.array(vals)/1000,color=['#666666','#0072B2','#D55E00']);axs[2].bar_label(bars,fmt='%.2f',fontsize=9)
    axs[2].set(xticks=xx,xticklabels=['Full','R24','R48'],ylabel='Scalar unknowns (thousands)',title='c  Two unknowns per BOX vertex',ylim=(0,280))
    fig.suptitle('VTU every accepted step; loop timing excludes startup; one observation per configuration',fontsize=11)
    save(fig,'Figure3_cost_accuracy')
    performance=[]
    for d in sorted((OUT/'runs').iterdir()):
        p=load(d/'plan.json');f=load(d/'finished.json');s=f['solver'];log=(d/'solver.log').read_text()
        performance.append(dict(run=d.name,host='frodo',cells_x=p['nx_cells'],cells_z=p['nz_cells'],vertices=p['DOF_vertices'],scalar_unknowns=p['scalar_unknowns'],physical_end_s=s['stop']['last_accepted_time_s'],accepted_steps=s['stop']['accepted_steps'],wall_s=s['wall_s'],cpu_s=s['cpu_s'],peak_rss_kib=s['peak_rss_kib'],newton_iterations=s['newton_iterations'],linear_iterations='',linear_status='not emitted; unknown',logged_failure_retry_events=sum(not line.lstrip().startswith('--') for line in s['failure_or_retry_log_lines']),solver_output_bytes=s['solver_output_bytes'],archive_bytes=f['archive']['bytes'],post_archive_wall_s=f['post_wall_s'],first10_loop_s=loop_time(log,10)))
    writecsv(OUT/'performance_all_runs.csv',performance)
    table='| Run | Wall s | CPU s | Peak MiB | Accepted | Newton | Output MB |\n|---|---:|---:|---:|---:|---:|---:|\n'
    for r in performance: table+=f"| {r['run']} | {r['wall_s']:.3f} | {r['cpu_s']:.3f} | {r['peak_rss_kib']/1024:.2f} | {r['accepted_steps']} | {r['newton_iterations']} | {r['solver_output_bytes']/1e6:.2f} |\n"
    report('MODEL_CONFIGURATION','''# 模型配置与执行来源

原 full-domain 数值参照为 Case 36877，600×200、W=0.6000000000000001 m、H=0.2 m、nx=28；R24/R48 分别24×200、48×200，Wr=W/28=0.021428571428571432 m、nx=1。物理波数均293.2153143350473 m⁻¹。原垂向分级、物性、底部压力、蒸发、浓度阈值、幅值1e-6和本征剖面不变。M0为Constant，M1为Eigenprofile。

复用 `model/transport_v1`，未修改或编译C++。ELF SHA256为 `d78ce5122cf69e4dee5a4d12bcc726779959f304b343246009cb45bb37e894cc`。输入与剖面SHA逐条存于 `runs/*/plan.json`；原始参照和源码固定于 `reference_hashes.json`。原始600×200结果未覆盖。

所有新运行从零启动，dt=227.3717729539139 s；短窗300步/68211.5318861742 s，完整窗600步/136423.0637723484 s。输出间隔由50改为1以保存每接受步场，不改变积分设置。完整参照仅13帧可重算加权场诊断，短窗共同7帧；守恒预算和原始表面记录均逐步比较。时间仅校正PVD有限有效位到原接受时刻，匹配容差1e-7 s，不进行时间插值。不同网格节点不直接相减。

## 执行修订和失败保留

预登记 `preregistration.json` 保留运行前原文，其中solver_executed=false及sal是当时计划，不代表最终状态。sal无空闲物理核，数组19194在未启动任何求解前取消；证据见 `sal_cancelled_before_start.json`。成功运行迁移至frodo，2条独立任务各一物理核（0/1，不重叠），BLAS线程1，每条12h上限。

首次frodo r01四条启动在MPI初始化处失败，无接受状态；原因是运行库包含orted.real却缺路径适配的orted包装器。保持原ELF/物理输入，添加包装器后在独立r02重新执行，失败目录完整保留。一次r02文件暂存语法错误使管理器文件缺失，没有启动求解器；相应日志/claim另存，未覆盖。见 `launcher_failure_diagnosis.json`、`staging_failure_retained.json`、`failed_initial_launcher/`。这些是启动环境故障，不能混计为数值不收敛或删除。

成功目录 `/temp/zhixin/paper2_wp5_20261103_r02`。先完成4条短窗，登记2条同机全宽10步计时控制，再通过精度、加密和性能门槛后，仅执行R24完整窗M0/M1。总8条成功运行、2420接受步、2428帧场；原始压缩包回传至本地忽略目录 `tmp/single_wavelength_validation_20261103`，全部成员SHA已验证，远端保留。没有2T/3T、幅值修改或其他物理案例。

## 统计约定

独立 `single_wavelength_post.py` 用横向梯形长度权重求均值、方差和标准差；Fourier幅值为2/W乘去均值场与exp(-ikx)的加权积分绝对值。本征投影为去横向均值浓度在A(z)cos(kx)上的体积加权投影，除以该基函数平方积分，不额外乘Sw。盐库存为单位宽度∫phi Sw x dV，水库存为∫phi Sw dV；盐量单位是归一化长度m，并非kg或mol。

通量来自求解器实际守恒SCVF核，未将VTK速度代用为守恒通量。逐步记录局部最大值、浓度阈值和接受步预算。原始sigmatop及冻结onset提取器保持不变；加权sigma只服务网格一致比较。初始库存按积分比较且配对减去初始库存差。节点算术统计的解析采样偏差见 `initial_sampling_bias.csv`。
''')
    report('SHORT_WINDOW_VALIDATION','''# 短窗口验证

4条从零运行全部达到68211.5318861742 s，各300接受步，通过全部预登记检查。配置/解析检查49项通过；科学回归测试9项通过。机器明细见 `short_validation.json`、`short_differences.csv`、`short_paired_transport.csv`。

| M1诊断：相对固定全宽参照的最大差异 | R24 | R48 |
|---|---:|---:|
| 加权表面sigma / Fourier | 0.20011% | 0.78858% |
| 归一化本征投影 | 0.10609% | 0.41917% |
| 深度sigma剖面加权L2 | 0.13559% | 0.53512% |

R24→R48逐接受步加密差异：sigma/Fourier最大0.59314%，本征投影0.31440%，均低于预登记3%。M1加权均值差异9.23e-13（相对），库存1.99e-13（相对）。全宽参照每波长21.43单元，R24更接近其离散尺度；R48更细但离该参照略远，不能解释为R48连续模型精度更低。

背景横向噪声、端点对称、奇Neumann模态、三深度守恒通量和盐水闭合均通过。配对初始库存差低于1e-15 m，额外累计通量仍unresolved。短窗onset本次可解析，但没有把短记录中不可解析onset预设为失败条件。

短窗精度和同机性能门槛通过后选择R24进入完整窗口，未执行R48完整窗口。短窗通过不代表显著非线性输运得到验证。
''')
    report('NUMERICAL_ACCURACY','''# 精度与可分辨性

运行前 `preregistration.json` 冻结工程门槛：均值/库存/局部最大值相对1e-6；相对全宽sigma、Fourier、投影5%；R24/R48模态3%；深度sigma加权L2为5%；盐残差1e-14 m、水1e-11 m；瞬时切面通量1e-16 m/s + 1e-5相对，累计通量1e-12 m + 1e-5相对；初始库存差1e-15 m；奇Neumann能量比1e-6；端点差1e-9；完整窗onset差5%；加速至少2倍。这些是预先指定的工程筛选容差，不是通用物理阈值。

R24完整窗M1相对固定600×200参照：加权sigma/Fourier最大0.24413%，本征投影0.14897%，深度sigma剖面L2为0.17200%，加权表面均值4.649e-13、垂向均值剖面5.797e-13（均相对差异）；局部表面最大值1.681e-7（相对）。全部通过。

两腿与三切面的最大单条瞬时通量差2.6501e-20 m/s，累计通量差2.8214e-17 m，区域库存绝对差1.1960e-16 m。盐预算最大残差3.1852e-17 m，水1.0138e-14 m。M1端点差最大5.558e-13，奇Neumann能量比3.426e-14；背景近零谱不作相对能量解释，记录null。

必须分开三类结论：A，单条平均输运与固定网格参照一致；B，配对额外累计输运R24与全宽差最大3.3965e-17 m；C，R24额外输运自身最大2.9891e-17 m，与配对闭合包络最大4.0909e-17 m同阶，且低于预登记绝对筛选下限1e-14 m。因此C为unresolved，B不计算近零分母相对误差。闭合包络不是置信区间，更不是全部离散误差。

配对恒等式为 (I1-I0)-(I1(0)-I0(0))+Delta C=R1-R0；最大初始库存差6.7763e-19 m，恒等式误差小于1e-18 m。图2显示这一未分辨信号，不能把它画成有效输运。U_dt、U_z仍unknown；没有时间或垂向网格收敛证明，没有定义t_eff。横向两级测试只能检验本窗口敏感性，不建立全连续模型收敛阶。

初始节点sigma受端点重复权重影响：R24为0.72e-6，加权值为1e-6/sqrt(2)=0.70710678e-6；旧全宽节点值约0.70769286e-6。不能将原节点sigma直接相减视为物理误差。onset仍按旧协议提取，结论是该检测协议在本案例的再现性，不是LSA临界时刻或输运有效时刻的验证。
''')
    report('PERFORMANCE_BENCHMARK','# 性能实测\n\n'+table+'''
所有成功计算在frodo各一物理核，独立任务同时最多2条；无跨机器耗时相除。完整自由度与档案尺寸见 `performance_all_runs.csv`。所有线性求解器迭代数为unknown，因为原ELF不输出，未填0。日志中的Newton.RetryTimeStepReductionFactor是参数回显，不是拒绝步；成功运行没有记录到求解失败/重试，接受时间步保持指定dt。原始启动失败另见模型配置报告。

## 公平短时计时控制

F600-M0/M1各从零只跑10步，是为满足同机性能比较额外登记的有界控制，非新增物理案例。与R24/R48对应前10步同机、同ELF、同输出间隔1、同物理时刻。求解循环时间排除启动；不含后处理/归档，后者在表格CSV单列。M0：60.637 s / 1.4592 s = 41.555倍（R24），/3.5467 s =17.097倍（R48）；M1：64.653 s /1.5101 s =42.814倍（R24），/3.649 s =17.718倍（R48）。各只有一次测量，不提供虚构置信区间或稳态全窗加速估计。

R24完整窗M0/M1墙钟81.224/81.675 s，另后处理与归档9.039/10.031 s；逐步输出约530/553 MB。全宽10步峰值约435 MiB，R24完整窗约44 MiB，不能据此保证全宽长窗峰值。

历史同frodo全宽M0墙钟2492.687 s，对本次R24完整窗比为30.689倍，仅为历史端到端参照：旧输出间隔50、本次1，CPU亲和性、负载和测时包装不同。来源 `analysis/unperturbed_36877_20261102/evidence/finished.json`。完整M1没有同机全窗控制，故其全窗加速比未测量；不能拿sal约9523 s除以frodo81.7 s。约42倍只属于本表前10步口径。
''')
    report('FULL_WINDOW_VALIDATION','''# 完整基准窗口验证

R24 M0/M1从t=0重新执行到136423.0637723484 s，各600接受步、601帧；无checkpoint拼接。全部完整窗门槛通过。结果见 `full_validation.json`、`full_differences.csv`、`full_paired_transport.csv` 及 Figures 1–3。

沿用冻结的Paper1风格原始sigma提取器，M1数值检测时刻t_detect=50632.5621 s，完整域参照50394.3336 s，相差238.2285 s（0.47273%）。R24的plateau区间44337.4957–56842.9432 s，与完整域区间重叠。这一宽plateau不应写成高精度物理onset。R24端点/最小sigma比1.186756，仍为微小扰动增长；线性onset可再现不等于已观察强非线性转捩。

全窗weighted sigma/Fourier差异上限0.24413%，本征投影0.14897%，剖面sigma L2为0.17200%；背景最大加权sigma约1.585e-13。平均浓度、整体及区域库存、三控制面守恒通量、阈值记录和闭合均通过。没有达到溶解度停止条件，阈值时间为本窗右删失；不能推断未来不会达到。

配对额外累计输运仍unresolved（最大2.989e-17 m），未建立t_eff。所有场检查密度、饱和度、输入初始场及原始节点sigma一致性。缩域的对称性在本窗口内保持，并不排除更长时间的对称性破缺。

选择R24作为本案例至T的低成本已验证方案。本轮没有执行R48全窗、2T/3T、其他物理案例或Phase2。
''')
    report('VALIDITY_LIMITATIONS','''# 适用范围和回退条件

| 研究阶段 | 本轮证据 | 判断 |
|---|---|---|
| 早期单模态衰减、检测及随后的微小增长 | 全窗与固定全域参照一致 | Case36877条件性GO |
| 弱非线性单模态发展 | 数学对称子空间允许，但本窗扰动仍微小 | 可作后续受控筛选；未独立验证显著非线性 |
| 显著对流盐输运 | 配对额外输运低于可分辨量级 | inconclusive，不能宣称替代验证通过 |
| 多模态竞争、跨波长合并 | 缩域排除大量允许模式 | 必须完整域或另行验证的多波长域 |

侧面Neumann不是周期边界。精确初始化在反射对称子空间内，端点值相等才可按Wr重复；奇Neumann模态、端点差及镜像误差都要继续监测。即使缩域保持对称，完整域仍可能对排除模态不稳定，因此不能只靠缩域自身无异常来确认充分性。

后续full-domain fallback触发项：多模态增长或相位竞争、横向非均匀土壤/边界、相邻指状流不等价/合并、奇模态或端点差超过登记容差、全域抽查和缩域结果显著分离、阈值时刻对网格/时间步/窗口不稳、额外输运接近误差包络。回退时保持物理参数/扰动协议，逐项做全域对照。这些是工程适用性门槛，非通用物理阈值。

不能把R24直接用于其他物理案例。原600单元下nx=3对应200单元/波长（case25748），nx=2对应300单元/波长（case27710）；强行使用24单元会显著粗化。应从约R200/R300的相当分辨率起步，再做受控加密；成本优势也不能沿用nx28的约42倍。nx=1没有域缩减。未来还须维持相关土壤参数联合支撑、闭合分支和初始化来源，不能为方便而替换物理波数。

全宽参照自身不是连续极限真值，垂向和时间误差未知；未模拟溶解度后的结晶动力学。本轮没有新的实验验证，不能把数值等价扩展成土壤盐壳预测能力。
''')
    report('NEXT_STAGE_GO_NO_GO','''# WP5 决策：本案例早期筛选条件性GO，普适非线性替代尚未验证

1. **数学依据：有条件成立。** 同物理波数、均质横向系数/边界、cos初始化形成反射对称子空间；Neumann不等于periodic，完整扰动空间并未保留。
2. **onset及后续增长：本窗口可再现。** t_detect差0.47273%，weighted sigma/Fourier最大差0.24413%。验证的是原检测协议和微小单模态增长，不是新的LSA临界时间。
3. **平均盐量与输运：单条量通过，额外输运未分辨。** 均值/库存/三切面通量/预算通过；M1−M0累计差仍约1e-17 m，不能定义t_eff。
4. **分辨率选择：本案例R24足够。** R24→R48短窗模态变化最大0.59314%，R24完整窗全部门槛通过，无证据要求再跑R48完整窗。其他案例必须按物理波长保留相当分辨率。
5. **实测加速：同frodo前10步R24为41.55–42.81倍，R48为17.10–17.72倍。** R24完整窗每腿约81–82 s；同机M1全窗加速未测量，不能用跨主机时间比代替。
6. **更长非线性筛选：方法上可设计，尚非已验证替代。** 下一轮应预登记有界窗口、时间/垂向误差检查、全域抽查和停止规则；本轮没有启动2T/3T。
7. **必须恢复完整域：** 横向异质性、多模态竞争、跨波长合并、对称性破缺、阈值预测不稳或缩域与全域明显分离；缩域保持对称本身不足以证明全域稳定。
8. **三种instability–transport–solubility区域：设计准备GO，自动批量研究NO-GO。** 当前只有一个物理案例、微小增长和未分辨额外输运，无法确立三种区域。先选有来源保障的代表案例，按波长保持网格尺度，做有界延窗与全域控制，并实测输运误差；不启动Phase2或机器学习。

## 本轮交付与复现

8份报告、3组PNG/PDF、逐步/共有时刻CSV、全部输入/日志/哈希和8条成功运行证据已保存。4条初始MPI启动失败及1次暂存失败保留。原始2428帧归档于忽略目录，SHA记录于 `short_archive.json`、`bench_archive.json`、`full_archive.json`，远端r02保留。原Paper1与完整域基准未改。

只做本地后处理的复现入口：`python scripts/analyze_single_wavelength.py short`，`python scripts/analyze_single_wavelength.py full`，`python scripts/test_single_wavelength.py`，`python scripts/report_single_wavelength.py`。需先按归档清单恢复忽略的profiles/原始场；不要为复现报告调用stage/launch。`preregistration.json`是不可回写的运行前快照，最终状态由finished和验证JSON给出。

下一项具体工作是形成分辨率匹配的有界非线性筛选方案及完整域回退抽查预算。科学结论保持：成本障碍显著降低，但有效额外盐输运尚未被本案例证据识别。
''')
    print('Wrote seven reports, three PNG/PDF figure groups, and eight-run performance CSV.')

if __name__=='__main__': main()
