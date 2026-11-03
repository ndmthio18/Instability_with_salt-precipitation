"""Evidence-bound WP6 first-stage report and figures; no compute submission."""
import csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from single_wavelength_post import widths
from analyze_wp6 import ROOT,OUT,rows,nums
CACHE=ROOT/'tmp/wp6_execution_20261103'
def load(p):return json.loads(Path(p).read_text())
def write(name,s):(OUT/name).write_text(s.strip()+'\n',encoding='utf-8')
def save(fig,name):
    for ext in ['png','pdf']:fig.savefig(OUT/(name+'.'+ext),dpi=200,bbox_inches='tight')
    plt.close(fig)
def main():
    b=load(OUT/'bridge_checks.json');e=load(OUT/'extend_checks.json');assert b['passed'] and e['numerical_checks_passed']
    w={r:rows(OUT/'runs'/f'WP6-{r}-extend/post/weighted_diagnostics.csv') for r in ['M0','M1']}
    on=rows(OUT/'runs/WP6-M1-extend/weighted_online.csv');t=nums(w['M1'],'time_s');a=nums(on,'weighted_fourier');h2=nums(on,'harmonic2_top');h3=nums(on,'harmonic3_top');fraction=(h2*h2+h3*h3)/(a*a+h2*h2+h3*h3)
    sig=nums(w['M1'],'weighted_top_sigma');g=sig/min(sig);m=e['metrics']['M1'];zero=e['metrics']['M0'];onset=e['onsets']['M1'];pair=e['paired_summary'][:7]
    plt.rcParams.update({'font.size':10,'axes.spines.right':False,'axes.spines.top':False,'pdf.fonttype':42})
    fig,ax=plt.subplots(2,2,figsize=(10,7),layout='constrained')
    ax[0,0].plot(t/3600,sig*1e6,color='#0072B2');ax[0,0].axvline(136423.0637723484/3600,c='grey',ls=':',label='WP5 endpoint');ax[0,0].set(ylabel='Weighted sigma (10$^{-6}$ mol/mol)',title='a  Seeded surface perturbation');ax[0,0].legend()
    ax[0,1].plot(t/3600,g,color='#0072B2');ax[0,1].set(ylabel='sigma / minimum sigma',title='b  Amplification from observed minimum')
    for role,c in [('M0','#222222'),('M1','#0072B2')]:
        ax[1,0].plot(nums(w[role],'time_s')/3600,nums(w[role],'weighted_top_mean'),c=c,label=role+' mean');ax[1,0].plot(nums(w[role],'time_s')/3600,nums(w[role],'domain_max'),'--',c=c,label=role+' max')
    ax[1,0].axhline(.1,color='#D55E00',ls=':',label='Solubility limit');ax[1,0].set(ylabel='Salt mole fraction (mol/mol)',title='c  Actual concentration; no extrapolation',ylim=(0,.105));ax[1,0].legend(fontsize=8)
    ax[1,1].semilogy(t/3600,fraction,color='#009E73');ax[1,1].axhline(.01,c='grey',ls=':',label='Review trigger, not a physical boundary');ax[1,1].set(ylabel='Energy fraction in top harmonics 2 + 3',title='d  Deviation from a single top harmonic');ax[1,1].legend(fontsize=8)
    for x in ax.flat:x.set_xlabel('Physical time (h)')
    save(fig,'Figure1_growth_solubility')
    paired=rows(OUT/'extend_paired_transport.csv');bud=rows(OUT/'runs/WP6-M1-extend/transport_budget.csv');meta=load(OUT/'runs/WP6-M1-extend/transport_metadata.json');cuts=meta['cuts'][:7]
    fig,ax=plt.subplots(2,2,figsize=(11,7),layout='constrained')
    for j,c in enumerate(cuts):
        rr=[r for r in paired if r['region']==c['name']];label=f"{200*(1-c['actual_z_over_H']):.2f} mm";color=plt.get_cmap('tab10')(j)
        ax[0,0].plot(nums(rr,'time_s')/3600,nums(rr,'delta_C_m')*1e15,label=label,color=color)
        ax[0,1].plot(nums(rr,'time_s')/3600,nums(rr,'relative_delta_C_initial_inventory'),label=label,color=color)
        bb=[r for r in bud if r['region']==c['name']]
        ax[1,0].plot(nums(bb,'time_s')/3600,nums(bb,'cut_salt_adv_down_m_s')*1e11,color=color,label=label+' adv')
        ax[1,0].plot(nums(bb,'time_s')/3600,nums(bb,'cut_salt_diff_down_m_s')*1e11,'--',color=color)
        ax[1,1].plot(nums(bb,'time_s')/3600,nums(bb,'salt_residual_m')*1e15,color=color,label=label)
    ax[0,0].set(ylabel='Paired cumulative salt difference (10$^{-15}$ m)',title='a  All seven cuts; downward positive');ax[0,0].legend(fontsize=8,ncol=2)
    ax[0,1].set(ylabel='Delta C / initial M0 regional inventory',title='b  Relative excess transport')
    ax[1,0].set(ylabel='Signed salt flux (10$^{-11}$ m/s)',title='c  M1: advection solid, diffusion dashed')
    ax[1,1].set(ylabel='Salt budget residual (10$^{-15}$ m)',title='d  M1 closure; not total error')
    for x in ax.flat:x.set_xlabel('Physical time (h)')
    fig.suptitle('Nonzero differences alone do not establish effective transport',fontsize=12)
    save(fig,'Figure2_seven_plane_transport')
    fields=np.load(CACHE/'runs/WP6-M1-extend/post/representative_fields.npz');xs=fields['xs'];z=fields['z_m'];keys=sorted([k for k in fields.files if k.startswith('x_')],key=lambda k:int(k.split('_')[1]));keys=keys[-2:];anom=[]
    for k in keys:
        x=fields[k];mean=x@widths(xs)/(xs[-1]-xs[0]);anom.append((x-mean[:,None])*1e6)
    limit=max(float(np.max(abs(x))) for x in anom);mask=(.2-z)<=.005
    fig,ax=plt.subplots(1,2,figsize=(10,4.4),layout='constrained')
    for a0,arr,k in zip(ax,anom,keys):
        im=a0.pcolormesh(xs/(xs[-1]-xs[0]),(.2-z[mask])*1000,arr[mask],cmap='RdBu_r',vmin=-limit,vmax=limit,shading='nearest',rasterized=True);a0.set(xlabel='x / wavelength',ylabel='Depth below surface (mm)',title=f"M1 at accepted step {k.split('_')[1]}");a0.invert_yaxis()
    fig.colorbar(im,ax=ax,label='x minus lateral weighted mean (10$^{-6}$ mol/mol)');fig.suptitle('Actual nodal concentration anomaly, common color scale; top 5 mm')
    save(fig,'Figure3_spatial_mode')
    newsummary=dict(end_time_s=float(t[-1]),weighted_growth_since_min=float(g[-1]),raw_growth_since_min=onset['amplification_since_min'],max_harmonic23_energy_fraction=float(max(fraction)),end_harmonic23_energy_fraction=float(fraction[-1]),odd_neumann_energy_fraction=m['odd_energy_max'],relative_sigma_end=m['relative_sigma_end'],max_shallow_relative_excess=max(abs(r['relative_delta_C_end']) for r in pair[3:]),any_persistent_candidate=any(r['persistent_candidate'] for r in pair),solubility_events={r:e['metrics'][r]['stop'] for r in ['M0','M1']})
    (OUT/'scientific_summary.json').write_text(json.dumps(newsummary,indent=2)+'\n')
    table='| 实际深度mm | 末端Delta C (m) | /初始M0库存 | 窗口最大闭合包络(m) | 分类 |\n|---|---:|---:|---:|---|\n'
    for c,r in zip(cuts,pair):table+=f"| {200*(1-c['actual_z_over_H']):.3f} | {r['delta_C_end_m']:.5e} | {r['relative_delta_C_end']:.5e} | {r['max_residual_envelope_m']:.5e} | {r['classification']} |\n"
    write('RESULTS_ZH.md',f'''# WP6首阶段科学结果

依据command/7完成新诊断、T桥接及通过后的2T配对。所有原WP5结果与Paper1未修改。本轮停止在2T；没有4T/8T、全域、多波长或nx60求解。

## 增长、结构与浓度

M0/M1均到达{t[-1]:.10f} s，各{len(t)-1}接受步。M1加权表面sigma从窗口最低点到末端放大 **{g[-1]:.6f}倍**；冻结Paper1风格原sigma协议放大{onset['amplification_since_min']:.6f}倍，提取t_detect={onset['onset_refined_s']:.6f} s。二者统计定义不同，不能把末端原节点sigma直接与其他网格相减。

末端相对扰动sigma/mean={m['relative_sigma_end']:.6e}；表面第2/3谐波能量比例最大{max(fraction):.6e}，奇Neumann模态能量比最大{m['odd_energy_max']:.6e}。实际空间场见Figure3，相同色标显示T与2T上部5 mm异常场。谱与空间证据用于判断是否出现新结构，不把非零高次谐波自动称为成熟非线性指状流。

M0末端表面均值/max为{zero['final_surface_mean']:.10f}/{zero['final_surface_max']:.10f}；M1为{m['final_surface_mean']:.10f}/{m['final_surface_max']:.10f} mol/mol。M1全域max={m['final_domain_max']:.10f}。两腿本窗均未观测实际0.1溶解度事件，事件时间右删失；不能给出两腿阈值时间差，更不能从末端浓度差推算结晶延迟。

背景M0最大端点差为{zero['max_endpoint_difference']:.6e} mol/mol，超过冻结对称性筛查的绝对下限1e-12；背景最大sigma为{zero['max_sigma']:.6e}，远小于M1信号，但不能事后提高门槛宣称所有GO条件通过。核心数值/守恒检查通过，域代表性门槛标记review_required，4T暂不准入；这不是已经发现物理长波失稳的证据。

## 七个预登记面的配对额外输运

{table}

Delta C按向下为正；已减去实测初始库存差检查配对库存恒等式。表中非零值不能单独作为物理效应：需看预登记持续门槛、预算包络、时间/网格/域误差。U_dt、U_z、U_x、U_domain均unknown，未用0代替，没有定义t_eff。各腿总通量、adv/diff有符号值及逐时归一化结果保存在原budget与`extend_paired_transport.csv`；Figure2展示全部7面，不只挑浅层。

顶/侧盐边界为零，D_I+Delta C=R1−R0。正Delta C不一定表示实际净向下输运，也可能表示背景向上迁移减弱；需同时读取单腿signed total/adv/diff。库存和通量单位来自phi Sw x归一化守恒量，并非kg或mol。自由域底部反应通量不是独立物理验证。

## 本轮证据边界

本轮是单案例、单canonical初值、单波长受限对称子空间的延长观测。T桥接证明新增诊断没有改变旧解，不等于2T输运已完成数值收敛或更宽域代表性验证。阈值未触发不是无限时间稳定结论；微小增长也不是有效盐再分配标签。
''')
    print(json.dumps(newsummary,indent=2))
if __name__=='__main__':main()
