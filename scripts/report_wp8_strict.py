"""Generate command-13 evidence reports including a failed full-horizon gate."""
import json
from wp8_strict_campaign import ROOT,OUT,T,save,sha
from analyze_wp8 import rows
def table(headers,data):return '| '+' | '.join(headers)+' |\n|'+'|'.join(['---']*len(headers))+'|\n'+'\n'.join('| '+' | '.join(map(str,r))+' |' for r in data)+'\n'
def write(name,s):(OUT/name).write_text(s.strip()+'\n',encoding='utf-8')
def main():
    gates={k:json.loads((OUT/k/'validation.json').read_text()) for k in ['T','4T']};r=json.loads((OUT/'results.json').read_text());failure=json.loads((OUT/'failure_diagnosis.json').read_text());tr=rows(OUT/'strict_transport_final.csv');summary=r['summary'];end=r['valid_background_end_s']
    costs=[]
    for stage,g in gates.items():
        c=g['cost'];costs.append(dict(stage=stage,completed=g['passed'],accepted_steps=g['events'].get('times') and len(g['events']['times'])-1 if g['passed'] else failure['accepted_steps'],accepted_end_s=T if stage=='T' else end,wall_with_monitor_s=c['solver_and_live_audit_wall_s'],wall_with_archive_s=c['total_wall_s'],archive_MiB=c['archive']['bytes']/1024**2,raw_fields_MiB=c['raw_field_bytes']/1024**2,solver_peak_RSS_MiB=c['peak_child_rss_kib']/1024,wall_cap_s=900,archive_cap_MiB=250,timeout=False,core_hours_including_archive=c['total_wall_s']/3600))
    save(OUT/'cost_audit.json',dict(runs=costs,one_core_per_run=True,sequential=True,total_core_hours_including_archive=sum(x['core_hours_including_archive'] for x in costs),wrapper_peak_memory='not_measured',new_runs=2,retries=0))
    write('EXECUTION_AND_GATES_ZH.md',f'''# WP8 严格背景分阶段执行与门槛

用户指令 `command/13.md`，审阅提交`14e15c7edb32a6a98faf7d915052590d35039ed6`。本轮只运行两条新背景：先0–T，通过后单独从同一初值运行0–4T。没有从T重启，没有新R96、多波长、完整域或8T。源文件、本构、原初值幅值、dt/2和原七面坐标均未修改。

## 预检与来源

复用`transport_wp8_background_v1`的r03可执行文件，ELF SHA256为`8cf2d5a410bf467059c3f2fb036a7ffa0f24c15c486d20e231a9e590f0195407`。本地manifest与编译source_sha256逐项一致；远端worker在启动前再次逐文件核对源码、ELF和输入SHA。原10步桥接通过，初始浓度/压力公式、物理参数变化白名单与实际七面坐标在本次再次检查。新输入仅Name和必要TEnd改变。

frodo是项目既有计算主机；各一次单核、单rank，运行前CPU空闲采样和内存/磁盘快照保存在各阶段launch.json。900 s单次上限、求解器2 GiB虚拟地址硬上限、每端250 MiB压缩原场归档上限均沿用原协议。4T计划固定T验收文件SHA，只有原T gate通过才启动。T gate文件保持冻结，没有为后续失败改写。

## T：PASS

1200步到{T:.10f} s，1201帧。压力、浓度横向差严格为0，仍按原1e-14阈值检查；原一格背景的3.88994e-14失败记录未改变。与原背景共有1201时刻，浓度剖面差4.51462e-15 mol/mol，库存差3.38813e-19 m，净盐通量差2.14629e-19 m/s。盐/水预算最大残差1.73714e-16 / 4.27242e-14 m，全部29项冻结验收通过。

## 4T：FAILED，已接受前缀可单独核查

接受4413步后，第4414步迭代失败，最后物理时刻{end:.10f} s = {end/T:.7f}T。4414帧已归档核验，未达到4T。不是溶解度停止，也不是墙钟或资源超限。

与已验收T运行有1200个严格共有时刻，所有共有VTU文件SHA相同，库存/通量轨迹逐值相同。T末点因独立T运行的终点截步与4T运行的普通步相差约2.06e-9 s而超过32 ULP配对约定，未插值或强配。

已接受前缀的每帧压力/浓度横向差均为0；盐/水残差分别<=9.92147e-16 / 6.84660e-14 m；独立盐库存重构差6.90332e-19 m，累计通量再积分差0。失败迭代没有被当作接受步：4413条Newton成功记录与4413个时间步一一对应，失败尾部没有新成功记录、时间步或VTU。

完整运行的`passed=false`保持不变。`accepted_prefix_passed=true`只表示保存的接受状态、时间、阈值标记、严格对称性、预算和来源可用于有限窗口分析，明确排除“成功退出/完成目标”的检查，**不把失败运行升级成4T通过**。初次失败验收记录另存validation_initial.json。所有新比较限制在这个前缀内。

## 溶解度与适用边界

到最后接受时刻，背景表面/全域最大浓度为0.05296016313357256 mol/mol，未达0.1。在线事件标记与独立场最大值逐帧一致；所有接受时刻和dt自洽，未隐藏接受过阈状态。停止记录为exception_or_incomplete。新背景的阈值时间只能记为数值中止导致的右删失于最后接受时刻，不能记成4T未触阈。本轮没有实际越阈事件，因此不能新宣称验证了真实越阈情况下的停止行为；该代码字节与既有停止实现保持相同。
''')
    write('FAILURE_DIAGNOSIS_ZH.md',f'''# 第4414步失败诊断（只读）

完整4T验证失败，未自动重跑、未改源码或容差。失败试探时刻为{failure['proposed_unaccepted_time_s']:.10f} s，不能作为已接受物理数据。失败步没有VTU，也没有进入累计通量。

从第2次Newton迭代开始，水残差最大1.08577e-18，低于固定2e-18内部判据；相对位移最大1.416795e-11，低于1e-10。盐残差在2.125894e-20–4.854729e-20之间停滞，最低值仍比2e-20高约6.29%，因此迭代0–20后按原实现拒绝该步。已接受步最大迭代数17，日志显示接近末段时迭代次数增多。

这与固定绝对盐残差判据遭遇舍入/残差评估下限相容；不是已证实的物理解发散，也不是浓度触阈。但小Newton位移本身不能证明失败试探解满足独立盐预算，因此不将其提升为有效状态，也不只为“跑通”把阈值提高。

下一项具体工程工作：在新版本中核查分量尺度、残差求和精度、Jacobian与舍入误差，设计有明确后向误差解释的停止准则，并以独立区域盐水预算约束验收。先重现原T与末段失败机制的有界测试，再登记新的T/4T预算与来源。原1e-14均匀性及物理桥接/守恒门槛不能改写。暂不以R96昂贵计算替代这个背景求解器问题。

逐迭代证据在newton_trace.csv及failed_step_iterations.csv，机器可读诊断在failure_diagnosis.json；这里只定位问题，没有实施修复或追加运行。
''')
    mt=table(['运行','末时刻平均差 mol/mol','末时刻局部最大差 mol/mol','共有帧'],[(k,f"{v['last']['delta_mean_top']:.8e}",f"{v['last']['delta_surface_max']:.8e}",v['shared_states']) for k,v in summary.items()])
    tt=table(['运行','实际深度mm','Delta I_c m','Delta C m','Delta I/I_b(t) %','Delta R m'],[(x['run'],f"{(1-float(x['actual_z_over_H']))*200:.6f}",f"{float(x['delta_I_corrected_m']):.7e}",f"{float(x['delta_C_m']):.7e}",f"{100*float(x['delta_I_over_current']):.7g}",f"{float(x['delta_R_m']):.3e}") for x in tr if x['region'] in ['upper_0.995000','upper_0.990000','upper_0.975000']])
    write('STRICT_RESPONSE_ZH.md',f'''# 严格背景前缀上的均值与守恒响应

有效比较截止{end:.10f} s（{end/T:.7f}T），每条二维轨迹4414个共有接受状态。未生成3.6775T之后的背景数据，未计算虚构4T对照。原WP8的4T M1−M0数值仍保留，但不能与本表的较短窗口混用。

{mt}
图A为均值差，图B分别给出表面平均值与局部最大值，图D为末态浅层深度剖面。M0平均偏移很小并不等于M0完全无横向结构：末态R24/R48 M0表面最大值高于严格背景3.07226e-6 / 4.38908e-7 mol/mol。

M1平均差在末态分别−1.13629e-7与−4.14207e-8，跨网格绝对差7.22081e-8 mol/mol；同一时刻M1−M0的跨网格差7.21212e-8，换成严格背景没有使均值响应更稳健。T–3T两网格较接近，但这不填补时间步、垂向和域宽误差。

## 七面库存和通量

qbar=phi*<Sw*x>；二维端仍使用已核验的真实乘积平均，严格背景Sw和x横向恒定，因此其协方差严格为0。库存按实际BOX体积和同一物理切面比较。Delta I_c扣除各自初始库存差，C为接受步SCVF净盐通量积分，向下离开上部区域为正，核查Delta I_c+Delta C=Delta R。库存单位m、通量m/s，不是kg盐。

{tt}
全部28个配对区域的预算恒等式最大误差{r['max_budget_identity_error_m']:.8e} m。CSV压缩文件strict_transport.csv.gz保留七面全部共有时刻的背景/二维库存、初始修正差、adv/diff/total符号、累计量、残差以及相对初始/当时背景库存的比例。

最浅层M1相对背景库存减少：R24约0.000115584%，R48约0.0000564448%。其量值高于对应闭合残差，可作为离散轨迹的可分辨再分配；尚非网格稳健或物理上有意义的净盐输出。所有四条轨迹、七个上部区域，在本共同窗口内的库存仍逐步增加；这里“减少”指相对背景较低，并不指自身随时间下降。

最浅面背景净J约−6.48162e-12 m/s，R24 M1约−6.48111e-12、R48 M1约−6.48153e-12，均仍向上。M0差量约1e-15或更小，常与预算残差同量级，不能据非零数值宣称有效输运。严格背景解决参照可用性，不自动解决小差量可分辨性。

## 偏离时间只作诊断

departure_diagnostics.csv预先采用1e-12、1e-10、1e-8 mol/mol三档最大深度平均差，要求连续5个接受状态，不作物理onset或t_eff标签。R24 M0超过1e-12的首个诊断区间(478390.2103,478503.8962] s，R48为(500672.6440,500786.3299] s；两者在有效窗口内均未持续超过1e-10。不存在一个与阈值无关的唯一“开始偏离”时刻。
''')
    ct=table(['阶段','接受步','末时刻/T','运行+在线检查s','含归档s','压缩MiB','原场MiB','求解器峰值RSS MiB'],[(x['stage'],x['accepted_steps'],f"{x['accepted_end_s']/T:.7f}",f"{x['wall_with_monitor_s']:.3f}",f"{x['wall_with_archive_s']:.3f}",f"{x['archive_MiB']:.3f}",f"{x['raw_fields_MiB']:.3f}",f"{x['solver_peak_RSS_MiB']:.3f}") for x in costs])
    write('COST_AND_PROVENANCE_ZH.md',f'''# 成本与来源

{ct}
两阶段顺序各占一核，含在线检查和压缩归档总墙钟约{sum(x['wall_with_archive_s'] for x in costs):.3f} s，合计约{sum(x['core_hours_including_archive'] for x in costs):.6f}核时；不含本地分析和网络传输等待。没有重试成本，也没有超时。4T是数值失败，不能因为成本低于预算改写其状态。

每次900 s墙钟上限、每端250 MiB**压缩原场归档**上限均满足；此上限并非展开原场大小上限。各阶段原场仍在frodo，压缩包已回传本地忽略的tmp/wp8_strict_background_20261103并逐成员核对SHA，1201+4414=5615帧全部核验。求解器实施2 GiB地址空间硬上限；表中是求解器子进程峰值RSS，监测Python进程峰值内存未单独记录，不能将未测值填0。

源和ELF未变；stage plan保留command/13与原protocol SHA、源码/可执行文件/input SHA、资源快照，以及4T启动所依赖的T gate SHA。原WP8及更早源和结果按旧provenance再核查；Paper1冻结handoff核验通过，未改Paper1。

复现后处理：先按archive_receipt与SHA恢复两个阶段的忽略缓存，然后运行validate_wp8_strict.py、analyze_wp8_strict.py、diagnose_wp8_strict.py、plot_wp8_strict.py、report_wp8_strict.py、verify_wp8_strict_delivery.py。4T验收预期返回失败；这是被保存的科学结果，不应为了让脚本退出0重跑模拟。已有run目录和launch记录会阻止重复提交。
''')
    write('SCIENTIFIC_REVIEW_ZH.md',f'''# 科学复核与下一步

**T验证通过；完整4T验证失败。** 可用严格背景前缀到{end:.3f} s={end/T:.4f}T；未自动修改求解器、放宽门槛或重跑。

1. **是否严格均匀？** 截至最后接受状态，4414帧压力和浓度横向差都为0。不能说已经保持到4T，因为没有到达4T。
2. **能否复现早期背景？** 可以。T窗口与旧一格背景的剖面差4.515e-15、净盐通量差2.147e-19 m/s。独立4T运行前T的1200个严格共有场SHA及预算完全相同；旧背景原1e-14均匀性失败记录仍保留。
3. **M0何时偏离？** 以最大深度平均差>1e-12且持续5步计，R24/R48首次分别在约478504 s（3.5075T）/500786 s（3.6708T）；以1e-10计本窗口均未识别。偏离时间取决于观测量/阈值，不是新的物理onset。
4. **M1是否有可测零模态响应？** 有。共同末态平均表面差−1.1363e-7/−4.1421e-8 mol/mol，明显大于各自M0平均差，但仍是存在离散误差的数值证据。
5. **是否比M1−M0更网格稳健？** 没有。严格背景差的跨网格差7.2208e-8，原M1−M0差的跨网格差7.2121e-8，几乎没有改善。晚期网格敏感性仍在。
6. **近表面盐库存是否真的减少？** 相对严格背景，最浅层M1库存低约4.394e-12/2.146e-12 m，预算和累计通量相容；相对当时库存仅约0.0001156%/0.0000564%。自身库存仍随时间增加，净盐通量仍向上，不支持“已形成显著向下净输运”。
7. **是否影响局部最大浓度？** 有影响且方向不同于平均响应：M1表面最大值反而高于背景约1.976e-4/1.075e-4 mol/mol。均值降低不能推断局部触阈推迟。背景及当前共同窗口无首次0.1越阈，阈值时间差不可识别。
8. **还需哪些验证？** 首先诊断并验证背景求解器的晚期盐残差停止判据，再按独立来源和预算重做完整T/4T门槛；之后R48–R96、晚期dt/2、垂向加密及域宽验证仍必要。2T既有时间步证据不能外推到该晚期窗口。本轮没有启动这些计算。

t_eff保持unidentified/null；不生成代理训练标签，不宣称结晶延迟。失败前缀的分析不替代完整4T验证，也不覆盖已有R24/R48晚期不收敛结论。
''')
    print('five reports and cost audit written; full 4T gate remains FAILED')
if __name__=='__main__':main()
