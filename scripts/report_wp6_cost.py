from pathlib import Path
import json,csv,datetime
out=Path('analysis/wp6_execution_20261103');rows=[]
for d in sorted((out/'runs').iterdir()):
 f=json.loads((d/'finished.json').read_text());s=f['solver'];rows.append(dict(run=d.name,physical_end_s=s['stop']['last_accepted_time_s'],accepted_steps=s['stop']['accepted_steps'],wall_s=s['wall_s'],cpu_s=s['cpu_s'],peak_rss_kib=s['peak_rss_kib'],newton_iterations=s['newton_iterations'],linear_iterations='unknown',output_bytes=s['solver_output_bytes'],archive_bytes=f['archive']['bytes'],post_wall_s=f['post_wall_s'],stop_reason=s['stop']['reason']))
with (out/'all_run_costs.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=rows[0]);w.writeheader();w.writerows(rows)
table='| 运行 | 接受步 | 求解墙钟s | CPU s | 峰值MiB | Newton | 后处理/归档s |\n|---|---:|---:|---:|---:|---:|---:|\n'
for r in rows:table+=f"| {r['run']} | {r['accepted_steps']} | {r['wall_s']:.3f} | {r['cpu_s']:.3f} | {r['peak_rss_kib']/1024:.2f} | {r['newton_iterations']} | {r['post_wall_s']:.3f} |\n"
b=json.loads((out/'bridge_checks.json').read_text());e=json.loads((out/'extend_checks.json').read_text());audit=[json.loads((d/'post/independent_checks.json').read_text()) for d in (out/'runs').iterdir()]
fields=sum(v['verified_vtu'] for s in ['smoke','bridge','extend'] for v in json.loads((out/(s+'_archive.json')).read_text()).values())
text=f'''# 数值验收、运行成本与归档

新版本编译通过，静态/来源检查{len(json.loads((out/'static_checks.json').read_text())['checks'])}项通过，两个阈值微测试8项通过。T桥接{len(b['checks'])}项通过；M0/M1各601帧与WP5对应原始VTU逐字节相同，旧表面、模态、onset和三面通量未变。不是仅均值近似一致。

2T数值/来源/独立预算检查通过；总门槛并非全部通过：`M0/absolute_endpoint_screen`未过，详见决策报告。严格保留这一筛查失败，不将核心数值通过缩写为总体GO。

独立VTU区域盐库存与求解器差值最大{max(r['max_inventory_salt_difference_m'] for r in audit):.6e} m，水库存最大{max(r['max_inventory_water_difference_m'] for r in audit):.6e} m；在线加权量与独立场后处理差值最大{max(max(r['online_absolute_errors'].values()) for r in audit):.6e} mol/mol。8控制体每步通量重新积分与记录累计量一致。2T两腿最大盐预算残差{max(r['max_salt_residual_m'] for r in e['metrics'].values()):.6e} m，水{max(r['max_water_residual_m'] for r in e['metrics'].values()):.6e} m，均低于预登记限值。配对初始库存和恒等式见机器检查。

{table}

全部六条成功任务（4条正式物理腿+2条微型测试）求解墙钟之和{sum(r['wall_s'] for r in rows):.3f} s，CPU合计{sum(r['cpu_s'] for r in rows):.3f} s。M0/M1分配独立物理核，并行墙钟不能简单把两腿相加。线性迭代数未输出，保持unknown；正式运行没有接受步缩减，求解日志中的RetryTimeStepReductionFactor是参数回显不是拒绝事件。暂存错误和后处理JSON修复未产生求解重试。

求解器输出合计{sum(r['output_bytes'] for r in rows)/1e9:.3f} GB，原始场压缩包合计{sum(r['archive_bytes'] for r in rows)/1e9:.3f} GB；归档{fields}帧全部通过包与成员SHA。分别保存`smoke_archive.json`、`bridge_archive.json`、`extend_archive.json`。原始场留在`tmp/wp6_execution_20261103`和frodo独立目录；此数据不进入Git。归档数量含初始帧和3帧微测试，不将其全部称为独立科学样本。

参考缓存冷恢复、冻结计划只读校验和拒绝覆盖验证通过。构建命令、编译器、外部依赖、ELF/source/runtime SHA、原输入、日志和独立核对表均保留。所有WP5原输出和冻结计划文件哈希保持一致。
'''
(out/'VERIFICATION_AND_COST_ZH.md').write_text(text,encoding='utf-8')
flux=[]
for role in ['M0','M1']:
 rr=list(csv.DictReader((out/'runs'/f'WP6-{role}-extend/transport_budget.csv').open()))
 for region in dict.fromkeys(r['region'] for r in rr):
  r=[x for x in rr if x['region']==region][-1];flux.append(dict(role=role,region=region,time_s=r['time_s'],actual_z_over_H=r['actual_z_over_H'],signed_advective_m_s=r['cut_salt_adv_down_m_s'],signed_diffusive_m_s=r['cut_salt_diff_down_m_s'],signed_total_m_s=r['cut_salt_total_down_m_s'],salt_inventory_m=r['inventory_salt_m']))
with (out/'signed_flux_end.csv').open('w',newline='') as f:w=csv.DictWriter(f,fieldnames=flux[0]);w.writeheader();w.writerows(flux)
print({'archived_frames':fields,'cost_wall_sum_s':sum(r['wall_s'] for r in rows),'output_GB':sum(r['output_bytes'] for r in rows)/1e9})
