"""Collect live/final evidence and report for the one no-imposed-perturbation test."""
import base64,csv,hashlib,io,json,subprocess,math
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'analysis/unperturbed_36877_20261102'
launch=json.loads((OUT/'launch.json').read_text());plan=json.loads((OUT/'plan.json').read_text())
code='''import base64,json,hashlib,csv,io
from pathlib import Path
r=Path(REMOTE);files={}
names=['stop_record.json','surface_events.csv','transport_budget.csv','transport_metadata.json','finished.json','worker.log','solver_pid.json']
for name in names:
 p=r/name
 if p.exists():
  b=p.read_bytes()
  if name.endswith('.csv') and b and not b.endswith(b'\\n'):b=b[:b.rfind(b'\\n')+1]
  files[name]={'sha256':hashlib.sha256(b).hexdigest(),'base64':base64.b64encode(b).decode()}
log=r/'solver.log'
if log.exists():
 with log.open('rb') as f:
  f.seek(max(0,log.stat().st_size-12000));tail=f.read().decode(errors='replace')
else:tail=''
print(json.dumps({'files':files,'log_tail':tail}))
'''.replace('REMOTE',repr(launch['root']))
bridge='/mnt/d/3OutputFile/1Research/Article/22Surrogate_model/paper/2/scripts/frodo_bridge.py'
p=subprocess.run(['wsl','python3',bridge,'run','python3 -'],input=code.encode(),capture_output=True,timeout=60)
if p.returncode:raise RuntimeError(p.stderr.decode())
j=json.loads(p.stdout);dest=OUT/'evidence';dest.mkdir(exist_ok=True)
for name,rec in j['files'].items():
 b=base64.b64decode(rec['base64']);assert hashlib.sha256(b).hexdigest()==rec['sha256'];(dest/name).write_bytes(b)
(dest/'solver_tail.log').write_text(j['log_tail'],encoding='utf-8')
events=list(csv.DictReader((dest/'surface_events.csv').open()))
if not events:raise SystemExit('Waiting for first complete diagnostic row')
stop=json.loads((dest/'stop_record.json').read_text());finished=json.loads((dest/'finished.json').read_text()) if (dest/'finished.json').exists() else None
baseline=max(1e-16,float(events[0]['sigma_top_node']))
threshold=max(1e-8,1000*baseline);onset=None;streak=0
for row in events:
 if float(row['sigma_top_node'])>threshold:streak+=1
 else:streak=0
 if streak==5 and onset is None:onset=float(events[events.index(row)-4]['time_s'])
last=events[-1]
bud=list(csv.DictReader((dest/'transport_budget.csv').open()));initial={r['region']:r for r in bud if r['step']=='0'}
res={s:max(abs(float(r[s+'_residual_m']))/float(initial[r['region']]['inventory_'+s+'_m']) for r in bud) for s in ['salt','water']}
metrics={'completed':finished is not None,'returncode':finished['returncode'] if finished else None,'stop':stop,'last_time_s':float(last['time_s']),'max_sigma':max(float(r['sigma_top_node']) for r in events),'last_sigma':float(last['sigma_top_node']),'last_top_max':float(last['x_top_max']),'last_top_mean':float(last['x_top_area_mean']),'visibility_threshold':threshold,'first_sustained_visible_s':onset,'max_relative_inventory_residual':res,'legacy_detection_s':plan['legacy_detect_s'],'wall_s':finished['wall_s'] if finished else None}
(OUT/'metrics.json').write_text(json.dumps(metrics,indent=2))
# A self-contained diagnostic figure, with log scale explicitly labeled.
points=[]
for e in events:
 t=float(e['time_s']);s=max(float(e['sigma_top_node']),1e-18)
 points.append(f'{70+660*t/plan["horizon_s"]:.3f},{330-280*(math.log10(s)+18)/16:.3f}')
svg=['<svg xmlns="http://www.w3.org/2000/svg" width="800" height="400" viewBox="0 0 800 400">',
 '<rect width="800" height="400" fill="white"/>','<g font-family="Arial" font-size="12" fill="#222">',
 '<text x="70" y="22">Case 36877: no imposed perturbation (live diagnostic)</text>']
for exponent in [-18,-14,-10,-8,-6,-2]:
 y=330-280*(exponent+18)/16
 svg.append(f'<path d="M70 {y}H730" stroke="#ddd"/><text x="12" y="{y+4}">1e{exponent}</text>')
svg.extend(['<text x="70" y="375">Time (days); top-row concentration standard deviation (log scale)</text>',
 f'<text x="70" y="351">0</text><text x="690" y="351">{plan["horizon_s"]/86400:.3f}</text>',
 '<polyline fill="none" stroke="#1565a8" stroke-width="2" points="'+' '.join(points)+'"/>',
 '<path d="M70 155 H730" stroke="#c54a36" stroke-dasharray="5 4"/>',
 '<text x="460" y="148">Predeclared visibility threshold: 1e-8</text>','</g></svg>'])
(OUT/'sigma_history.svg').write_text('\n'.join(svg),encoding='utf-8')
status='已结束' if finished else '计算中'
verdict=(f'在 {onset:.3f} s 首次满足预设的持续横向波动判据。' if onset is not None else '截至当前输出，没有达到预设的持续横向波动判据。')
report=f'''# Case 36877 无外加扰动二维试算

状态：{status}。{verdict}

- 物理时刻：{metrics['last_time_s']:.6f} s / {plan['horizon_s']:.6f} s。
- 最新表面浓度节点标准差：{metrics['last_sigma']:.8e}；全过程最大值：{metrics['max_sigma']:.8e}。
- 最新顶部局部最大浓度：{metrics['last_top_max']:.8g}；面积均值：{metrics['last_top_mean']:.8g}。
- 盐/水累计残差相对各区域初始库存最大值：{res['salt']:.6e} / {res['water']:.6e}。
- 停止记录：`{stop['reason']}`；进程退出码：{metrics['returncode']}。

采用 case 36877 的原始 600×200 网格、物性、边界及 dt 上限 227.3717729539139 s，终点与旧运行一致（约 1.579 天）。将 Eigenprofile 改为 Constant，幅值和随机标准差设为 0；使用已验证诊断的 Paper 2 二进制，每 50 步保存场，新增局部阈值停止。输入、源码和 ELF 哈希见 `plan.json`。原扰动运行的 t_detect 约为 50394.334 s，仅作历史参照，不是同二进制的严格单因素配对。

预先登记的可见横向波动判据：顶部节点标准差超过 max(1e-8,1000×初始噪声底限)，连续 5 个接受步；噪声底限至少取 1e-16。该阈值用于判断数值噪声是否放大到可分辨水平，不等价于 LSA 临界时间、原 std-minimum 检测器或有效盐输运。

结果只能回答此网格、容差、时间窗口中的自发增长。未观测到增长不能证明数学稳定，也不能排除更长时间出现。若观测到增长，仍须核查保存的空间场与增长轨迹，才能判断指状结构及输运强度。

逐步波动图见 `sigma_history.svg`；纵轴为对数，零值仅为绘图放在 1e-18。

远端原始输出：`{launch['root']}`。求解器 PID 见 `evidence/solver_pid.json`。使用 `python scripts/collect_unperturbed_36877.py` 更新本报告及轻量诊断；该命令不启动或重跑求解器。
'''
if (OUT/'field_checks.json').exists():
    fields=json.loads((OUT/'field_checks.json').read_text())
    report+=f"\n已对全部 {fields['frames_checked']} 帧 VTU 核验结束 manifest 的 SHA-256，并计算每个高度的横向波动。所有保存场的最大横向浓度极差为 {fields['max_lateral_x_span']:.6e}，最大横向标准差为 {fields['max_lateral_x_std']:.6e}；没有可分辨的指状浓度结构。PVD 时间标签经过格式化取整，精确接受时刻以预算/事件 CSV 为准。\n"
(OUT/'RESULT_ZH.md').write_text(report,encoding='utf-8');print(json.dumps(metrics,indent=2))
