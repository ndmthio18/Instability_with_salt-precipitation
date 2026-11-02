"""One extended high-permeability sandy case: prepare/launch, collect, or watch."""
import argparse,base64,csv,hashlib,json,subprocess,sys,time
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
from onset_semantics import FrozenOnsets, standard_fields
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'analysis/highk_unperturbed_42929_20261102'
REMOTE='/temp/zhixin/paper2_transport_20261102/highk_unperturbed_42929_20261102_r01'
BRIDGE='/mnt/d/3OutputFile/1Research/Article/22Surrogate_model/paper/2/scripts/frodo_bridge.py'
def remote(code):
 p=subprocess.run(['wsl','python3',BRIDGE,'run','python3 -'],input=code.encode(),capture_output=True,timeout=90)
 if p.returncode:raise RuntimeError(p.stderr.decode(errors='replace'))
 return json.loads(p.stdout)
def launch():
 OUT.mkdir(exist_ok=False)
 audit=ROOT/'analysis/phase0_20261102'
 inventory={r['run_key']:r for r in json.loads((audit/'remote_frodo_inventory.json').read_text())['runs']}
 candidates=[]
 onsets=FrozenOnsets()
 for row in csv.DictReader((audit/'run_audit.csv').open(encoding='utf-8-sig')):
  row.update(standard_fields(onsets.run_by_key(row['run_key'])))
  if row['soil_group']=='sandy' and row['onset_class']=='resolved_onset' and row['run_key'] in inventory:
   row['K']=float(inventory[row['run_key']]['parameters']['Problem.Permeability']);candidates.append(row)
 candidates.sort(key=lambda r:-r['K']);selected=candidates[0];assert selected['case_id']=='42929'
 with (OUT/'screening.csv').open('w',newline='') as f:
  names=['case_id','run_key','K','onset_s','dumux_t_detect_s','lsa_t_crit_s','lsa_t_seed_theta010_s','lsa_t_sat_reference_s','amplification_since_min','effective_saturation_top'];w=csv.DictWriter(f,fieldnames=names,extrasaction='ignore');w.writeheader();w.writerows(candidates)
 original=inventory[selected['run_key']];v=original['parameters'].copy();end=10*selected['dumux_t_detect_s']
 changes={'Problem.Name':'unperturbed_42929','Problem.PerturbationType':'Constant','Problem.Amplitude':'0','Problem.StandardDeviationSaltMoleFraction':'0','TimeLoop.TEnd':format(end,'.17g'),'TimeLoop.OutputInterval':'100','Transport.StopAtSolubility':'true','Transport.SolubilityMoleFraction':'0.1'}
 v.update(changes);v.pop('Problem.EigenprofilePath',None);sections={}
 for k,val in v.items():
  section,key=k.split('.',1);sections.setdefault(section,{})[key]=val
 text='\n\n'.join('['+s+']\n'+'\n'.join(k+' = '+v for k,v in entries.items()) for s,entries in sections.items())+'\n'
 build=json.loads((ROOT/'analysis/transport_build_20261102/attempt_02/result.json').read_text())
 for n,h in build['source_sha256'].items():assert hashlib.sha256((ROOT/'model/transport_v1/source'/n).read_bytes()).hexdigest()==h
 plan={'case_id':42929,'soil_label':'high-permeability sandy; coarse grain size not verified','selection':selected,'original_parameters':original['parameters'],'parameters':v,'changes':changes,'horizon_multiple_of_legacy_detect':10,'horizon_s':end,'old_detect_s':float(selected['onset_s']),'executable':build['executable'],'source_sha256':build['source_sha256'],'input_sha256':hashlib.sha256(text.encode()).hexdigest(),'remote':REMOTE,'criterion':'top sigma > 1e-8 for 5 accepted steps; this is visibility, not LSA onset','wall_timeout_s':86400,'old_small_file_hashes':original['small_file_hashes'],'old_path':original['path']}
 plan['legacy_onset_reference']=standard_fields(selected)
 (OUT/'plan.json').write_text(json.dumps(plan,indent=2));(OUT/'params.input').write_text(text)
 worker='''import subprocess,json,hashlib,os,time
from pathlib import Path
r=Path(__file__).resolve().parent;p=json.loads((r/'plan.json').read_text());exe=Path(p['executable']['path'])
assert hashlib.sha256(exe.read_bytes()).hexdigest()==p['executable']['sha256']
assert hashlib.sha256((r/'params.input').read_bytes()).hexdigest()==p['input_sha256']
start=time.monotonic()
with (r/'solver.log').open('w') as f:
 proc=subprocess.Popen(['nice','-n','10',str(exe),'params.input'],cwd=r,stdout=f,stderr=subprocess.STDOUT,env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1'))
 (r/'solver_pid.json').write_text(json.dumps({'pid':proc.pid}))
 try:rc=proc.wait(timeout=p['wall_timeout_s'])
 except subprocess.TimeoutExpired:proc.terminate();proc.wait();rc='wall_timeout'
 result={'returncode':rc,'wall_s':time.monotonic()-start,'stop':json.loads((r/'stop_record.json').read_text()) if (r/'stop_record.json').exists() else None}
 result['files']={f.name:{'sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'bytes':f.stat().st_size} for f in r.iterdir() if f.is_file()}
 (r/'finished.json').write_text(json.dumps(result,indent=2))
'''
 (OUT/'worker.py').write_text(worker)
 code='''import subprocess,json,hashlib,sys,xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
p=PLAN;old=Path(p['old_path']);checks={}
for n,h in p['old_small_file_hashes'].items():
 q=old/n
 assert hashlib.sha256(q.read_bytes()).hexdigest()==h,n
 checks[n]=h
last=sorted(old.glob('*.vtu'))[-1];arrays={a.attrib.get('Name'):np.fromstring(a.text or '',sep=' ') for a in ET.parse(last).getroot().iter('DataArray') if a.attrib.get('format')=='ascii' and a.attrib.get('Name') in ['rhoW','moleFraction']}
v=p['original_parameters'];x=arrays['moleFraction'];mw=float(v['Problem.MolarMassH2O']);ms=float(v['Problem.MolarMassNaCl']);w=x*ms/((1-x)*mw+x*ms)
rho=float(v['Problem.PureWaterDensity'])+1000*w*(float(v['Problem.DensityLinearCoefficient'])+float(v['Problem.DensityQuadraticCoefficient'])*w)
error=float(np.max(np.abs(rho-arrays['rhoW'])));assert error<1e-8,error
r=Path(p['remote']);r.mkdir(exist_ok=False)
(r/'plan.json').write_text(json.dumps(p,indent=2));(r/'params.input').write_text(INPUT);(r/'worker.py').write_text(WORKER)
provenance={'old_small_file_hashes_verified':checks,'old_final_field_sha256':hashlib.sha256(last.read_bytes()).hexdigest(),'old_field_rho_max_error':error}
(r/'legacy_check.json').write_text(json.dumps(provenance,indent=2))
with (r/'worker.log').open('w') as f:
 proc=subprocess.Popen([sys.executable,str(r/'worker.py')],stdin=subprocess.DEVNULL,stdout=f,stderr=subprocess.STDOUT,start_new_session=True)
print(json.dumps({'worker_pid':proc.pid,'root':str(r),'legacy_check':provenance}))
'''.replace('PLAN',repr(plan)).replace('INPUT',repr(text)).replace('WORKER',repr(worker))
 result=remote(code);(OUT/'launch.json').write_text(json.dumps(result,indent=2));print(json.dumps(result,indent=2))
def collect():
 code='''import json,base64
from pathlib import Path
r=Path(REMOTE);files={}
for name in ['stop_record.json','surface_events.csv','transport_budget.csv','transport_metadata.json','finished.json','legacy_check.json','solver_pid.json','worker.log']:
 p=r/name
 if p.exists():
  b=p.read_bytes()
  if name.endswith('.csv') and b and not b.endswith(b'\\n'):b=b[:b.rfind(b'\\n')+1]
  files[name]=base64.b64encode(b).decode()
log=r/'solver.log'
with log.open('rb') as f:f.seek(max(0,log.stat().st_size-6000));tail=f.read().decode(errors='replace')
print(json.dumps({'files':files,'tail':tail}))
'''.replace('REMOTE',repr(REMOTE))
 data=remote(code);e=OUT/'evidence';e.mkdir(exist_ok=True)
 for n,b in data['files'].items():(e/n).write_bytes(base64.b64decode(b))
 (e/'solver_tail.log').write_text(data['tail'])
 if not (e/'surface_events.csv').exists():return {'completed':False,'status':'initializing'}
 rows=list(csv.DictReader((e/'surface_events.csv').open()));plan=json.loads((OUT/'plan.json').read_text());last=rows[-1]
 streak=0;hit=None
 for i,r in enumerate(rows):
  streak=streak+1 if float(r['sigma_top_node'])>1e-8 else 0
  if streak==5 and hit is None:hit=float(rows[i-4]['time_s'])
 done=json.loads((e/'finished.json').read_text()) if (e/'finished.json').exists() else None
 budget=list(csv.DictReader((e/'transport_budget.csv').open()));initial={r['region']:r for r in budget if r['step']=='0'}
 residual={s:max(abs(float(r[s+'_residual_m']))/float(initial[r['region']]['inventory_'+s+'_m']) for r in budget) for s in ['salt','water']}
 result={'completed':done is not None,'returncode':done['returncode'] if done else None,'time_s':float(last['time_s']),'horizon_s':plan['horizon_s'],'sigma_max':max(float(r['sigma_top_node']) for r in rows),'sigma_last':float(last['sigma_top_node']),'visible_growth_first_s':hit,'x_top_max':float(last['x_top_max']),'relative_budget_residual_max':residual,'stop':json.loads((e/'stop_record.json').read_text())}
 (OUT/'metrics.json').write_text(json.dumps(result,indent=2))
 report=f'''# 高渗透率砂质案例 42929：延长无扰动试算

状态：{'已结束' if done else '计算中'}。最新物理时间 {result['time_s']:.3f} s / {plan['horizon_s']:.3f} s。

按已访问、已有 resolved onset 的砂质运行筛选，case 42929 的 K 最高，为 2.47926725055e-11 m²（旧 case 36877 为 1.26138608e-11）。数据库标签为 sandy，没有粒径证据，因此本报告称高渗透率砂质候选，不将其认定为实测粗砂。筛选表见 `screening.csv`。

旧扰动运行 t_detect={plan['old_detect_s']:.6f} s，约 {plan['old_detect_s']/3600:.4f} h，属于 std-minimum 后增长的检测，并非已形成成熟指状结构。新无扰动运行采用其 1200×200 网格、物性、边界和 dt 上限 29.36067352995103 s，延长到 10 t_detect={plan['horizon_s']:.6f} s（{plan['horizon_s']/3600:.3f} h）。LSA 阈值参考约 11.042 day，不用它代替实际二维阈值；新运行若 max(x) 达到 0.1 则提前停止。

旧输入/profile/manifest 哈希和末帧 rho(x) 一致性已核对，见 `launch.json`。新运行使用已验证诊断的 Paper 2 ELF；历史二进制身份不完整，不能称两次运行只有扰动一个差异。此次先利用已有检测时刻确定窗口，未另行启动有扰动腿或完整 Pilot A/B。

当前表面节点 sigma={result['sigma_last']:.8e}，全过程最大={result['sigma_max']:.8e}；顶部 max(x)={result['x_top_max']:.8g}。预登记可见波动判据为 sigma>1e-8 连续 5 个接受步，首次满足时刻={hit}。该判据不等价于 LSA onset 或有效输运。

盐/水最大累计库存相对残差：{residual['salt']:.5e}/{residual['water']:.5e}。停止状态 `{result['stop']['reason']}`，退出码 {result['returncode']}。

计算终点与判据在运行前登记。无增长仅适用于这个窗口；若出现增长，应进一步核查空间场和预算，不能只看 std 判据认定成熟对流。

远端：`{REMOTE}`。运行入口 `python scripts/run_highk_unperturbed.py launch` 拒绝已有目录；`collect` 收集结果；`watch` 每 60 s 收集一次，最多 25 h。
'''
 (OUT/'RESULT_ZH.md').write_text(report,encoding='utf-8');return result
def watch():
 for _ in range(1500):
  try:
   r=collect();(OUT/'watch_status.json').write_text(json.dumps(r,indent=2))
   if r['completed']:
    p=ROOT/'PROGRESS.md';t=p.read_text(encoding='utf-8');now=datetime.now(ZoneInfo('Europe/Berlin')).strftime('%Y-%m-%d %H:%M')
    t+=f"\n- {now} Codex: 高渗透率砂质 case 42929 延长无扰动试算结束，退出码 {r['returncode']}，物理时刻 {r['time_s']:.3f} s，最大 sigma={r['sigma_max']:.3e}；结果见 `analysis/highk_unperturbed_42929_20261102/RESULT_ZH.md`。\n";p.write_text(t,encoding='utf-8');return
  except Exception as ex:(OUT/'watch_last_error.txt').write_text(repr(ex))
  time.sleep(60)
if __name__=='__main__':
 ap=argparse.ArgumentParser();ap.add_argument('mode',choices=['launch','collect','watch']);a=ap.parse_args()
 if a.mode=='launch':launch()
 elif a.mode=='collect':print(json.dumps(collect(),indent=2))
 else:watch()
