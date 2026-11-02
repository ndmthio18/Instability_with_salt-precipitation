"""Bounded local watcher for the one authorized remote calculation."""
import json,subprocess,sys,time
from pathlib import Path
from datetime import datetime
from zoneinfo import ZoneInfo
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'analysis/unperturbed_36877_20261102'
end=time.monotonic()+13*3600
while time.monotonic()<end:
    p=subprocess.run([sys.executable,str(ROOT/'scripts/collect_unperturbed_36877.py')],capture_output=True,text=True,timeout=90)
    now=datetime.now(ZoneInfo('Europe/Berlin')).strftime('%Y-%m-%d %H:%M')
    status={'checked_berlin':now,'collector_returncode':p.returncode,'stderr':p.stderr[-2000:]}
    if p.returncode==0:
        m=json.loads((OUT/'metrics.json').read_text());status.update(last_time_s=m['last_time_s'],completed=m['completed'])
        if m['completed']:
            progress=ROOT/'PROGRESS.md';text=progress.read_text(encoding='utf-8')
            item='在 frodo 完成 case 36877 无外加扰动二维试算，检查数值噪声是否在原运行窗口内增长。'
            if m['returncode']==0:text=text.replace('- [ ] '+item,'- [x] '+item)
            result='达到持续横向波动判据' if m['first_sustained_visible_s'] is not None else '未达到持续横向波动判据'
            text+=f"\n- {now} Codex: case 36877 无外加扰动试算结束，退出码 {m['returncode']}，达到 {m['last_time_s']:.3f} s，{result}，最大 sigma={m['max_sigma']:.3e}；结果见 `analysis/unperturbed_36877_20261102/RESULT_ZH.md`。\n"
            progress.write_text(text,encoding='utf-8')
            (OUT/'watch_status.json').write_text(json.dumps(status,indent=2));break
    (OUT/'watch_status.json').write_text(json.dumps(status,indent=2))
    time.sleep(60)
else:
    (OUT/'watch_timeout.json').write_text(json.dumps({'reason':'13-hour local monitoring bound reached; inspect remote worker'}))
