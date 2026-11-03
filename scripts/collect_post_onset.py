"""WP4 read-only collection. Never invokes or alters a solver/run."""
import base64,hashlib,json,subprocess,sys
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'analysis/post_onset_amplification_20261103'; OUT.mkdir(exist_ok=True)
CACHE=ROOT/'tmp/post_onset_amplification_20261103/series'; CACHE.mkdir(parents=True,exist_ok=True)
def sha(b):return hashlib.sha256(b).hexdigest()
def remote(code,timeout=180):
 p=subprocess.run(['wsl','python3','/mnt/d/3OutputFile/1Research/Article/22Surrogate_model/paper/2/scripts/frodo_bridge.py','run','env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 python3 -'],input=code.encode(),capture_output=True,timeout=timeout)
 if p.returncode:raise RuntimeError(p.stderr.decode(errors='replace'))
 return json.loads(p.stdout)
def main():
 inv=json.loads((ROOT/'analysis/phase0_20261102/remote_frodo_inventory.json').read_text())['runs']
 targets=[];manifest=[]
 for r in inv:
  names=[n for n in r['files'] if n.endswith('_sigmatop.csv')]
  if not names:continue
  key=r['run_key']; dest=CACHE/(sha(key.encode())[:20]+'.csv')
  if not dest.exists() or sha(dest.read_bytes())!=r['sigmatop_sha256']:targets.append(dict(key=key,path=r['path']+'/'+names[0],expected=r['sigmatop_sha256']))
  manifest.append(dict(run_key=key,path=r['path']+'/'+names[0],sha256=r['sigmatop_sha256'],local=str(dest.relative_to(ROOT)).replace('\\','/')))
 # One bounded read per named immutable CSV, no remote writes.
 code='import os,json,base64,hashlib,pathlib\nresult=[]\n'
 code+='targets='+repr(targets)+'\n'
 code+='for t in targets:\n b=pathlib.Path(t["path"]).read_bytes(); h=hashlib.sha256(b).hexdigest(); assert h==t["expected"],t["key"]; result.append(dict(key=t["key"],data=base64.b64encode(b).decode()))\n'
 code+='print(json.dumps(dict(load=os.getloadavg(),cpu_count=os.cpu_count(),series=result)))\n'
 if targets:
  result=remote(code)
  for r in result['series']:(CACHE/(sha(r['key'].encode())[:20]+'.csv')).write_bytes(base64.b64decode(r['data']))
  (OUT/'collection_resource_snapshot.json').write_text(json.dumps({k:v for k,v in result.items() if k!='series'},indent=2))
 for r in manifest:assert sha((ROOT/r['local']).read_bytes())==r['sha256']
 (OUT/'series_source_manifest.json').write_text(json.dumps(manifest,indent=2))
 p=subprocess.run(['ssh','-o','BatchMode=yes','-o','ConnectTimeout=8','cxhpc','hostname'],capture_output=True,timeout=20)
 (OUT/'cxhpc_access.json').write_text(json.dumps(dict(returncode=p.returncode,stderr=p.stderr.decode(errors='replace'),stdout=p.stdout.decode(errors='replace'),bypass=False),indent=2))
 print(json.dumps(dict(verified_series=len(manifest),newly_read=len(targets),cxhpc_returncode=p.returncode)))
if __name__=='__main__':main()
