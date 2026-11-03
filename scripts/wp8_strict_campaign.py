"""Command 13: sequential gated T then independent continuous 4T, no restart."""
import base64,hashlib,json,subprocess,sys,tarfile
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'analysis/wp8_strict_background_20261103';CACHE=ROOT/'tmp/wp8_strict_background_20261103'
REMOTE='/temp/zhixin/paper2_wp8_strict_background_20261103';BUILD='/temp/zhixin/paper2_wp8_background_build_20261103_r03';BRIDGE='/mnt/d/3OutputFile/1Research/Article/22Surrogate_model/paper/2/scripts/frodo_bridge.py';T=136423.0637723484
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(p,v):p.parent.mkdir(parents=True,exist_ok=True);p.write_text(json.dumps(v,indent=2)+'\n')
def remote(code):
    r=subprocess.run(['wsl','python3',BRIDGE,'run','env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python3 -'],input=code.encode(),capture_output=True,timeout=60)
    if r.returncode:raise RuntimeError(r.stderr.decode(errors='replace'))
    return json.loads(r.stdout)
def prepare(stage):
    assert stage in ['T','4T'];d=OUT/stage;assert not d.exists()
    if stage=='4T':assert json.loads((OUT/'T/validation.json').read_text())['passed'], 'T gate required'
    prior=ROOT/'analysis/wp8_zero_mode_20261103';protocol=json.loads((prior/'background_protocol.json').read_text());assert protocol['native_short_validation']['passed']
    manifest=json.loads((ROOT/'model/transport_wp8_background_v1/source_manifest.json').read_text());build=json.loads((prior/'build/r03_attempt_01/result.json').read_text());assert build['source_sha256']=={k:v['sha256'] for k,v in manifest.items()}
    for k,v in manifest.items():assert sha(ROOT/'model/transport_wp8_background_v1/source'/k)==v['sha256']
    base=ROOT/'analysis/wp7_shortwave_validation_20261103/runs/WP7-1D-background-T/params.input';text=base.read_text().replace('Name = WP7-1D-background-T','Name = WP8-strict-'+stage).replace('TEnd = 136423.0637723484','TEnd = '+str(T*(1 if stage=='T' else 4)))
    d.mkdir(parents=True);(d/'params.input').write_text(text)
    plan=dict(stage=stage,reviewed_commit='14e15c7edb32a6a98faf7d915052590d35039ed6',authorization_sha256=sha(ROOT/'command/13.md'),previous_protocol_sha256=sha(prior/'background_protocol.json'),previous_failed_background_preserved=True,executable=build['executable']['path'],elf_sha256=build['executable']['sha256'],remote_source=BUILD+'/source',source_hashes=build['source_sha256'],input_sha256=sha(d/'params.input'),gates=protocol['gates'],wall_cap_s=900,archive_cap_mib=250,memory_gib=2,cpus=1,T_s=T,horizon_s=T*(1 if stage=='T' else 4),restart=False,previous_T_gate_sha256=sha(OUT/'T/validation.json') if stage=='4T' else None)
    save(d/'plan.json',plan);print('prepared',stage)
def launch(stage):
    d=OUT/stage;assert not (d/'launch.json').exists();p=json.loads((d/'plan.json').read_text());assert sha(d/'params.input')==p['input_sha256']
    if stage=='4T':assert json.loads((OUT/'T/validation.json').read_text())['passed'] and p['previous_T_gate_sha256']==sha(OUT/'T/validation.json')
    files={'params.input':base64.b64encode((d/'params.input').read_bytes()).decode(),'worker.py':base64.b64encode((ROOT/'scripts/wp8_strict_worker.py').read_bytes()).decode()}
    code=r'''import base64,json,os,subprocess,sys,time
from pathlib import Path
r=Path(REMOTE)/STAGE;assert not r.exists();r.mkdir(parents=True)
for n,b in FILES.items():(r/n).write_bytes(base64.b64decode(b))
def snap():return {x.split()[0]:list(map(int,x.split()[1:])) for x in Path('/proc/stat').read_text().splitlines() if x.startswith('cpu')}
a=snap();time.sleep(1);b=snap();usage={}
for k in b:
 q=[v-u for u,v in zip(a[k],b[k])];usage[k]=1-(q[3]+q[4])/sum(q) if sum(q) else 0
cpu=min(os.sched_getaffinity(0),key=lambda i:usage['cpu'+str(i)]);mem=int(next(x.split()[1] for x in Path('/proc/meminfo').read_text().splitlines() if x.startswith('MemAvailable:')))*1024;free=os.statvfs(r).f_bavail*os.statvfs(r).f_frsize
assert mem>2*1024**3 and free>1024**3
p=PLAN;p['cpu']=cpu;(r/'plan.json').write_text(json.dumps(p,indent=2)+'\n');res=dict(cpu=cpu,prior_cpu_usage=usage['cpu'+str(cpu)],memory_available_bytes=mem,free_disk_bytes=free,load=os.getloadavg());(r/'resources.json').write_text(json.dumps(res,indent=2))
with (r/'worker.log').open('w') as f:worker=subprocess.Popen([sys.executable,str(r/'worker.py'),str(r)],stdout=f,stderr=subprocess.STDOUT,start_new_session=True,stdin=subprocess.DEVNULL)
print(json.dumps(dict(pid=worker.pid,resources=res,plan=p)))
'''.replace('REMOTE',repr(REMOTE)).replace('STAGE',repr(stage)).replace('FILES',repr(files)).replace('PLAN',repr(p))
    r=remote(code);save(d/'plan.json',r.pop('plan'));save(d/'launch.json',r);print(json.dumps(r))
def status(stage):
    print(json.dumps(remote('import json\nfrom pathlib import Path\nr=Path('+repr(REMOTE+'/'+stage)+')\nprint(json.dumps({n:json.loads((r/n).read_text()) for n in ["monitor.json","finished.json","stop_record.json"] if (r/n).exists()}))'),indent=2))
def collect(stage):
    d=OUT/stage;cache=CACHE/stage;cache.mkdir(parents=True,exist_ok=True)
    r=remote('''import json,hashlib,base64
from pathlib import Path
r=Path('''+repr(REMOTE+'/'+stage)+''');assert (r/'finished.json').exists();files={}
for p in r.rglob('*'):
 if p.is_file() and p.suffix in ['.json','.csv','.log','.npz']:
  b=p.read_bytes();files[str(p.relative_to(r))]=dict(data=base64.b64encode(b).decode(),sha256=hashlib.sha256(b).hexdigest())
print(json.dumps(files))''')
    for n,v in r.items():
        b=base64.b64decode(v['data']);assert hashlib.sha256(b).hexdigest()==v['sha256'];dest=(cache if n.endswith('.npz') else d)/n;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(b)
    finished=json.loads((d/'finished.json').read_text());arc=cache/'fields.tar.gz'
    if not arc.exists() or sha(arc)!=finished['archive']['sha256']:
        local='/mnt/d/'+str(arc.relative_to(Path('D:/'))).replace('\\','/')
        subprocess.run(['wsl','python3',BRIDGE,'pull',finished['archive']['path'],local],check=True,timeout=120)
    assert sha(arc)==finished['archive']['sha256'];manifest=json.loads((d/'post/field_manifest.json').read_text());expected={v['name']:v['sha256'] for v in manifest};seen=set()
    with tarfile.open(arc,'r|gz') as tar:
        for m in tar:
            if m.isfile():assert m.name in expected and hashlib.file_digest(tar.extractfile(m),'sha256').hexdigest()==expected[m.name];seen.add(m.name)
    assert seen==set(expected);save(d/'archive_receipt.json',dict(verified_frames=len(seen),archive_sha256=sha(arc),profile_sha256=sha(cache/'post/profiles.npz'),local_archive=str(arc)));print('collected and verified',stage,len(seen))
if __name__=='__main__':globals()[sys.argv[1]](sys.argv[2])
