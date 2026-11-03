"""Command 15: ONE continuous corrected strict background, t=0 to 4T. No restart/replay/retry."""
import base64,json,subprocess,sys,tarfile,hashlib
from pathlib import Path
import wp8_strict_campaign as c

ROOT=c.ROOT
OUT=ROOT/'analysis/wp8_full4T_20261103'
CACHE=ROOT/'tmp/wp8_full4T_20261103'
REMOTE='/temp/zhixin/paper2_wp8_full4T_20261103'
STAGE='run'
PROTOCOL=ROOT/'analysis/wp8_newton_correction_20261103/full4T_protocol.json'
PENDING=ROOT/'analysis/wp8_newton_correction_20261103/pending_full4T.input'
AUTH=ROOT/'command/15.md'
MODEL='transport_wp8_background_v2'
T=c.T
sha=c.sha
save=c.save

def remote(code,timeout=900):
    r=subprocess.run(['wsl','python3',c.BRIDGE,'run','env OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 python3 -'],
                     input=code.encode(),capture_output=True,timeout=timeout)
    if r.returncode:raise RuntimeError(r.stderr.decode(errors='replace'))
    return json.loads(r.stdout)

def _protocol():
    p=json.loads(PROTOCOL.read_text())
    assert p['decision']=='GO_FOR_USER_AUTHORIZED_FULL_4T_VALIDATION_ONLY' and not p['launched']
    assert p['restart'] is False and p['replay'] is False and p['original_initial_state'] is True
    return p

def prepare():
    d=OUT/STAGE;assert not d.exists(),'one authorized run only'
    p=_protocol()
    assert sha(PENDING)==p['input_sha256']
    man=json.loads((ROOT/'model'/MODEL/'source_manifest.json').read_text())
    assert {k:v['sha256'] for k,v in man.items()}==p['source_hashes']
    for k,v in man.items():assert sha(ROOT/'model'/MODEL/'source'/k)==v['sha256'],k
    hist=ROOT/'analysis/wp8_strict_background_20261103'
    histT=json.loads((hist/'T/validation.json').read_text());hist4=json.loads((hist/'4T/validation.json').read_text())
    assert histT['passed'] and not hist4['passed'] and hist4['accepted_prefix_passed']
    d.mkdir(parents=True);(d/'params.input').write_bytes(PENDING.read_bytes())
    assert sha(d/'params.input')==p['input_sha256']
    plan=dict(stage=STAGE,command='command/15.md',authorization_sha256=sha(AUTH),
        reviewed_commit=p['reviewed_commit'],protocol_sha256=sha(PROTOCOL),
        executable=p['executable']['path'],elf_sha256=p['executable']['sha256'],
        remote_source='/temp/zhixin/paper2_wp8_background_v2_build_20261103/source',
        source_hashes=p['source_hashes'],source_manifest_sha256=p['source_manifest_sha256'],
        input_sha256=p['input_sha256'],gates=p['gates'],criterion=p['criterion'],
        dt_s=p['dt_s'],horizon_s=p['horizon_s'],T_s=T,steps_cap=p['max_steps'],
        checkpoints_steps=p['checkpoints_steps'],threshold=p['threshold'],
        wall_cap_s=p['cost']['wall_hard_cap_s'],archive_cap_mib=p['cost']['archive_cap_mib_per_end'],
        memory_gib=p['cost']['solver_address_space_gib'],cpus=p['cost']['cpus'],
        minimum_free_disk_gib=p['cost']['minimum_free_disk_gib'],
        restart=False,replay=False,diagnostic_replay=False,original_initial_state=True,
        historical_prefix_steps=4413,historical_T_gate_sha256=sha(hist/'T/validation.json'),
        historical_4T_prefix_sha256=sha(hist/'4T/validation.json'),
        no_automatic_retry=True,full_4T_authorized=True)
    save(d/'plan.json',plan);print(json.dumps(plan,indent=2))

LAUNCH_CODE = r'''import base64,hashlib,json,os,subprocess,sys,time
from pathlib import Path
r=Path(REMOTE)/STAGE;assert not r.exists();r.mkdir(parents=True)
for n,b in FILES.items():(r/n).write_bytes(base64.b64decode(b))
p=PLAN
assert hashlib.sha256((r/'params.input').read_bytes()).hexdigest()==p['input_sha256']
exe=Path(p['executable']);assert hashlib.sha256(exe.read_bytes()).hexdigest()==p['elf_sha256']
src=Path(p['remote_source'])
for n,h in p['source_hashes'].items():assert hashlib.sha256((src/n).read_bytes()).hexdigest()==h,n
def snap():return {x.split()[0]:list(map(int,x.split()[1:])) for x in Path('/proc/stat').read_text().splitlines() if x.startswith('cpu')}
a=snap();time.sleep(1);b=snap();usage={}
for k in b:
 q=[v-u for u,v in zip(a[k],b[k])];usage[k]=1-(q[3]+q[4])/sum(q) if sum(q) else 0
cpu=min(os.sched_getaffinity(0),key=lambda i:usage['cpu'+str(i)])
mem=int(next(x.split()[1] for x in Path('/proc/meminfo').read_text().splitlines() if x.startswith('MemAvailable:')))*1024
free=os.statvfs(r).f_bavail*os.statvfs(r).f_frsize
assert mem>2*1024**3 and free>p['minimum_free_disk_gib']*1024**3
p['cpu']=cpu;(r/'plan.json').write_text(json.dumps(p,indent=2)+'\n')
res=dict(cpu=cpu,prior_cpu_usage=usage['cpu'+str(cpu)],memory_available_bytes=mem,free_disk_bytes=free,load=os.getloadavg(),live_check=True)
(r/'resources.json').write_text(json.dumps(res,indent=2))
with (r/'worker.log').open('w') as f:
 w=subprocess.Popen([sys.executable,str(r/'worker.py'),str(r)],stdout=f,stderr=subprocess.STDOUT,start_new_session=True,stdin=subprocess.DEVNULL)
print(json.dumps(dict(pid=w.pid,resources=res,plan=p)))
'''

def launch():
    d=OUT/STAGE;assert not (d/'launch.json').exists(),'already launched; no retry under command 15'
    p=json.loads((d/'plan.json').read_text());assert sha(d/'params.input')==p['input_sha256']
    assert sha(AUTH)==p['authorization_sha256']
    files={'params.input':base64.b64encode((d/'params.input').read_bytes()).decode(),
           'worker.py':base64.b64encode((ROOT/'scripts/wp8_strict_worker.py').read_bytes()).decode()}
    code=LAUNCH_CODE.replace('REMOTE',repr(REMOTE)).replace('STAGE',repr(STAGE)).replace('FILES',repr(files)).replace('PLAN',repr(p))
    r=remote(code,timeout=180);save(d/'plan.json',r.pop('plan'));save(d/'launch.json',r);print(json.dumps(r,indent=2))

STATUS_CODE=('import json\nfrom pathlib import Path\nr=Path(%r)\n'
             'print(json.dumps({n:json.loads((r/n).read_text()) for n in '
             '["monitor.json","finished.json","stop_record.json","resources.json"] if (r/n).exists()}))')

COLLECT_CODE='''import json,hashlib,base64
from pathlib import Path
r=Path(%r);assert (r/'finished.json').exists();files={}
for p in r.rglob('*'):
 if p.is_file() and p.suffix in ['.json','.csv','.log','.npz']:
  b=p.read_bytes();files[str(p.relative_to(r))]=dict(data=base64.b64encode(b).decode(),sha256=hashlib.sha256(b).hexdigest())
print(json.dumps(files))'''

def status():
    print(json.dumps(remote(STATUS_CODE%(REMOTE+'/'+STAGE),timeout=180),indent=2))

def collect():
    d=OUT/STAGE;cache=CACHE/STAGE;cache.mkdir(parents=True,exist_ok=True)
    r=remote(COLLECT_CODE%(REMOTE+'/'+STAGE),timeout=1800)
    for n,v in r.items():
        b=base64.b64decode(v['data']);assert hashlib.sha256(b).hexdigest()==v['sha256']
        dest=(cache if n.endswith('.npz') else d)/n;dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(b)
    finished=json.loads((d/'finished.json').read_text());arc=cache/'fields.tar.gz'
    if not arc.exists() or sha(arc)!=finished['archive']['sha256']:
        local='/mnt/d/'+str(arc.relative_to(Path('D:/'))).replace('\\','/')
        subprocess.run(['wsl','python3',c.BRIDGE,'pull',finished['archive']['path'],local],check=True,timeout=1800)
    assert sha(arc)==finished['archive']['sha256']
    manifest=json.loads((d/'post/field_manifest.json').read_text());expected={v['name']:v['sha256'] for v in manifest};seen=set()
    with tarfile.open(arc,'r|gz') as tar:
        for m in tar:
            if m.isfile():
                assert m.name in expected and hashlib.file_digest(tar.extractfile(m),'sha256').hexdigest()==expected[m.name];seen.add(m.name)
    assert seen==set(expected)
    save(d/'archive_receipt.json',dict(verified_frames=len(seen),archive_sha256=sha(arc),
        profile_sha256=sha(cache/'post/profiles.npz'),local_archive=str(arc),
        remote_archive=finished['archive']['path'],archive_bytes=finished['archive']['bytes']))
    print('collected and verified',len(seen))

if __name__=='__main__':globals()[sys.argv[1]]()
