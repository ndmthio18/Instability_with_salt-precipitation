"""Run and retrieve a one-step conservative local-flux bridge smoke on frodo."""
import base64,csv,hashlib,json,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'analysis/wp7_shortwave_validation_20261103/local_flux_smoke'
SOURCE=ROOT/'analysis/wp7_shortwave_validation_20261103/runs/WP7-R48-M1-bench10'
BRIDGE='/mnt/d/3OutputFile/1Research/Article/22Surrogate_model/paper/2/scripts/frodo_bridge.py'
REMOTE='/temp/zhixin/paper2_wp7_localflux_smoke_20261103'
DT=113.68588647695695
def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def remote(code,timeout=120):
    q=subprocess.run(['wsl','python3',BRIDGE,'run','python3 -'],input=code.encode(),capture_output=True,timeout=timeout)
    if q.returncode:raise RuntimeError(q.stderr.decode(errors='replace'))
    return json.loads(q.stdout)
def prepare():
    assert not OUT.exists();OUT.mkdir(parents=True)
    text=(SOURCE/'params.input').read_text().replace('Name = WP7-R48-M1-bench10','Name = WP7-local-flux-smoke').replace('TEnd = 1136.8588647695694',f'TEnd = {DT}')
    marker='[Transport]\n';assert marker in text;text=text.replace(marker,marker+'WriteLocalFluxes = true\n',1)
    (OUT/'params.input').write_text(text);(OUT/'eigenprofile.dat').write_bytes((SOURCE/'eigenprofile.dat').read_bytes())
    plan=dict(id='WP7-local-flux-smoke',input_sha256=sha(OUT/'params.input'),profile_sha256=sha(OUT/'eigenprofile.dat'),
      executable_sha256='504cef7a8e95c9a3c6fe1a59c763cba43f5aceeeb169479d8781e0e7f4360387',TEnd_s=DT,accepted_steps=1,
      preregistration_sha256=sha(ROOT/'analysis/wp7_shortwave_validation_20261103/local_flux_smoke_E04.json'))
    (OUT/'plan.json').write_text(json.dumps(plan,indent=2)+'\n')
def execute():
    files={p.name:base64.b64encode(p.read_bytes()).decode() for p in [OUT/'params.input',OUT/'eigenprofile.dat',OUT/'plan.json']}
    code=r'''import base64,hashlib,json,os,subprocess,time
from pathlib import Path
r=Path(REMOTE);assert not r.exists();r.mkdir();files=FILES
for n,s in files.items():(r/n).write_bytes(base64.b64decode(s))
p=json.loads((r/'plan.json').read_text());exe=Path('/temp/zhixin/paper2_wp7_localflux_build_20261103/build-gcc-release/appl/simpleRichardsnc/paper2_transport')
assert hashlib.sha256(exe.read_bytes()).hexdigest()==p['executable_sha256']
def snap():return {x.split()[0]:list(map(int,x.split()[1:])) for x in Path('/proc/stat').read_text().splitlines() if x.startswith('cpu')}
a=snap();time.sleep(1);b=snap();usage={}
for k in b:
 d=[v-u for u,v in zip(a[k],b[k])];usage[k]=1-(d[3]+d[4])/sum(d) if sum(d) else 0
cpu=min(range(os.cpu_count()),key=lambda i:usage['cpu'+str(i)])
def bind():os.sched_setaffinity(0,{cpu})
env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1',MKL_NUM_THREADS='1')
with (r/'solver.log').open('w') as log:q=subprocess.run([str(exe),'params.input'],cwd=r,stdout=log,stderr=subprocess.STDOUT,env=env,preexec_fn=bind,timeout=120)
result=dict(returncode=q.returncode,cpu=cpu,affinity_usage=usage['cpu'+str(cpu)],elf_sha256=hashlib.sha256(exe.read_bytes()).hexdigest())
(r/'execution.json').write_text(json.dumps(result,indent=2));print(json.dumps(result))
'''.replace('REMOTE',repr(REMOTE)).replace('FILES',repr(files))
    result=remote(code,180);assert result['returncode']==0;(OUT/'execution.json').write_text(json.dumps(result,indent=2)+'\n')
def collect():
    code=r'''import base64,hashlib,json
from pathlib import Path
r=Path(REMOTE);files={}
for n in ['transport_budget.csv','local_vertical_flux.csv','transport_metadata.json','surface_events.csv','weighted_online.csv','stop_record.json','solver.log','execution.json']:
 p=r/n;b=p.read_bytes();files[n]=dict(data=base64.b64encode(b).decode(),sha256=hashlib.sha256(b).hexdigest())
print(json.dumps(files))
'''.replace('REMOTE',repr(REMOTE))
    files=remote(code)
    for n,rec in files.items():
        data=base64.b64decode(rec['data']);assert hashlib.sha256(data).hexdigest()==rec['sha256'];(OUT/n).write_bytes(data)
    budget={(r['step'],r['region']):r for r in csv.DictReader((OUT/'transport_budget.csv').open())}
    local={}
    for row in csv.DictReader((OUT/'local_vertical_flux.csv').open()):local.setdefault((row['step'],row['region']),[]).append(row)
    width=0.6/28;mapping={'water':('water_down_integrated_m2_s','cut_water_down_m_s'),'advective':('salt_adv_down_integrated_m2_s','cut_salt_adv_down_m_s'),'diffusive':('salt_diff_down_integrated_m2_s','cut_salt_diff_down_m_s'),'total':('salt_total_down_integrated_m2_s','cut_salt_total_down_m_s')}
    errors=[]
    for key,rows in local.items():
        for component,(column,target) in mapping.items():
            observed=sum(float(r[column]) for r in rows)/width;expected=float(budget[key][target]);absolute=abs(observed-expected);limit=1e-20+1e-12*abs(expected)
            errors.append(dict(step=int(key[0]),region=key[1],component=component,observed_m_s=observed,expected_m_s=expected,absolute_error_m_s=absolute,limit_m_s=limit,passed=absolute<=limit))
    stop=json.loads((OUT/'stop_record.json').read_text());result=dict(passed=all(r['passed'] for r in errors) and stop['reason']=='requested_horizon' and stop['accepted_steps']==1,
      comparisons=len(errors),max_absolute_error_m_s=max(r['absolute_error_m_s'] for r in errors),max_error_to_limit=max(r['absolute_error_m_s']/r['limit_m_s'] for r in errors),rows=sum(len(v) for v in local.values()),stop=stop)
    (OUT/'bridge_errors.csv').write_text('step,region,component,observed_m_s,expected_m_s,absolute_error_m_s,limit_m_s,passed\n'+'\n'.join(','.join(str(r[k]) for k in ['step','region','component','observed_m_s','expected_m_s','absolute_error_m_s','limit_m_s','passed']) for r in errors)+'\n')
    (OUT/'bridge_result.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2));assert result['passed']
if __name__=='__main__':prepare();execute();collect()
