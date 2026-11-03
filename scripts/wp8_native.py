"""Only isolated build collection and a bounded ten-step native constraint bridge.

No 4T/production entry point is provided. Existing remote directories refuse rerun.
"""
import base64,hashlib,json,subprocess,sys,csv
from pathlib import Path
import numpy as np
from analyze_wp8 import ROOT,OUT,RUNS,pairs,rows
from single_wavelength_post import widths,params
BRIDGE='/mnt/d/3OutputFile/1Research/Article/22Surrogate_model/paper/2/scripts/frodo_bridge.py'
REV=sys.argv[2] if len(sys.argv)>2 else 'r03'
assert REV in ['r01','r02','r03']
SUFFIX='' if REV=='r01' else '_'+REV
BUILD='/temp/zhixin/paper2_wp8_background_build_20261103'+SUFFIX;REMOTE='/temp/zhixin/paper2_wp8_background_smoke10_20261103'+SUFFIX;DEST=OUT/('background_smoke10'+SUFFIX);CACHE=ROOT/('tmp/wp8_zero_mode_20261103/background_smoke10'+SUFFIX)
BUILD_RESULT=OUT/'build'/('attempt_01' if REV=='r01' else REV+'_attempt_01')/'result.json'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def remote(code,timeout=150):
    r=subprocess.run(['wsl','python3',BRIDGE,'run','python3 -'],input=code.encode(),capture_output=True,timeout=timeout)
    if r.returncode:raise RuntimeError(r.stderr.decode(errors='replace'))
    return json.loads(r.stdout)
def status():
    r=remote('import json\nfrom pathlib import Path\np=Path('+repr(BUILD+'/attempt_01/result.json')+')\nprint(p.read_text() if p.exists() else json.dumps({"status":"building"}))');print(json.dumps({k:r.get(k) for k in ['status','configure_exit','build_exit','executable']},ensure_ascii=True))
def collectbuild():
    s=(ROOT/'scripts/collect_wp7_build.py').read_text().replace('analysis/wp7_shortwave_validation_20261103/build','analysis/wp8_zero_mode_20261103/build').replace('ROOT=Path(__file__).resolve().parents[1]','ROOT=Path('+repr(str(ROOT))+')')
    launch='attempt_01_launch.json' if REV=='r01' else REV+'_attempt_01_launch.json'
    s=s.replace("sorted(OUT.glob('attempt_*_launch.json'))[-1]","(OUT/"+repr(launch)+")").replace("dest=OUT/Path(launch['attempt']).name",'dest=Path('+repr(str(BUILD_RESULT.parent))+')')
    exec(compile(s,'collect_wp8_build','exec'),{'__file__':__file__})
def smoke():
    build=json.loads(BUILD_RESULT.read_text());assert build['build_exit']==0
    assert not DEST.exists();DEST.mkdir();CACHE.mkdir(parents=True)
    text=(ROOT/RUNS['B']/'params.input').read_text().replace('Name = WP7-1D-background-T','Name = WP8-constrained-smoke10').replace('TEnd = 136423.0637723484','TEnd = 1136.8588647695694')
    (DEST/'params.input').write_text(text)
    protocol=json.loads((OUT/'background_protocol.json').read_text())
    plan=dict(executable_sha256=build['executable']['sha256'],input_sha256=sha(DEST/'params.input'),source_manifest_sha256=sha(ROOT/'model/transport_wp8_background_v1/source_manifest.json'),short_steps=10,wall_limit_s=120,gates=protocol['gates'],long_simulation_authorized=False)
    (DEST/'plan.json').write_text(json.dumps(plan,indent=2)+'\n')
    files={p.name:base64.b64encode(p.read_bytes()).decode() for p in DEST.iterdir() if p.is_file()}
    code=r'''import base64,json,hashlib,subprocess,time,os,tarfile
from pathlib import Path
r=Path(REMOTE);assert not r.exists();r.mkdir()
for n,b in FILES.items():(r/n).write_bytes(base64.b64decode(b))
p=json.loads((r/'plan.json').read_text());exe=Path(BUILD)/'build-gcc-release/appl/simpleRichardsnc/paper2_transport';assert hashlib.sha256(exe.read_bytes()).hexdigest()==p['executable_sha256']
def snap():return {x.split()[0]:list(map(int,x.split()[1:])) for x in Path('/proc/stat').read_text().splitlines() if x.startswith('cpu')}
a=snap();time.sleep(1);b=snap();usage={}
for k in b:
 d=[v-u for u,v in zip(a[k],b[k])];usage[k]=1-(d[3]+d[4])/sum(d) if sum(d) else 0
cpu=min(os.sched_getaffinity(0),key=lambda i:usage['cpu'+str(i)])
mem=int(next(x.split()[1] for x in Path('/proc/meminfo').read_text().splitlines() if x.startswith('MemAvailable:')))*1024
assert mem>2*1024**3
def bind():os.sched_setaffinity(0,{cpu})
t=time.time();result=dict(cpu=cpu,prior_cpu_usage=usage['cpu'+str(cpu)],memory_available_bytes=mem,elf_sha256=p['executable_sha256'])
try:
 with (r/'solver.log').open('w') as f:q=subprocess.run([str(exe),'params.input'],cwd=r,stdout=f,stderr=subprocess.STDOUT,env=dict(os.environ,OMP_NUM_THREADS='1',OPENBLAS_NUM_THREADS='1'),preexec_fn=bind,timeout=120)
 result['returncode']=q.returncode
except subprocess.TimeoutExpired:result['timeout']=True
result['wall_s']=time.time()-t
arc=r/'fields.tar.gz'
with tarfile.open(arc,'w:gz') as tar:
 for f in sorted(r.glob('*.vtu')):tar.add(f,arcname=f.name)
result['archive']=dict(path=str(arc),sha256=hashlib.sha256(arc.read_bytes()).hexdigest(),bytes=arc.stat().st_size)
(r/'execution.json').write_text(json.dumps(result,indent=2));print(json.dumps(result))
'''.replace('REMOTE',repr(REMOTE)).replace('FILES',repr(files)).replace('BUILD',repr(BUILD))
    result=remote(code);(DEST/'execution.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result))
def collect():
    data=remote('''import base64,json,hashlib
from pathlib import Path
r=Path('''+repr(REMOTE)+''');out={}
for p in r.iterdir():
 if p.is_file() and p.suffix in ['.json','.csv','.log','.vtu','.gz']:
  b=p.read_bytes();out[p.name]=dict(sha256=hashlib.sha256(b).hexdigest(),bytes=len(b),data=base64.b64encode(b).decode())
print(json.dumps(out))''')
    receipt={}
    for n,r in data.items():
        b=base64.b64decode(r.pop('data'));assert hashlib.sha256(b).hexdigest()==r['sha256'];dest=CACHE if n.endswith(('.vtu','.gz')) else DEST;(dest/n).write_bytes(b);receipt[n]=r
    (DEST/'transfer_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print('smoke fields/diagnostics archived and verified')
def validate():
    import xml.etree.ElementTree as ET
    pp=sorted(CACHE.glob('*.vtu'));events=rows(DEST/'surface_events.csv');old=rows(ROOT/RUNS['B']/'surface_events.csv');assert len(pp)==11
    maximum=0.;pmaximum=0.;initial=None;profiles=[];storage_errors=[];budget=rows(DEST/'transport_budget.csv');par=params(DEST/'params.input');phi=float(par['Problem.Porosity'])
    for frame,p in enumerate(pp):
        a={v.attrib['Name']:np.fromstring(v.text or '',sep=' ') for v in ET.parse(p).iter('DataArray') if v.attrib.get('Name') in ['Coordinates','moleFraction','sw','pressure','pw']}
        assert 'pw' in a
        xyz=a['Coordinates'].reshape(-1,3);z=xyz[:,1];xs=np.unique(xyz[:,0]);zs=np.unique(z);wx=widths(xs);wz=widths(zs);order=np.argsort(np.searchsorted(zs,z)*len(xs)+np.searchsorted(xs,xyz[:,0]));xx=a['moleFraction'][order].reshape(len(zs),len(xs));sw=a['sw'][order].reshape(xx.shape);profiles.append(xx@wx/wx.sum())
        for row in [r for r in budget if int(r['step'])==frame]:
            mask=zs>=float(row['actual_z_over_H'])*.2
            salt0=phi*np.sum((sw*xx@wx/wx.sum())[mask]*wz[mask]);water0=phi*np.sum((sw@wx/wx.sum())[mask]*wz[mask]);storage_errors.append((abs(salt0-float(row['inventory_salt_m'])),abs(water0-float(row['inventory_water_m']))))
        for h in np.unique(z):
            maximum=max(maximum,float(np.ptp(a['moleFraction'][z==h])))
            for key in ['pressure','pw']:
                if key in a:pmaximum=max(pmaximum,float(np.ptp(a[key][z==h])))
        if '00000' in p.name:initial=a
    assert initial is not None
    initial_x_error=float(abs(initial['moleFraction']-float(par['Problem.SaltMoleFraction'])).max());xx0=float(par['Problem.SaltMoleFraction']);mass=xx0*float(par['Problem.MolarMassNaCl'])/((1-xx0)*float(par['Problem.MolarMassH2O'])+xx0*float(par['Problem.MolarMassNaCl']));rho=float(par['Problem.PureWaterDensity'])+1000*mass*(float(par['Problem.DensityLinearCoefficient'])+float(par['Problem.DensityQuadraticCoefficient'])*mass);z0=initial['Coordinates'].reshape(-1,3)[:,1];expected=1e5-(float(par['Problem.WaterLevel'])+z0)*rho*float(par['Problem.GravitationalAcceleration']);initial_p_error=float(abs(initial['pw']-expected).max())
    bdata=dict(np.load(ROOT/'tmp/wp8_zero_mode_20261103/B/profiles_spectra.npz'));profile_error=float(abs(np.array(profiles)-bdata['x_mean'][:11]).max());oldpar=params(ROOT/RUNS['B']/'params.input');changed=[k for k in par if par[k]!=oldpar.get(k)]
    oldbudget=rows(ROOT/RUNS['B']/'transport_budget.csv');errors={};reint=[]
    for key in ['inventory_salt_m','inventory_water_m','cut_salt_total_down_m_s','salt_residual_m','water_residual_m']:
        err=[]
        for region in {r['region'] for r in budget}:
            aa=[r for r in budget if r['region']==region];bb=[r for r in oldbudget if r['region']==region];pairs_=pairs([float(r['time_s']) for r in aa],[float(r['time_s']) for r in bb]);assert len(pairs_)==11
            err.extend(abs(float(aa[i][key])-float(bb[j][key])) for i,j in pairs_)
        errors[key]=max(err)
    for region in {r['region'] for r in budget}:
        cc=0.
        for row in [r for r in budget if r['region']==region]:
            if int(row['step']):cc+=float(row['dt_s'])*float(row['cut_salt_total_down_m_s'])
            reint.append(abs(cc-float(row['cumulative_cut_salt_m'])))
    stop=json.loads((DEST/'stop_record.json').read_text());execution=json.loads((DEST/'execution.json').read_text());salt=max(abs(float(r['salt_residual_m'])) for r in budget);water=max(abs(float(r['water_residual_m'])) for r in budget)
    checks=dict(exit=execution.get('returncode')==0,ten_steps=stop['accepted_steps']==10,requested_horizon=stop['reason']=='requested_horizon',exact_uniform_concentration=maximum==0,exact_uniform_pressure=pmaximum==0,salt_closure=salt<=1e-12,water_closure=water<=1e-10,inventory_bridge=errors['inventory_salt_m']<=1e-12,flux_bridge=errors['cut_salt_total_down_m_s']<=1e-17)
    checks.update(initial_x=initial_x_error<=1e-14,initial_pressure=initial_p_error<=1e-8,profile_bridge=profile_error<=1e-10,storage_reconstruction=max(x[0] for x in storage_errors)<=1e-14,flux_reintegration=max(reint)<=1e-14,inputs_preserved=set(changed)=={'Problem.Name','TimeLoop.TEnd'})
    checks={k:bool(v) for k,v in checks.items()}
    result=dict(passed=all(checks.values()),checks=checks,max_horizontal_x_range=maximum,max_horizontal_pressure_range=pmaximum,max_salt_residual_m=salt,max_water_residual_m=water,bridge_absolute_errors=errors,stop=stop,wall_s=execution['wall_s'],long_run_validated=False)
    result.update(initial_x_error=initial_x_error,initial_pressure_error_Pa=initial_p_error,profile_max_error=profile_error,storage_reconstruction_salt_error_m=max(x[0] for x in storage_errors),storage_reconstruction_water_error_m=max(x[1] for x in storage_errors),flux_reintegration_error_m=max(reint),changed_inputs=changed)
    (DEST/'validation.json').write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2));assert result['passed']
if __name__=='__main__':globals()[sys.argv[1]]()
