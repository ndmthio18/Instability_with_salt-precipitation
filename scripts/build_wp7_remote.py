"""Stage only Paper 2's isolated source and start configure/build, never the app.

The remote worker records commands, toolchain and binary identity. Each invocation
keeps a separate attempt directory. No ctest, solver launch or shared-source edit.
"""
import base64, hashlib, json, subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
SOURCE=ROOT/'model/transport_wp7_localflux_v1/source'
files={str(p.relative_to(SOURCE)).replace('\\','/'):{'data':base64.b64encode(p.read_bytes()).decode(),
       'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in SOURCE.rglob('*') if p.is_file()}
worker=r'''
import json,subprocess,hashlib,datetime,sys,os
from pathlib import Path
root=Path('/temp/zhixin/paper2_wp7_localflux_build_20261103');attempt=Path(sys.argv[1])
report={'new_simulations':False,'commands':[],'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
def cmd(args,name):
 report['commands'].append(args)
 with (attempt/name).open('w') as f:
  p=subprocess.run(args,stdout=f,stderr=subprocess.STDOUT,env=dict(os.environ,OPENBLAS_NUM_THREADS='1',OMP_NUM_THREADS='1'))
 return p.returncode
def capture(args):
 p=subprocess.run(args,capture_output=True,text=True);return {'exit_code':p.returncode,'stdout':p.stdout,'stderr':p.stderr}
stack=Path('/home/zhixin/dumux/dumux')
deps=[p for p in stack.iterdir() if p.is_dir() and (p.name.startswith('dune-') or p.name=='dumux') and (p/'build-cmake').is_dir()]
report['dependencies']={str(p):{'git':capture(['git','-C',str(p),'rev-parse','HEAD']),
 'module':(p/'dune.module').read_text() if (p/'dune.module').exists() else ''} for p in deps}
report['compiler']=capture(['/usr/bin/c++','--version']);report['cmake']=capture(['cmake','--version'])
prefix=';'.join(str(p/'build-cmake') for p in deps)
build=root/'build-gcc-release'
configure=['cmake','-S',str(root/'source'),'-B',str(build),'-DCMAKE_BUILD_TYPE=Release',
 '-DCMAKE_CXX_COMPILER=/usr/bin/c++','-DCMAKE_EXPORT_COMPILE_COMMANDS=ON','-DCMAKE_PREFIX_PATH='+prefix]
report['configure_exit']=cmd(configure,'configure.log')
if report['configure_exit']==0:
 report['build_exit']=cmd(['cmake','--build',str(build),'--target','paper2_transport','--parallel','2'],'build.log')
 exe=build/'appl/simpleRichardsnc/paper2_transport'
 if report['build_exit']==0:
  report['executable']={'path':str(exe),'size':exe.stat().st_size,'sha256':hashlib.sha256(exe.read_bytes()).hexdigest()}
  report['elf_dynamic']=capture(['readelf','-d',str(exe)])
  report['compile_commands']=json.loads((build/'compile_commands.json').read_text())
report['source_sha256']={str(p.relative_to(root/'source')):hashlib.sha256(p.read_bytes()).hexdigest() for p in (root/'source').rglob('*') if p.is_file()}
report['finished_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
(attempt/'result.json').write_text(json.dumps(report,indent=2))
'''
remote=r'''
import json,base64,hashlib,subprocess,sys
from pathlib import Path
root=Path('/temp/zhixin/paper2_wp7_localflux_build_20261103');source=root/'source'
root.mkdir(exist_ok=True)
assert not (root/'build_started.json').exists(), 'WP7 build already launched; collect first'
import os
mem=int(next(x.split()[1] for x in Path('/proc/meminfo').read_text().splitlines() if x.startswith('MemAvailable:')))*1024
assert os.getloadavg()[0]<8 and mem>8*1024**3
(root/'build_started.json').write_text(json.dumps(dict(load=os.getloadavg(),memory=mem)))
if root.resolve()!=Path('/temp/zhixin/paper2_wp7_localflux_build_20261103'):raise RuntimeError('Unexpected build root')
attempt=root/('attempt_%02d'%(len(list(root.glob('attempt_*')))+1));attempt.mkdir()
files=FILES
for name,rec in files.items():
 p=source/name
 if not p.resolve().is_relative_to(source.resolve()):raise ValueError(name)
 b=base64.b64decode(rec['data']);assert hashlib.sha256(b).hexdigest()==rec['sha256']
 p.parent.mkdir(parents=True,exist_ok=True);p.write_bytes(b)
worker=attempt/'build_worker.py';worker.write_text(WORKER)
(attempt/'staged_source_sha256.json').write_text(json.dumps({k:v['sha256'] for k,v in files.items()},indent=2))
log=(attempt/'worker.log').open('w')
p=subprocess.Popen([sys.executable,str(worker),str(attempt)],stdout=log,stderr=subprocess.STDOUT,start_new_session=True,stdin=subprocess.DEVNULL)
print(json.dumps({'attempt':str(attempt),'pid':p.pid,'action':'configure_and_build_only'}))
'''.replace('FILES',repr(files)).replace('WORKER',repr(worker))
bridge='/mnt/d/3OutputFile/1Research/Article/22Surrogate_model/paper/2/scripts/frodo_bridge.py'
r=subprocess.run(['wsl','python3',bridge,'run','python3 -'],input=remote.encode(),capture_output=True,timeout=60)
if r.returncode:raise RuntimeError(r.stderr.decode(errors='replace'))
print(r.stdout.decode())
out=ROOT/'analysis/wp7_shortwave_validation_20261103/build';out.mkdir(parents=True,exist_ok=True)
j=json.loads(r.stdout);(out/(Path(j['attempt']).name+'_launch.json')).write_text(json.dumps(j,indent=2),encoding='utf-8')
