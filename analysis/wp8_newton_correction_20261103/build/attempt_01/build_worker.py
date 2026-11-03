
import json,subprocess,hashlib,datetime,sys,os
from pathlib import Path
root=Path('/temp/zhixin/paper2_wp8_background_v2_build_20261103');attempt=Path(sys.argv[1])
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
 report['build_exit']=cmd(['cmake','--build',str(build),'--target','paper2_transport','--parallel','1'],'build.log')
 exe=build/'appl/simpleRichardsnc/paper2_transport'
 if report['build_exit']==0:
  report['executable']={'path':str(exe),'size':exe.stat().st_size,'sha256':hashlib.sha256(exe.read_bytes()).hexdigest()}
  report['elf_dynamic']=capture(['readelf','-d',str(exe)])
  report['compile_commands']=json.loads((build/'compile_commands.json').read_text())
report['source_sha256']={str(p.relative_to(root/'source')):hashlib.sha256(p.read_bytes()).hexdigest() for p in (root/'source').rglob('*') if p.is_file()}
report['finished_utc']=datetime.datetime.now(datetime.timezone.utc).isoformat()
(attempt/'result.json').write_text(json.dumps(report,indent=2))
