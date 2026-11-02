"""Static/provenance checks only. Does not load or execute the C++ program."""
import ast,difflib,hashlib,json,re,subprocess
from datetime import datetime,timezone
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];MODEL=ROOT/'model/transport_v1';OUT=ROOT/'analysis/transport_build_20261102'
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
up=json.loads((MODEL/'upstream_manifest.json').read_text())
shared=json.loads((OUT/'shared_source_after_build.json').read_text())
latest=sorted(OUT.glob('attempt_*/result.json'))[-1];build=json.loads(latest.read_text())
src=MODEL/'source';app=src/'appl/simpleRichardsnc'
main=(app/'main.cc').read_text(encoding='utf-8');model=(app/'model.hh').read_text(encoding='utf-8');diag=(app/'transportdiagnostics.hh').read_text()
checks={}
checks['upstream_snapshot_bytes_intact']=all(sha(MODEL/'upstream'/n)==v['sha256'] for n,v in up['files'].items())
checks['shared_remote_source_unchanged']=all(shared.get(n)==v['sha256'] for n,v in up['files'].items())
actual={str(p.relative_to(src)).replace('\\','/'):sha(p) for p in src.rglob('*') if p.is_file()}
checks['compiled_source_matches_local']=actual==build['source_sha256']
checks['configure_and_link_success']=build['configure_exit']==0 and build.get('build_exit')==0
checks['no_application_execution_in_build']=build['new_simulations'] is False and all(c[0]=='cmake' for c in build['commands'])
checks['shared_flux_and_storage_kernels']=all(x in model for x in ['return normalizedStorage(problem, volVars);','const auto f = transportFlux(','static TransportFlux transportFlux']) and all(x in diag for x in ['Residual::normalizedStorage(','Residual::transportFlux('])
checks['storage_stays_normalized']='storage[Indices::conti1EqIdx] = problem.porosity() * volVars.sw() * volVars.moleFraction();' in model
checks['accepted_dt_integration_and_order']=main.index('const Scalar acceptedDt = timeLoop->timeStepSize();')>main.index('nonLinearSolver.solve(x, *timeLoop);') and 'transport.write(timeLoop->time(), acceptedDt, x);' in main
checks['serial_and_output_overwrite_guards']='comm().size() != 1' in diag and 'Refusing to overwrite diagnostic:' in diag
checks['threshold_stop_and_interval']='!timeLoop->finished() && !transport.shouldStop()' in main and 'surfaceLow_=initialized_ ? lastTime_ : 0' in diag
checks['darcy_visualization_includes_krw']='problem->krw(volVars.sw()) / problem->viscosity()' in main and main.index('updateCellVelocity();')<main.index('vtkWriter.write(0.0);')
checks['bottom_reaction_explicitly_labeled']='not an independent boundary-flux check' in diag
pyfiles=[ROOT/'scripts'/n for n in ['prepare_transport_source.py','build_transport_remote.py','collect_transport_build.py','prepare_transport_run.py','check_transport_static.py','frodo_bridge.py']]
for p in pyfiles:ast.parse(p.read_text(encoding='utf-8'),filename=str(p))
checks['python_syntax']=True
# Reproducible delta from the preserved current upstream snapshot.
patch=[]
for name in sorted(set(up['files'])|set(actual)):
 a=MODEL/'upstream'/name;b=src/name
 aa=a.read_text(encoding='utf-8').splitlines(keepends=True) if a.exists() else []
 bb=b.read_text(encoding='utf-8').splitlines(keepends=True) if b.exists() else []
 patch.extend(difflib.unified_diff(aa,bb,fromfile='a/'+name if a.exists() else '/dev/null',tofile='b/'+name if b.exists() else '/dev/null'))
patchpath=MODEL/'transport_v1.patch';patchpath.write_text(''.join(patch),encoding='utf-8',newline='\n')
# git apply --check is read-only, but operates relative to a temporary copy of
# the baseline so the parent repository discovery cannot redirect its scope.
import tempfile,shutil
with tempfile.TemporaryDirectory(prefix='paper2_patch_check_') as d:
 shutil.copytree(MODEL/'upstream',Path(d)/'source')
 r=subprocess.run(['git','apply','--check',str(patchpath.resolve())],cwd=Path(d)/'source',capture_output=True,text=True)
checks['patch_applies_to_preserved_baseline']=r.returncode==0
log=(latest.parent/'build.log').read_text(encoding='utf-8',errors='replace')
warnings=[line for line in log.splitlines() if 'warning:' in line]
report={'checked_utc':datetime.now(timezone.utc).isoformat(),'checks':checks,'passed':all(checks.values()),
 'patch_check_stderr':r.stderr,'source_sha256':actual,'patch_sha256':sha(patchpath),
 'warning_count':len(warnings),'diagnostics_header_warnings':[x for x in warnings if 'transportdiagnostics.hh:' in x],
 'project_warnings':[x for x in warnings if '/paper2_transport_20261102/source/' in x],
 'scope':'static assertions and provenance only; no runtime conservation, convergence, or feasibility claim',
 'simulation_started':False,'build_result_sha256':sha(latest)}
(OUT/'static_checks.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
print(json.dumps({'passed':report['passed'],'checks':checks,'warning_count':len(warnings),'diagnostics_header_warnings':report['diagnostics_header_warnings']},indent=2))
if not report['passed']:raise SystemExit(1)
