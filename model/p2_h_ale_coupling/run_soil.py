"""Bounded independent cases; B activation-only precedes C coupled experiments."""
import concurrent.futures
import json
import pathlib
import subprocess
import sys

root=pathlib.Path.cwd()
stage=sys.argv[1]
base=(root/'repo/model/p2_f_analytical_growth/integration/p2f_1d.input').read_text()
base=base.replace('p2f_analytical','p2h_ale').replace('TEnd = 60','TEnd = 2500')
if stage=='b':
    cases=[('activation_1um',80,64,.5,1e-6,False),('activation_2um',80,64,.5,2e-6,False),
           ('activation_5um',80,64,.5,5e-6,False),('activation_10um',80,64,.5,1e-5,False),
           ('activation_rollback',80,64,.5,1e-6,True)]
elif stage=='c':
    cases=[('nominal',80,64,.5,1e-6,False),('soil160',160,64,.5,1e-6,False),
           ('crust128',80,128,.5,1e-6,False),('dt_half',80,64,.25,1e-6,False),
           ('threshold2',80,64,.5,2e-6,False),('threshold5',80,64,.5,5e-6,False),
           ('threshold10',80,64,.5,1e-5,False),('coupled_rollback',80,64,.5,1e-6,True)]
else:
    raise ValueError(stage)
def run(case):
    name,n,crust,dt,h,inject=case
    path=root/('stage_'+stage)/name
    path.mkdir(parents=True,exist_ok=False)
    (path/'params.input').write_text(base)
    args=[str(root/'build/soil_solver'),'params.input','-Grid.Cells1',str(n),
          '-TimeLoop.DtInitial',str(dt),'-TimeLoop.MaxTimeStepSize',str(dt),
          '-P2H.CrustCells',str(crust),'-P2H.ActivationHeight',str(h),
          '-P2H.ActivationOnly',str(stage=='b').lower(),'-P2H.InjectRejectedTrial',str(inject).lower(),
          '-P2H.RejectAfterActive',str(stage=='c').lower()]
    (path/'command.json').write_text(json.dumps(args,indent=2))
    with (path/'solver.out').open('w') as out:
        result=subprocess.run(args,cwd=path,stdout=out,stderr=subprocess.STDOUT)
    (path/'exit_status.txt').write_text(str(result.returncode))
    return dict(case=name,exit_code=result.returncode)
with concurrent.futures.ThreadPoolExecutor(max_workers=4) as pool:
    results=list(pool.map(run,cases))
(root/('stage_'+stage)/'outcomes.json').write_text(json.dumps(results,indent=2))
print(json.dumps(results,indent=2))
# Physical supply failures remain labeled outcomes; the audit decides success.
