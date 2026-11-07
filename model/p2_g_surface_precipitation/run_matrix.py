import pathlib,subprocess,concurrent.futures,json
root=pathlib.Path.cwd()
base=(root/'source/model/p2_f_analytical_growth/integration/p2f_1d.input').read_text()
base=base.replace('p2f_analytical','surface_film').replace('TEnd = 60','TEnd = 1000')
cases=[('n40',40,.5,1e-6,True),('n80',80,.5,1e-6,True),('n160',160,.5,1e-6,True),
       ('dt_half',80,.25,1e-6,True),('threshold_half',80,.5,.5e-6,True),
       ('threshold_double',80,.5,2e-6,True),('no_precip_control',80,.5,1e-6,False)]
def run(case):
    name,n,dt,h,enabled=case
    path=root/'runs'/name
    path.mkdir(parents=True,exist_ok=False)
    (path/'params.input').write_text(base)
    args=[str(root/'build/solver'),'params.input','-Grid.Cells1',str(n),'-TimeLoop.DtInitial',str(dt),
          '-TimeLoop.MaxTimeStepSize',str(dt),'-P2G.ContinuumHeight',str(h),
          '-P2G.EnablePrecipitation',str(enabled).lower()]
    (path/'command.json').write_text(json.dumps(args,indent=2))
    with (path/'solver.out').open('w') as out:
        r=subprocess.run(args,cwd=path,stdout=out,stderr=subprocess.STDOUT)
    (path/'exit_status.txt').write_text(str(r.returncode))
    return name,r.returncode
with concurrent.futures.ThreadPoolExecutor(max_workers=7) as pool:
    outcomes=list(pool.map(run,cases))
(root/'results/outcomes.json').write_text(json.dumps(outcomes,indent=2))
if any(code for _,code in outcomes): raise SystemExit(2)
