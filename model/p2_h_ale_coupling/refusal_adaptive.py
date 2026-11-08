"""Real transport refusal vs precrust control; no partial donor mutation allowed."""
import concurrent.futures,json,os,pathlib,subprocess
assert os.environ.get('SLURM_JOB_ID')
root=pathlib.Path.cwd();folder=root/'refusal';folder.mkdir(exist_ok=False)
base=(root/'repo/model/p2_f_analytical_growth/integration/p2f_1d.input').read_text().replace('p2f_analytical','p2h_ale')
cases=[dict(name='depth_refused',n=80,h_um=1,depth_mm=.000001,mode='adaptive_multilayer',precrust=False,tend=400),
       dict(name='precrust_reference',n=80,h_um=1,depth_mm=.000001,mode='adaptive_multilayer',precrust=True,tend=400)]
(folder/'cases.json').write_text(json.dumps(cases,indent=2))
def run(c):
    path=folder/c['name'];path.mkdir();(path/'params.input').write_text(base)
    args=[str(root/'build/soil_solver'),'params.input','-TimeLoop.TEnd','400',
        '-TimeLoop.DtInitial','.5','-TimeLoop.MaxTimeStepSize','.5','-P2H.DonorMode','adaptive_multilayer',
        '-P2H.MaxDonorDepth','1e-9','-P2H.PrecrustOnly',str(c['precrust']).lower()]
    (path/'command.json').write_text(json.dumps(args,indent=2))
    with (path/'solver.out').open('w') as f:code=subprocess.run(args,cwd=path,stdout=f,stderr=subprocess.STDOUT,timeout=300).returncode
    (path/'exit_status.txt').write_text(str(code));return dict(case=c['name'],exit_code=code)
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:outcomes=list(pool.map(run,cases))
(folder/'outcomes.json').write_text(json.dumps(outcomes,indent=2));print(outcomes)
