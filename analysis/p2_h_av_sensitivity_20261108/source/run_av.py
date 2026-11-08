"""Registered P2-H av sensitivity matrix; Slurm execution only."""
import concurrent.futures,json,os,pathlib,subprocess,sys,time

root=pathlib.Path.cwd()
stage=sys.argv[1]
assert os.environ.get('SLURM_JOB_ID')
base=(root/'repo/model/p2_f_analytical_growth/integration/p2f_1d.input').read_text()
base=base.replace('p2f_analytical','p2h_ale')
base=base.split('[P2F]')[0]+'[Transport]'+base.split('[Transport]',1)[1]
cases=[]
def add(name,n,av,**extra):
    cases.append(dict(name=name,n=n,av=av,ale=extra.pop('ale',64),dt=extra.pop('dt',.5),**extra))
if stage in ('smoke','smoke_r02','smoke_exact'):
    for av in (.1,1.,10.): add(f'smoke_n80_av{av:g}',80,av,tend=720)
    add('smoke_n80_av1_rollback',80,1.,tend=700,inject=True)
elif stage in ('primary','primary_exact'):
    for n in (80,160):
        for av in (.1,1.,10.): add(f'n{n}_av{av:g}',n,av)
elif stage in ('dt_half','dt_half_exact'):
    for n in (80,160):
        for av in (.1,1.,10.): add(f'n{n}_av{av:g}_dt2',n,av,dt=.25)
elif stage in ('ale_refine','ale_refine_exact'):
    for av in (.1,1.,10.): add(f'n80_av{av:g}_ale128',80,av,ale=128)
elif stage == 'ale256_exact':
    add('n80_av10_ale256',80,10.,ale=256)
else:
    raise ValueError(stage)
folder=root/stage
folder.mkdir(exist_ok=False)
(folder/'cases.json').write_text(json.dumps(cases,indent=2))
def run(case):
    path=folder/case['name']
    path.mkdir()
    (path/'params.input').write_text(base)
    options={
        'Grid.Cells1':case['n'],
        'TimeLoop.DtInitial':case['dt'],
        'TimeLoop.MaxTimeStepSize':case['dt'],
        'TimeLoop.TEnd':case.get('tend',2500),
        'P2H.ActivationHeight':3.5e-6,
        'P2H.DonorMode':'adaptive_multilayer',
        'P2H.MaxDonorDepth':.005,
        'P2H.MinimumSupplySaturation':.99,
        'P2H.CrustCells':case['ale'],
        'P2H.ReactionCoefficient':1e-3,
        'P2H.AvMultiplier':case['av'],
        'P2H.PrecipitationMode':'implicit_surface',
        'P2H.GridDiagnostics':True,
        'P2H.InjectRejectedTrial':case.get('inject',False),
    }
    build='build_exact' if stage.endswith('exact') else 'build'
    args=[str(root/build/'soil_solver'),'params.input']
    for key,value in options.items():
        args += ['-'+key,str(value).lower() if isinstance(value,bool) else str(value)]
    (path/'command.json').write_text(json.dumps(args,indent=2))
    started=time.monotonic()
    with (path/'solver.out').open('w') as output:
        try:
            code=subprocess.run(args,cwd=path,stdout=output,stderr=subprocess.STDOUT,timeout=5400).returncode
        except subprocess.TimeoutExpired:
            code=124
    result=dict(case=case['name'],exit_code=code,elapsed_s=time.monotonic()-started)
    (path/'exit_status.txt').write_text(str(code))
    print(json.dumps(result),flush=True)
    return result
with concurrent.futures.ThreadPoolExecutor(max_workers=min(len(cases),int(os.environ.get('SLURM_CPUS_PER_TASK','1')))) as pool:
    outcomes=list(pool.map(run,cases))
(folder/'outcomes.json').write_text(json.dumps(outcomes,indent=2))
