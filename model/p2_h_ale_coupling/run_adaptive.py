"""Registered matrix; run only under the supplied Slurm scripts."""
import concurrent.futures,json,os,pathlib,subprocess,sys
root=pathlib.Path.cwd(); stage=sys.argv[1]
assert os.environ.get('SLURM_JOB_ID')
base=(root/'repo/model/p2_f_analytical_growth/integration/p2f_1d.input').read_text()
base=base.replace('p2f_analytical','p2h_ale').replace('TEnd = 60','TEnd = 2500')
cases=[]
def add(name,n=80,h=1,depth=5,mode='adaptive_multilayer',**extra):
    cases.append(dict(name=name,n=n,h_um=h,depth_mm=depth,mode=mode,**extra))
if stage=='smoke':
    for n in (80,160): add(f'smoke_n{n}',n,tend=400)
elif stage=='controlled':
    for n in (40,80,160,320):
        for h in (1,2,5,10):
            for d in (2,5,10): add(f'control_n{n}_h{h}_d{d}',n,h,d,controlled=True)
    for n in (80,160): add(f'uniform_n{n}',n,5,5,controlled=True,grading=1)
    add('insufficient_depth',80,1,.000001,controlled=True)
elif stage=='matrix':
    for n in (40,80,160,320):
        add(f'old_n{n}',n,mode='top_only')
        add(f'pre_n{n}',n,precrust=True)
        for h in (1,2,5,10):
            for d in (2,5,10): add(f'new_n{n}_h{h}_d{d}',n,h,d)
    for h in (2,5,10): add(f'old_n80_h{h}',80,h,mode='top_only')
    for n in (160,320): add(f'dt_half_n{n}',n,dt=.25)
    add('rollback_activation',inject=True)
    add('rollback_coupled',inject=True,reject_active=True)
else: raise ValueError(stage)
folder=root/stage; folder.mkdir(exist_ok=True)
(folder/'cases.json').write_text(json.dumps(cases,indent=2))
def run(c):
    path=folder/c['name']; path.mkdir(exist_ok=False)
    (path/'params.input').write_text(base)
    options={'Grid.Cells1':c['n'],'TimeLoop.DtInitial':c.get('dt',.5),
        'TimeLoop.MaxTimeStepSize':c.get('dt',.5),'TimeLoop.TEnd':c.get('tend',2500),
        'P2H.DonorMode':c['mode'],'P2H.MaxDonorDepth':c['depth_mm']*1e-3,
        'P2H.ActivationHeight':c['h_um']*1e-6,'P2H.CrustCells':64,
        'P2H.ControlledDonorTest':c.get('controlled',False),
        'P2H.PrecrustOnly':c.get('precrust',False),
        'P2H.InjectRejectedTrial':c.get('inject',False),'P2H.RejectAfterActive':c.get('reject_active',False)}
    if 'grading' in c: options['Grid.Grading1']=c['grading']
    args=[str(root/'build/soil_solver'),'params.input']
    for k,v in options.items(): args+=['-'+k,str(v).lower() if isinstance(v,bool) else str(v)]
    (path/'command.json').write_text(json.dumps(args,indent=2))
    with (path/'solver.out').open('w') as out:
        try: code=subprocess.run(args,cwd=path,stdout=out,stderr=subprocess.STDOUT,timeout=1800).returncode
        except subprocess.TimeoutExpired: code=124
    (path/'exit_status.txt').write_text(str(code))
    result=dict(case=c['name'],exit_code=code)
    print(json.dumps(result),flush=True); return result
with concurrent.futures.ThreadPoolExecutor(max_workers=int(os.environ.get('SLURM_CPUS_PER_TASK','1'))) as pool:
    results=list(pool.map(run,cases))
(folder/'outcomes.json').write_text(json.dumps(results,indent=2))
