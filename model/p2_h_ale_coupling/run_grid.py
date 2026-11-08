"""Real DuMux runs, immutable case directories, explicit physical meshes."""
import concurrent.futures, json, os, pathlib, subprocess, sys, time
root=pathlib.Path.cwd(); stage=sys.argv[1]
assert os.environ.get('SLURM_JOB_ID')
base=(root/'repo/model/p2_f_analytical_growth/integration/p2f_1d.input').read_text()
base=base.replace('p2f_analytical','p2h_ale')
base=base.split('[P2F]')[0]+'[Transport]'+base.split('[Transport]',1)[1]
cases=[]
def add(name,n,**kw): cases.append(dict(name=name,n=n,**kw))
if stage=='fixed_smoke':
    add('fixed320_20s',320,tend=20,mode='implicit_surface')
    add('projection80_20s',80,tend=20,mode='projection')
elif stage=='fixed_tests':
    for mode in ('projection','implicit_surface'):
        add('undersaturated_'+mode,80,tend=20,mode=mode,salt=.01)
    add('projection40_regression',40,mode='projection')
elif stage=='fixed_pre':
    for n in (40,80,160,320): add(f'fixed{n}',n,mode='implicit_surface')
    for n in (160,320,640): add(f'fixedA{n}',n,grading=-.97**(80/n),mode='implicit_surface')
    for deep,fine in ((40,80),(80,80),(160,40),(160,80),(160,160),(320,80)):
        add(f'fixedB{deep}_top{fine}',f'{deep} {fine}',positions='0 0.198 0.2',grading='1 1',mode='implicit_surface')
    for n in (160,320): add(f'fixed_half{n}',n,dt=.25,mode='implicit_surface')
elif stage=='fixed_coupled':
    for n in (40,80,160,320): add(f'fixed{n}',n,pre=False,mode='implicit_surface')
    for n in (160,320): add(f'fixed_half{n}',n,dt=.25,pre=False,mode='implicit_surface')
    add('fixed80_rollback',80,pre=False,mode='implicit_surface',inject=True)
elif stage=='smoke':
    add('historical320_10s',320,tend=10)
    add('familyB120_10s','40 80',positions='0 0.198 0.2',grading='1 1',tend=10)
elif stage=='diagnostic':
    for n in (40,80,160,320): add(f'historical{n}',n)
    for n in (40,80,160,320,640): add(f'familyA{n}',n,grading=-.97**(80/n))
    for deep in (40,80,160,320):
        add(f'familyB{deep}_top80',f'{deep} 80',positions='0 0.198 0.2',grading='1 1')
    for fine in (40,160): add(f'familyB160_top{fine}',f'160 {fine}',positions='0 0.198 0.2',grading='1 1')
    add('historical320_half',320,dt=.25)
elif stage=='ale_diagnostic':
    add('historical320_ale',320,pre=False)
    add('historical160_ale',160,pre=False,tend=710)
else: raise ValueError(stage)
folder=root/stage; folder.mkdir(exist_ok=True)
(folder/'cases.json').write_text(json.dumps(cases,indent=2))
def run(c):
    d=folder/c['name']; d.mkdir(exist_ok=False)
    (d/'params.input').write_text(base)
    options={'Grid.Cells1':c['n'],'TimeLoop.DtInitial':c.get('dt',.5),
        'TimeLoop.MaxTimeStepSize':c.get('dt',.5),'TimeLoop.TEnd':c.get('tend',2500),
        'P2H.DonorMode':'adaptive_multilayer','P2H.MaxDonorDepth':.005,
        'P2H.CrustCells':64,'P2H.PrecrustOnly':c.get('pre',True),'P2H.GridDiagnostics':True,
        'P2H.PrecipitationMode':c.get('mode','projection'),
        'P2H.InjectRejectedTrial':c.get('inject',False)}
    if 'salt' in c: options['Problem.SaltMoleFraction']=c['salt']
    if 'grading' in c: options['Grid.Grading1']=c['grading']
    if 'positions' in c: options['Grid.Positions1']=c['positions']
    args=[str(root/'build/soil_solver'),'params.input']
    for k,v in options.items(): args+=['-'+k,str(v).lower() if isinstance(v,bool) else str(v)]
    (d/'command.json').write_text(json.dumps(args,indent=2)); start=time.monotonic()
    with (d/'solver.out').open('w') as out:
        try: code=subprocess.run(args,cwd=d,stdout=out,stderr=subprocess.STDOUT,timeout=5000).returncode
        except subprocess.TimeoutExpired: code=124
    (d/'exit_status.txt').write_text(str(code))
    result=dict(case=c['name'],exit_code=code,elapsed_s=time.monotonic()-start)
    print(json.dumps(result),flush=True); return result
with concurrent.futures.ThreadPoolExecutor(max_workers=int(os.environ.get('SLURM_CPUS_PER_TASK',1))) as pool:
    results=list(pool.map(run,cases))
(folder/'outcomes.json').write_text(json.dumps(results,indent=2))
