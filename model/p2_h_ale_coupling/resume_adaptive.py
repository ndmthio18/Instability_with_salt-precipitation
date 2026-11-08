"""Finish only cases interrupted/unstarted by the primary Slurm TIMEOUT.

Never retry a numerical failure automatically, and never overwrite a prefix.
Commands must match the primary runner's recorded argument arrays exactly.
"""
import concurrent.futures,csv,json,os,pathlib,subprocess
assert os.environ.get('SLURM_JOB_ID')
root=pathlib.Path.cwd();folder=root/'matrix'
cases=json.loads((folder/'cases.json').read_text())
base=(root/'repo/model/p2_f_analytical_growth/integration/p2f_1d.input').read_text()
base=base.replace('p2f_analytical','p2h_ale').replace('TEnd = 60','TEnd = 2500')
def args(c):
    options={'Grid.Cells1':c['n'],'TimeLoop.DtInitial':c.get('dt',.5),
        'TimeLoop.MaxTimeStepSize':c.get('dt',.5),'TimeLoop.TEnd':c.get('tend',2500),
        'P2H.DonorMode':c['mode'],'P2H.MaxDonorDepth':c['depth_mm']*1e-3,
        'P2H.ActivationHeight':c['h_um']*1e-6,'P2H.CrustCells':64,
        'P2H.ControlledDonorTest':c.get('controlled',False),'P2H.PrecrustOnly':c.get('precrust',False),
        'P2H.InjectRejectedTrial':c.get('inject',False),'P2H.RejectAfterActive':c.get('reject_active',False)}
    if 'grading' in c:options['Grid.Grading1']=c['grading']
    command=[str(root/'build/soil_solver'),'params.input']
    for k,v in options.items():command+=['-'+k,str(v).lower() if isinstance(v,bool) else str(v)]
    return command
for c in cases:
    p=folder/c['name']
    if (p/'command.json').exists():
        assert json.loads((p/'command.json').read_text())==args(c)
        assert (p/'params.input').read_text()==base
state=subprocess.check_output(['sacct','-j','16371535','--noheader','--parsable2','--format=State']).decode().splitlines()[0].split('|')[0].strip()
metadata=dict(primary_job='16371535',primary_state=state,job=os.environ['SLURM_JOB_ID'],command_identity_checked=True,cases=[])
if state=='COMPLETED':
    metadata['action']='no continuation needed'
elif state=='TIMEOUT':
    pending=[]
    for c in cases:
        p=(folder/c['name']).resolve();assert p.parent==folder.resolve()
        code=int((p/'exit_status.txt').read_text()) if (p/'exit_status.txt').exists() else None
        if code is not None and code not in (-15,-9,137,143):continue
        item=dict(case=c['name'],previous_exit_code=code,last_accepted_time_s=None)
        if p.exists():
            if (p/'coupled_history.csv').exists():
                with (p/'coupled_history.csv').open() as f:
                    for row in csv.DictReader(f):
                        if all(v is not None and v!='' for v in row.values()):item['last_accepted_time_s']=float(row['time_s'])
            dest=root/'failures/16371535'/c['name'];dest.parent.mkdir(parents=True,exist_ok=True)
            assert not dest.exists();p.rename(dest);item['archived_prefix']=str(dest.relative_to(root))
        metadata['cases'].append(item);pending.append(c)
    # Start the potentially longest remaining timestep/rollback checks first;
    # case physics and commands are identical regardless of placement order.
    pending.sort(key=lambda c:(0 if c.get('dt',.5)<.5 else 1 if c.get('inject') else 2,c['name']))
    metadata['execution_order']=[c['name'] for c in pending]
    def run(c):
        p=folder/c['name'];p.mkdir(exist_ok=False);(p/'params.input').write_text(base)
        command=args(c);(p/'command.json').write_text(json.dumps(command,indent=2))
        with (p/'solver.out').open('w') as f:
            try:code=subprocess.run(command,cwd=p,stdout=f,stderr=subprocess.STDOUT,timeout=1800).returncode
            except subprocess.TimeoutExpired:code=124
        (p/'exit_status.txt').write_text(str(code));print(json.dumps(dict(case=c['name'],exit_code=code)),flush=True)
    (root/'resume_metadata.json').write_text(json.dumps(metadata,indent=2))
    with concurrent.futures.ThreadPoolExecutor(max_workers=int(os.environ['SLURM_CPUS_PER_TASK'])) as pool:list(pool.map(run,pending))
    outcomes=[dict(case=c['name'],exit_code=int((folder/c['name']/'exit_status.txt').read_text())) for c in cases]
    (folder/'outcomes.json').write_text(json.dumps(outcomes,indent=2));metadata['action']='replayed scheduler-interrupted cases from frozen initial conditions'
else:raise RuntimeError('Unexpected primary job state; diagnose without automatic replay: '+state)
(root/'resume_metadata.json').write_text(json.dumps(metadata,indent=2));print(json.dumps(metadata,indent=2))
