"""Registered calculations and gates. Only invoke in an HPCC Slurm allocation."""
import csv,hashlib,json,os,pathlib,subprocess,sys,time
root=pathlib.Path.cwd(); src=root/'repo/model/prat_analytical_boundary'
if not os.environ.get('SLURM_JOB_ID') or os.uname().nodename.startswith(('ibcln','login','test')): raise RuntimeError('compute allocation required')
phase=sys.argv[1]
def run(name,options):
    d=root/'runs'/name; d.mkdir(parents=True,exist_ok=False)
    (d/'params.input').write_text((src/'column.input').read_text())
    args=[str(root/'build/soil_solver'),'params.input']
    for k,v in options.items(): args+=['-'+k,str(v).lower() if isinstance(v,bool) else str(v)]
    (d/'command.json').write_text(json.dumps(args,indent=2))
    start=time.monotonic()
    with (d/'solver.out').open('w') as log:
        try: code=subprocess.run(args,cwd=str(d),stdout=log,stderr=subprocess.STDOUT,timeout=16*3600).returncode
        except subprocess.TimeoutExpired: code=124
    rows=list(csv.DictReader((d/'history.csv').open())) if (d/'history.csv').exists() else []
    last=rows[-1] if rows else {}; horizon=float(options.get('TimeLoop.TEnd',25000))
    result={'name':name,'exit_code':code,'job':os.environ['SLURM_JOB_ID'],'node':os.uname().nodename,'elapsed_s':time.monotonic()-start,
      'requested_s':horizon,'last':last,'max_salt_residual':max([abs(float(x['salt_residual_rel'])) for x in rows] or [1]),
      'max_water_residual':max([abs(float(x['water_residual_rel'])) for x in rows] or [1])}
    result['completed']=code==0 and abs(float(last.get('time_s',-1))-horizon)<1e-6
    result['coupling_pass']=result['completed'] and last.get('active')=='1' and max(result['max_salt_residual'],result['max_water_residual'])<=1e-8
    (d/'outcome.json').write_text(json.dumps(result,indent=2)); print(json.dumps(result),flush=True)
    return result
if phase=='verification':
    hydro=run('hydrostatic',{'TimeLoop.TEnd':100,'PratBoundary.ZeroEvaporation':True})
    short=run('short_evaporation',{'TimeLoop.TEnd':300})
    assert hydro['completed'] and short['completed']
    assert float(hydro['last']['h_m'])==0 and abs(float(hydro['last']['vapour_kg']))==0
    assert max(hydro['max_salt_residual'],hydro['max_water_residual'],short['max_salt_residual'],short['max_water_residual'])<=1e-8
    (root/'unit/verification.json').write_text(json.dumps({'status':'PASS','cases':[hydro,short]},indent=2))
elif phase=='pilot':
    assert json.loads((root/'unit/verification.json').read_text())['status']=='PASS'
    run('pilot_n80_dt0p5',{})
elif phase=='matrix':
    pilot=json.loads((root/'runs/pilot_n80_dt0p5/outcome.json').read_text())
    if not pilot['coupling_pass']: raise RuntimeError('pilot gate failed; matrix forbidden')
    index=int(sys.argv[2]); n,dt=[(40,.5),(160,.5),(320,.5),(80,.25)][index]
    with (root/'runs/pilot_n80_dt0p5/soil_profiles.csv').open() as f:
        anchor=sorted(set(float(r['z_m']) for r in csv.DictReader(f) if float(r['time_s'])==0))
    assert len(anchor)==81
    if n==40: nodes=anchor[::2]
    else:
        factor=n//80
        nodes=[a+(b-a)*i/factor for a,b in zip(anchor[:-1],anchor[1:]) for i in range(factor)]+[anchor[-1]]
    opts={'Grid.Positions1':' '.join(format(z,'.17g') for z in nodes),'Grid.Cells1':' '.join(['1']*n),'Grid.Grading1':' '.join(['1']*n),
          'TimeLoop.DtInitial':dt,'TimeLoop.MaxTimeStepSize':dt}
    run('nested_n{}_dt{}'.format(n,str(dt).replace('.','p')),opts)
elif phase=='rollback':
    # Same physical initial/reservoir composition; replay to first activation + 1 s.
    pilot=json.loads((root/'runs/pilot_n80_dt0p5/outcome.json').read_text())
    activation=float(pilot['last']['activation_s']); assert activation>0
    run('rollback_injection',{'TimeLoop.TEnd':activation+1,'PratBoundary.InjectRejectedTrial':True})
elif phase=='regression':
    executable=pathlib.Path('/WORK/hulm_work/zhixin/p2_h_progressive_activation_20261108/build/soil_solver')
    d=root/'runs/frozen_short'; d.mkdir(parents=True,exist_ok=False)
    (d/'params.input').write_text((src/'column.input').read_text())
    args=[str(executable),'params.input','-Problem.EvaporationBoundaryMode','p2h_ale',
          '-P2H.DonorMode','progressive_top_boundary_filling','-P2H.PrecipitationMode','implicit_surface',
          '-P2H.PrecrustOnly','true','-TimeLoop.TEnd','300']
    (d/'command.json').write_text(json.dumps(args,indent=2))
    (d/'binary.sha256').write_text(hashlib.sha256(executable.read_bytes()).hexdigest()+'  '+str(executable)+'\n')
    with (d/'solver.out').open('w') as f: code=subprocess.run(args,cwd=str(d),stdout=f,stderr=subprocess.STDOUT).returncode
    assert code==0
    old=list(csv.DictReader((d/'coupled_history.csv').open()))
    new=list(csv.DictReader((root/'runs/short_evaporation/history.csv').open()))
    assert len(old)==len(new)
    mapping={'time_s':'time_s','c_soil':'top_soil_mass_fraction','sw_soil':'top_soil_sw',
             'j_kg_m2_s':'evaporation_kg_m2_s','soil_salt_kg':'soil_salt_kg','soil_water_kg':'soil_water_kg'}
    errors={a:max(abs(float(x[a])-float(y[b])) for x,y in zip(new[1:],old[1:])) for a,b in mapping.items()}
    assert max(errors.values())<1e-12
    (root/'unit/frozen_regression.json').write_text(json.dumps({'status':'PASS','exit_code':code,'max_absolute_differences':errors,
         'job':os.environ['SLURM_JOB_ID'],'node':os.uname().nodename},indent=2))
else: raise ValueError(phase)
