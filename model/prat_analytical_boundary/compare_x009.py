"""Matched-concentration comparisons from accepted records only; HPCC compute job."""
import csv,gzip,hashlib,json,os,pathlib,shutil
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=pathlib.Path.cwd(); out=root/'x009_comparison'
assert os.environ.get('SLURM_JOB_ID') and not os.uname().nodename.startswith(('ibcln','login','test'))
(out/'figures').mkdir(parents=True,exist_ok=True); (out/'records').mkdir(exist_ok=True)
prefix=pathlib.Path('/WORK/hulm_work/zhixin')
registered=[
 ('analytical_kr0p0023',root/'runs/x009_analytical_n80_dt0p5',True,.0023),
 ('analytical_kr0p001',root/'runs/x009_analytical_kr0p001_n80_dt0p5',True,.001),
 ('ALE_reversible',prefix/'p2_h_flux_mechanism_x009_20261109/runs/x009_progressive_reservoir_n80_a64_dt0p5',False,.001),
 ('ALE_no_dissolution',prefix/'p2_h_x009_no_dissolution_20261109/runs/x009_progressive_n80_a64_dt0p5_no_dissolution',False,.001),
 ('ALE_reversible_dt_half',prefix/'p2_h_flux_mechanism_x009_20261109/runs/x009_progressive_reservoir_n80_a64_dt0p25_timecheck',False,.001)]
def sha(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()
def csvout(name,rows):
    with (out/name).open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
def load(p): return np.genfromtxt(str(p),delimiter=',',names=True)
ceq=.1*.05844/(.9*.01801528+.1*.05844); rho=998.203+1000*ceq*(.670804+.373854*ceq)
data={}; summary=[]; budgets=[]; hashes={}; specs=[]; trajectories=[]
for name,d,analytical,kr in registered:
    hist=d/('history.csv' if analytical else 'coupled_history.csv'); z=load(hist)
    active=z['active']>0; assert any(active)
    t=z['time_s']; h=z['h_m']; hm=h if analytical else z['h_mass_m']; zeros=np.zeros(len(t))
    m={'time_s':t,'active':z['active'],'h_m':h,'h_mass_m':hm,'solid_kg_m2':hm*2165*.57,
       'positive_internal_kg_m2':zeros if analytical else z['internal_positive_kg_m2'],
       'dissolution_kg_m2':zeros if analytical else z['internal_dissolution_kg_m2'],
       'net_internal_kg_m2':zeros if analytical else z['internal_solid_kg_m2'],
       'epsilon_min':zeros+.43 if analytical else z['epsilon_min'],
       'j_kg_m2_s':z['j_kg_m2_s'] if analytical else z['evaporation_kg_m2_s'],
       'soil_c':z['c_soil'] if analytical else z['top_soil_mass_fraction'],
       'soil_sw':z['sw_soil'] if analytical else z['top_soil_sw'],
       'salt_residual_rel':z['salt_residual_rel'] if analytical else z['salt_residual_kg']/z['soil_salt_kg'][0],
       'water_residual_rel':z['water_residual_rel'] if analytical else z['water_residual_kg']/z['soil_water_kg'][0]}
    activation=float(t[np.where(active)[0][0]]); m['activation_s']=activation
    m['precipitation_s']=float(t[np.where(hm>0)[0][0]])
    data[name]=m
    outcome=json.loads((d/'outcome.json').read_text())
    stop='requested_horizon' if abs(t[-1]-25000)<1e-6 and outcome['exit_code']==0 else 'solver_exception'
    with (d/'solver.out').open('rb') as f:
        f.seek(max(0,(d/'solver.out').stat().st_size-16384)); tail=f.read().decode('utf-8','replace')
    for line in tail.splitlines():
        if 'run failed:' in line: stop=line[line.index('run failed:')+11:].strip()
    summary.append({'case':name,'x_initial':.09,'x_reservoir':.09,'kr_m_s':kr,'dt_max_s':.25 if name.endswith('dt_half') else .5,
      'crust_cells':0 if analytical else 64,'porosity':'fixed 0.43' if analytical else 'evolves from 0.43',
      'internal_reaction':'absent by assumption' if analytical else ('positive-only counterfactual' if 'no_dissolution' in name else 'reversible'),
      'first_precipitation_s':m['precipitation_s'],'activation_s':activation,'last_accepted_s':float(t[-1]),
      'completed_25000':bool(abs(t[-1]-25000)<1e-6 and outcome['exit_code']==0),'exit_code':outcome['exit_code'],'stop_reason':stop,
      'final_h_um':float(h[-1]*1e6),'final_h_mass_um':float(hm[-1]*1e6),'final_solid_g_m2':float(m['solid_kg_m2'][-1]*1000),
      'final_internal_positive_g_m2':float(m['positive_internal_kg_m2'][-1]*1000),'final_dissolution_g_m2':float(m['dissolution_kg_m2'][-1]*1000),
      'final_epsilon_min':float(m['epsilon_min'][-1]),'final_j_kg_m2_s':float(m['j_kg_m2_s'][-1]),'final_soil_c':float(m['soil_c'][-1]),
      'max_salt_residual_rel':float(max(abs(m['salt_residual_rel']))),'max_water_residual_rel':float(max(abs(m['water_residual_rel']))),
      'job':outcome.get('job',outcome.get('slurm_job_id'))})
    if analytical:
        sflux=max(abs((np.diff(z['solid_kg_m2']+z['pore_salt_kg_m2'])-z['salt_out_kg_m2_s'][1:]*z['dt_s'][1:])*.001/z['soil_salt_kg'][0]))
        wflux=max(abs((np.diff(z['pore_water_kg_m2'])-(z['total_out_kg_m2_s'][1:]-z['salt_out_kg_m2_s'][1:]-z['j_kg_m2_s'][1:])*z['dt_s'][1:])*.001/z['soil_water_kg'][0]))
        dyn=active[1:]&active[:-1]; da=2*.57*(1-ceq)*kr*z['h_m'][:-1][dyn]/(.43*np.sqrt(.43)*1.3e-9)
        vref=ceq*z['j_kg_m2_s'][1:][dyn]/(2165*.57*(1-ceq))*da/(1+da)
        vr=max(abs(vref-z['v_m_s'][1:][dyn])/vref)
        events=list(csv.DictReader((d/'events.csv').open()))
        transfer=max(max(abs(float(e['value1'])),abs(float(e['value2']))) for e in events if e['event']=='activation_trial')
        budgets.append({'case':name,'salt_flux_residual_rel':float(sflux),'water_flux_residual_rel':float(wflux),'activation_transfer_residual_rel':transfer,
           'finite_da_velocity_rel_error':float(vr),'max_interface_deviation_active':float(max(abs(m['soil_c'][active]/ceq-1))),
           'min_sw_active':float(min(m['soil_sw'][active])),'max_pe_u':float(max(z['pe_u'])),'max_pe_move':float(max(z['pe_move'])),
           'budget_gate':'PASS' if max(sflux,wflux,transfer)<1e-12 and max(summary[-1]['max_salt_residual_rel'],summary[-1]['max_water_residual_rel'])<1e-8 else 'FAIL'})
    # Explicitly audit effective physical inputs, including command-line overrides.
    config={}; section=''
    for line in (d/'params.input').read_text().splitlines():
        line=line.split('#')[0].strip()
        if line.startswith('['): section=line.strip('[]')
        elif '=' in line:
            k,v=line.split('=',1); config[section+'.'+k.strip()]=v.strip()
    command=json.loads((d/'command.json').read_text())
    for i in range(2,len(command),2): config[command[i].lstrip('-')]=command[i+1]
    keys=['Grid.Positions0','Grid.Positions1','Grid.Cells0','Grid.Cells1','Grid.Grading1','Problem.WaterLevel','Problem.Permeability','Problem.Porosity',
      'Problem.Swr','Problem.PcN','Problem.PcAlpha','Problem.SweL','Problem.SaltMoleFraction','Problem.BottomBoundaryMode','Problem.SaltDiffusivity','Transport.SolubilityMoleFraction']
    specs.append(dict([('case',name)]+[(k,config.get(k,'')) for k in keys]))
    assert float(config['Problem.SaltMoleFraction'])==.09 and float(config['Grid.Cells1'])==80
    target=out/'records'/name; target.mkdir(parents=True,exist_ok=True)
    for filename in ['params.input','command.json','outcome.json','case_spec.json','run_driver.py','events.csv','coupling_events.csv','stop_record.json','last_accepted_soil.csv']:
        if (d/filename).exists(): shutil.copy2(str(d/filename),str(target/filename)); hashes[str(d/filename)]=sha(d/filename)
    for filename in [hist.name,'salt_budget.csv','reservoir_config.csv','crust_fields.csv']:
        if (d/filename).exists(): hashes[str(d/filename)]=sha(d/filename)
    executable=pathlib.Path(command[0]); hashes[str(executable)]=sha(executable)
    hashes[str(d/'solver.out')]=sha(d/'solver.out')
    (target/'solver_tail.txt').write_text(tail)
    if analytical:
        with (d/'solver.out').open('rb') as f,gzip.open(str(target/'solver.out.gz'),'wb') as g: shutil.copyfileobj(f,g)
    selected={0,len(t)-1,int(np.where(active)[0][0])}; next_time=0
    for i in range(len(t)):
        if t[i]>=next_time: selected.add(i); next_time=t[i]+50
    for i in sorted(selected):
        trajectories.append(dict([('case',name)]+[(k,float(v[i])) for k,v in m.items() if isinstance(v,np.ndarray)]))
csvout('case_summary.csv',summary); csvout('boundary_audits.csv',budgets); csvout('effective_inputs.csv',specs); csvout('sampled_histories.csv',trajectories)
prefix_checks=[]
for name in ['analytical_kr0p0023','analytical_kr0p001','ALE_no_dissolution']:
    m=data[name]; ref=data['ALE_reversible']; ix=(m['time_s']>0)&(m['time_s']<=min(m['activation_s'],ref['activation_s']))
    row={'case':name,'reference':'ALE_reversible'}
    for key in ['h_mass_m','soil_c','soil_sw','j_kg_m2_s']:
        row[key+'_max_absolute_difference']=float(max(abs(m[key][ix]-np.interp(m['time_s'][ix],ref['time_s'],ref[key]))))
    row['status']='PASS' if max(v for k,v in row.items() if k.endswith('_difference'))<1e-12 else 'FAIL'
    prefix_checks.append(row)
csvout('preactivation_regression.csv',prefix_checks)
assert all(all(row[k].split()==specs[0][k].split() for k in specs[0] if k!='case') for row in specs)
# Only interpolate between accepted records; never extend any failed run.
primary=list(data)[:4]; common_end=min(data[n]['time_s'][-1] for n in primary); common_start=max(data[n]['activation_s'] for n in primary)
targets=[common_start+v for v in [0,10,100,300,500] if common_start+v<=common_end]+[float(common_end)]
matched=[]
for target in sorted(set(targets)):
    for name,m in data.items():
        if not m['time_s'][0]<=target<=m['time_s'][-1]: continue
        row={'case':name,'physical_time_s':target}
        for k,v in m.items():
            if isinstance(v,np.ndarray) and k not in ['time_s','active']: row[k]=float(np.interp(target,m['time_s'],v))
        matched.append(row)
csvout('matched_physical_times.csv',matched)
labels={'analytical_kr0p0023':'Analytical, kr=2.3e-3','analytical_kr0p001':'Analytical, kr=1e-3',
 'ALE_reversible':'ALE reversible, kr=1e-3','ALE_no_dissolution':'ALE no dissolution, kr=1e-3','ALE_reversible_dt_half':'ALE reversible, dt/2'}
colors=['#0072B2','#D55E00','#009E73','#CC79A7']; styles=['-','--','-.',':']
plt.rcParams.update({'font.size':9,'pdf.fonttype':42})
def save(fig,name):
    fig.tight_layout(); fig.savefig(str(out/'figures'/(name+'.png')),dpi=200); fig.savefig(str(out/'figures'/(name+'.pdf'))); plt.close(fig)
fig,axs=plt.subplots(3,2,figsize=(10,10))
j_equilibrium=float(np.interp(common_start,data['ALE_reversible']['time_s'],data['ALE_reversible']['j_kg_m2_s']))
panels=[('h_m',1e6,'Geometric height (um)'),('h_mass_m',1e6,'Solid-equivalent height (um)'),('solid_kg_m2',1000,'Solid NaCl (g/m2)'),
        ('j_kg_m2_s',1000,'Evaporation (g/m2/s)'),('soil_c',1,'Soil interface mass fraction'),('epsilon_min',1,'Minimum crust porosity')]
for ax,(key,fac,ylabel) in zip(axs.flat,panels):
    for i,name in enumerate(primary):
        m=data[name]; ix=(m['time_s']>=common_start)&(m['time_s']<=common_end)
        values=(m[key][ix]/j_equilibrium-1)*1e6 if key=='j_kg_m2_s' else m[key][ix]*fac
        ax.plot(m['time_s'][ix]-common_start,values,color=colors[i],ls=styles[i],label=labels[name])
    if key=='j_kg_m2_s': ylabel='Evaporation change vs equilibrium film (ppm)'
    ax.set(xlabel='Time since common activation reference (s)',ylabel=ylabel); ax.grid(alpha=.2)
axs[0,0].legend(fontsize=8); axs[2,0].axhline(.98*ceq,color='black',ls=':',lw=.7)
fig.suptitle('x(initial)=x(reservoir)=0.09; common accepted window only\nFixed analytical porosity is an assumption, not an ALE field',fontsize=10)
save(fig,'comparison_common_window')
fig,axs=plt.subplots(1,2,figsize=(11,4.4))
for i,name in enumerate(primary):
    m=data[name]; ix=m['active']>0
    axs[0].plot(m['time_s'][ix],m['h_m'][ix]*1e6,color=colors[i],ls=styles[i],label=labels[name])
    axs[0].plot(m['time_s'][-1],m['h_m'][-1]*1e6,'o',color=colors[i],ms=4)
    if not name.startswith('analytical'):
        axs[1].plot(m['time_s'][ix],m['positive_internal_kg_m2'][ix]*1000,color=colors[i],ls='-',label=name+' precipitation')
        axs[1].plot(m['time_s'][ix],m['dissolution_kg_m2'][ix]*1000,color=colors[i],ls='--',label=name+' dissolution')
axs[0].set(xlabel='Actual physical time (s)',ylabel='Geometric height (um)',title='Endpoints occur at different times'); axs[0].legend(fontsize=7)
axs[1].set(xlabel='Actual physical time (s)',ylabel='Internal reaction inventory (g/m2)',title='Analytical: internal reaction absent by assumption'); axs[1].legend(fontsize=7)
save(fig,'censored_lifetimes_and_reactions')
fig,axs=plt.subplots(2,1,figsize=(8,6),sharex=True)
for i,name in enumerate(primary):
    m=data[name]
    for ax,key in zip(axs,['salt_residual_rel','water_residual_rel']):
        ax.semilogy(m['time_s'],np.maximum(abs(m[key]),1e-18),color=colors[i],ls=styles[i],label=labels[name]); ax.axhline(1e-8,color='black',ls=':',lw=.7)
axs[0].set(ylabel='Absolute relative salt residual'); axs[0].legend(fontsize=8)
axs[1].set(xlabel='Actual physical time (s)',ylabel='Absolute relative water residual'); fig.suptitle('Exact zero residuals displayed at 1e-18')
save(fig,'comparison_conservation')
pairing={'common_start_s':float(common_start),'common_end_s':float(common_end),'interpolation':'Linear within accepted-time brackets only, never beyond a stopped run',
 'evaporation_plot':'1e6*(j/j_equilibrium-1), ppm; avoids ambiguous automatic axis offset','j_equilibrium_kg_m2_s':j_equilibrium,
 'parameters':'Primary analytical kr=0.0023 differs from legacy kr=0.001; added analytical kr=0.001 is the matched coefficient control.',
 'limitations':'Different constitutive closures and stopping criteria remain. No-dissolution is a positive-only reaction counterfactual, not removal of internal precipitation.',
 'analytical_internal_inventory':'Exactly zero by model construction; no spatial internal fields generated.'}
(out/'comparison_protocol.json').write_text(json.dumps(pairing,indent=2))
for filename in ['compare_x009.py','run_x009.py','x009.slurm','compare_x009.slurm']:
    p=root/'repo/model/prat_analytical_boundary'/filename
    if p.exists(): hashes[str(p)]=sha(p)
for pattern in ['x009-*.err','x009-*.out']:
    for p in (root/'logs').glob(pattern):
        if os.environ['SLURM_JOB_ID'] in p.name: continue
        shutil.copy2(str(p),str(out/'records'/p.name)); hashes[str(p)]=sha(p)
(out/'provenance.json').write_text(json.dumps({'job':os.environ['SLURM_JOB_ID'],'node':os.uname().nodename,'sha256':hashes,'cases':[(n,str(d)) for n,d,a,k in registered]},indent=2))
(out/'checksums.json').write_text(json.dumps({str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file() and p.name!='checksums.json'},indent=2))
print(json.dumps({'summary':summary,'budgets':budgets,'protocol':pairing},indent=2))
