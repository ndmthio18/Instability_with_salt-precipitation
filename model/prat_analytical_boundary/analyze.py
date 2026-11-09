"""Independent accepted-step audits and plots; raw histories retained on HPCC."""
import csv,gzip,hashlib,json,os,pathlib,shutil,sys
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
root=pathlib.Path.cwd(); out=root/'deliverable'
if not os.environ.get('SLURM_JOB_ID'): raise RuntimeError('compute allocation required')
out.mkdir(exist_ok=True); (out/'figures').mkdir(exist_ok=True); (out/'logs').mkdir(exist_ok=True)
def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()
def table(path,rows):
    if not rows: return
    with path.open('w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)
def load(p): return np.genfromtxt(str(p),delimiter=',',names=True)
c=.1*.05844/(.9*.01801528+.1*.05844); A=.001
rho=998.203+1000*c*(.670804+.373854*c); eps=.43; deff=np.sqrt(eps)*1.3e-9
def gasfilm(w):
    T=293.15; P=1e5; mw=.01801528; ma=.02896
    n=[.11670521452767e4,-.72421316703206e6,-.17073846940092e2,.12020824702470e5,-.32325550322333e7,.14915108613530e2,-.48232657361591e4,.40511340542057e6,-.23855557567849,.65017534844798e3]
    sig=T+n[8]/(T-n[9]); aa=(sig+n[0])*sig+n[1]; bb=(n[2]*sig+n[3])*sig+n[4]; cc=(n[5]*sig+n[6])*sig+n[7]
    ps=1e6*(2*cc/(np.sqrt(bb*bb-4*aa*cc)-bb))**4
    m=w/(.05844*(1-w)); r=np.sqrt(m); osm=1-.3915*r/(1+1.2*r)+m*(.0765+.2664*np.exp(-2*r))+.00127*m*m
    activity=np.exp(-2*m*mw*osm); ya=.4*ps/P; ys=activity*ps/P
    gasmass=lambda y:mw*y/(mw*y+ma*(1-y))
    rhog=P*(ya*mw+(1-ya)*ma)/(8.314472*T)
    return rhog*(2.13e-5*(T/273.15)**1.8)/.0083*(gasmass(ys)-gasmass(ya))
summaries=[]; audits=[]; hashes={}; series={}; source_data=[]; activation_table=[]
for d in sorted((root/'runs').iterdir()):
    if not (d/'outcome.json').exists(): continue
    r=json.loads((d/'outcome.json').read_text()); z=load(d/'history.csv'); series[d.name]=z
    last=r['last']; active=z['active']>0
    salt=(np.diff(z['solid_kg_m2']+z['pore_salt_kg_m2'])-z['salt_out_kg_m2_s'][1:]*z['dt_s'][1:])*A/z['soil_salt_kg'][0]
    water=(np.diff(z['pore_water_kg_m2'])-(z['total_out_kg_m2_s'][1:]-z['salt_out_kg_m2_s'][1:]-z['j_kg_m2_s'][1:])*z['dt_s'][1:])*A/z['soil_water_kg'][0]
    e=list(csv.DictReader((d/'events.csv').open()))
    activations=[x for x in e if x['event']=='accepted_activation']
    trials=[x for x in e if x['event']=='activation_trial']
    fp=[float(x['time_s']) for x in e if x['event']=='accepted_first_precipitation']
    az=z[active][0] if any(active) else None
    ai=int(np.where(active)[0][0]) if any(active) else None
    activation_table.append({'case':d.name,'first_precipitation_s':fp[0] if fp else None,
       'accepted_activation_s':float(az['time_s']) if az is not None else None,
       'actual_activation_h_um':float(az['h_m'])*1e6 if az is not None else None,
       'activation_solid_kg_m2':float(az['solid_kg_m2']) if az is not None else None,
       'activation_pore_salt_kg_m2':float(az['pore_salt_kg_m2']) if az is not None else None,
       'activation_pore_water_kg_m2':float(az['pore_water_kg_m2']) if az is not None else None,
       'last_precursor_salt_flux_kg_m2_s':float(az['salt_out_kg_m2_s']) if az is not None else None,
       'first_growth_salt_flux_kg_m2_s':float(z[ai+1]['salt_out_kg_m2_s']) if ai is not None and ai+1<len(z) else None,
       'activation_flux_ratio':float(z[ai+1]['salt_out_kg_m2_s']/az['salt_out_kg_m2_s']) if ai is not None and ai+1<len(z) else None,
       'final_growth_evap_relative_to_eq_film':float(last['j_kg_m2_s'])/float(gasfilm(c))-1 if any(active) else None})
    max_activation=max([max(abs(float(x['value1'])),abs(float(x['value2']))) for x in trials] or [0])
    monotonic=bool(np.all(np.diff(z['time_s'])>0)); intervals=bool(np.max(abs(np.diff(z['time_s'])-z['dt_s'][1:]))<1e-9)
    ar={'case':d.name,'salt_flux_residual_rel':float(max(abs(salt))), 'water_flux_residual_rel':float(max(abs(water))),
        'global_salt_residual_rel':r['max_salt_residual'],'global_water_residual_rel':r['max_water_residual'],
        'activation_transfer_residual_rel':max_activation,'accepted_activation_count':len(activations),
        'strictly_increasing_time':monotonic,'dt_matches_accepted_interval':intervals,
        'max_interface_relative_error_active':float(max(abs(z['c_soil'][active]/c-1))) if any(active) else None,
        'max_pe_u':float(max(z['pe_u'])),'max_pe_move':float(max(z['pe_move']))}
    ar['zero_initial_inventories']=all(z[k][0]==0 for k in ['h_m','solid_kg_m2','pore_salt_kg_m2','pore_water_kg_m2','active'])
    ar['height_solid_identity_max_absolute_m']=float(max(abs(z['h_m']-z['solid_kg_m2']/(2165*(1-eps)))))
    ar['numerical_budget_gate']='PASS' if max(ar[k] for k in ['salt_flux_residual_rel','water_flux_residual_rel','activation_transfer_residual_rel'])<=1e-12 and max(r['max_salt_residual'],r['max_water_residual'])<=1e-8 and monotonic and intervals and ar['zero_initial_inventories'] else 'FAIL'
    dyn=(z['active'][1:]>0)&(z['active'][:-1]>0)
    da=2*(1-eps)*(1-c)*.0023*z['h_m'][:-1][dyn]/(eps*deff)
    vref=c*z['j_kg_m2_s'][1:][dyn]/(2165*(1-eps)*(1-c))*da/(1+da)
    vobs=z['v_m_s'][1:][dyn]
    ar['finite_da_velocity_relative_error']=float(max(abs(vref-vobs)/np.maximum(vref,1e-300))) if len(vref) else None
    j=z['j_kg_m2_s'][1:][dyn]; ct=z['c_top'][1:][dyn]
    ar['gasfilm_coupling_relative_error']=float(max(abs(gasfilm(ct)-j)/j)) if len(j) else None
    ar['top_kinetic_absolute_error']=float(max(abs(ct-c-2165*vobs/(2*rho*.0023)))) if len(j) else None
    audits.append(ar)
    summaries.append({'case':d.name,'requested_s':r['requested_s'],'accepted_s':float(last['time_s']),'completed':r['completed'],
       'activation_s':float(last['activation_s']),'h_final_um':float(last['h_m'])*1e6,'solid_final_g_m2':float(last['solid_kg_m2'])*1000,
       'evaporation_final_kg_m2_s':float(last['j_kg_m2_s']),'evaporated_kg_m2':float(last['vapour_kg'])/A,
       'salt_residual_rel':r['max_salt_residual'],'water_residual_rel':r['max_water_residual'],
       'coupling_gate':'PASS' if r['coupling_pass'] else 'FAIL','job':r['job'],'node':r['node']})
    selected={0,len(z)-1}; next_t=100
    for i,t in enumerate(z['time_s']):
        if t>=next_t: selected.add(i); next_t=t+100
        if i and z['active'][i]!=z['active'][i-1]: selected.update([i-1,i])
    for i in sorted(selected): source_data.append(dict([('case',d.name)]+[(k,float(z[k][i])) for k in z.dtype.names]))
    target=out/'records'/d.name; target.mkdir(parents=True,exist_ok=True)
    for name in ['outcome.json','command.json','params.input','events.csv','boundary_faces.csv','last_accepted_soil.csv']:
        if (d/name).exists(): shutil.copy2(str(d/name),str(target/name))
    for p in d.iterdir():
        if p.is_file(): hashes[str(p.relative_to(root))]=digest(p)
    with (d/'solver.out').open('rb') as f,gzip.open(str(out/'logs'/(d.name+'.out.gz')),'wb') as g: shutil.copyfileobj(f,g)
table(out/'convergence.csv',summaries); table(out/'diagnostic_mass_budgets.csv',audits); table(out/'sampled_history.csv',source_data)
table(out/'activation_comparison.csv',activation_table)
comparison=[]
names=['nested_n40_dt0p5','pilot_n80_dt0p5','nested_n160_dt0p5','nested_n320_dt0p5']
metrics=['h_m','solid_kg_m2','j_kg_m2_s','vapour_kg']
for a,b,tol,kind in [(names[i],names[i+1],.02,'grid') for i in range(3)]+[(names[1],'nested_n80_dt0p25',.01,'timestep')]:
    row={'kind':kind,'coarse':a,'fine':b,'tolerance':tol,'status':'NOT TESTED'}
    for k in metrics: row[k+'_relative_difference']=None
    if a in series and b in series:
        x,y=series[a][-1],series[b][-1]
        if abs(x['time_s']-25000)<1e-6 and abs(y['time_s']-25000)<1e-6 and x['active'] and y['active']:
            for k in metrics: row[k+'_relative_difference']=float(abs(x[k]-y[k])/max(abs(y[k]),1e-300))
            row['status']='PASS' if max(row[k+'_relative_difference'] for k in metrics)<=tol else 'FAIL'
        else: row['status']='FAIL'
    comparison.append(row)
table(out/'grid_timestep_comparisons.csv',comparison)
# Hydrostatic checks compare actual profiles at t=0 and the accepted endpoint.
hp=load(root/'runs/hydrostatic/soil_profiles.csv'); first=hp[hp['time_s']==0]; end=hp[hp['time_s']==max(hp['time_s'])]
hydro={k:float(max(abs(first[k]-end[k]))) for k in ['c','sw','pressure_pa']}
hydro['status']='PASS' if hydro['c']<1e-12 and hydro['sw']<1e-12 and hydro['pressure_pa']<1e-6 else 'FAIL'
rollback={'status':'NOT TESTED'}
if 'rollback_injection' in series:
    d=root/'runs/rollback_injection'; ev=(d/'events.csv').read_text(); z=series['rollback_injection']; p=series['pilot_n80_dt0p5']
    pa={float(x['time_s']):x for x in p}; before=float(p[-1]['activation_s'])-.5
    common=[(x,pa[float(x['time_s'])]) for x in z if x['time_s']<=before and float(x['time_s']) in pa]
    diff=max(abs(float(x[k]-y[k])) for x,y in common for k in ['soil_salt_kg','soil_water_kg','solid_kg_m2','pore_salt_kg_m2','pore_water_kg_m2'])
    rollback={'injected_rejection_observed':'rejected_injected rollback' in ev,'accepted_prefix_max_difference':diff,
      'last_accepted_s':float(z[-1]['time_s']),'status':'PASS' if 'rejected_injected rollback' in ev and diff<1e-12 and any(z['active']) else 'FAIL'}
unit={p.name:json.loads(p.read_text()) for p in (root/'unit').glob('*.json')}
assessment={'hydrostatic_profile':hydro,'rollback':rollback,'units':unit,'budget_audits':audits,'comparisons':comparison,
    'pilot_gate':next(s['coupling_gate'] for s in summaries if s['case']=='pilot_n80_dt0p5'),
    'quantitative_research':'NO-GO' if any(x['status']!='PASS' for x in comparison) else 'NUMERICAL GATES ONLY; physical closure remains unvalidated'}
(out/'assessment.json').write_text(json.dumps(assessment,indent=2))
production={k:z for k,z in series.items() if k.startswith(('pilot','nested'))}
colors=['#0072B2','#D55E00','#009E73','#CC79A7','#555555']; styles=['-','--','-.',':','-']
plt.rcParams.update({'font.size':10,'pdf.fonttype':42,'axes.spines.right':False,'axes.spines.top':False})
def save(fig,name):
    fig.tight_layout(); fig.savefig(str(out/'figures'/(name+'.png')),dpi=200); fig.savefig(str(out/'figures'/(name+'.pdf'))); plt.close(fig)
for num,key,factor,ylabel in [(1,'h_m',1e6,'Crust height (um)'),(2,'solid_kg_m2',1000,'Solid NaCl (g/m2)'),(3,'j_kg_m2_s',1000,'Evaporation (g/m2/s)')]:
    fig,ax=plt.subplots(figsize=(7,4))
    for i,(name,z) in enumerate(production.items()):
        label=name.replace('nested_','').replace('pilot_','')
        if z[-1]['time_s']<25000-1e-6: label+=' (stopped)'
        curve=z[z['dt_s']>0] if key=='j_kg_m2_s' else z
        ax.plot(curve['time_s'],curve[key]*factor,label=label,color=colors[i],ls=styles[i])
        ax.plot(z['time_s'][-1],z[key][-1]*factor,'o',color=colors[i],ms=4)
    ax.set(xlabel='Accepted physical time (s)',ylabel=ylabel); ax.legend(fontsize=8); ax.grid(alpha=.2)
    if key=='j_kg_m2_s':
        ax.set_ylim(bottom=0); ax.set_title('Fluxes shown only for accepted time intervals',fontsize=10)
    save(fig,'figure{}_{}'.format(num,key))
fig,axs=plt.subplots(2,1,figsize=(7,6),sharex=True)
for i,(name,z) in enumerate(production.items()):
    axs[0].plot(z['time_s'],z['c_soil']/c-1,label=name,color=colors[i],ls=styles[i])
    axs[1].plot(z['time_s'],z['sw_soil'],color=colors[i],ls=styles[i])
axs[0].axhline(.02,color='black',ls=':',lw=.8); axs[0].axhline(-.02,color='black',ls=':',lw=.8)
axs[0].set(ylabel='Soil interface C/Ceq - 1',title='Validity screens apply after activation'); axs[0].legend(fontsize=8)
axs[1].axhline(.99,color='black',ls=':',lw=.8); axs[1].set(xlabel='Accepted physical time (s)',ylabel='Soil interface saturation')
save(fig,'figure4_interface')
fig,axs=plt.subplots(2,1,figsize=(7,6),sharex=True)
for i,(name,z) in enumerate(production.items()):
    for ax,key in zip(axs,['salt_residual_rel','water_residual_rel']):
        ax.semilogy(z['time_s'],np.maximum(abs(z[key]),1e-18),label=name,color=colors[i],ls=styles[i]); ax.axhline(1e-8,color='black',ls=':',lw=.8)
axs[0].set(ylabel='Absolute relative salt residual'); axs[0].legend(fontsize=8)
axs[1].set(xlabel='Accepted physical time (s)',ylabel='Absolute relative water residual')
fig.suptitle('Exact zeros displayed at 1e-18',fontsize=9)
save(fig,'figure5_conservation')
for p in (root/'unit').iterdir():
    if p.is_file(): shutil.copy2(str(p),str(out/'logs'/p.name))
frozen=root/'runs/frozen_short'
if frozen.exists():
    target=out/'records/frozen_short'; target.mkdir(parents=True,exist_ok=True)
    for name in ['command.json','binary.sha256','params.input']:
        shutil.copy2(str(frozen/name),str(target/name))
    for p in frozen.iterdir():
        if p.is_file(): hashes[str(p.relative_to(root))]=digest(p)
    with (frozen/'solver.out').open('rb') as f,gzip.open(str(out/'logs/frozen_short.out.gz'),'wb') as g: shutil.copyfileobj(f,g)
for p in list((root/'logs').glob('*'))+list((root/'build').glob('*.log')):
    if p.is_file():
        if ('analysis-'+os.environ['SLURM_JOB_ID']) in p.name: continue  # still-open log is not immutable evidence
        hashes[str(p.relative_to(root))]=digest(p)
        with p.open('rb') as f,gzip.open(str(out/'logs'/(p.name+'.gz')),'wb') as g: shutil.copyfileobj(f,g)
for directory in ['repo/model/prat_analytical_boundary','repo/model/p1_r2_conservative_richards/source/appl']:
    for p in (root/directory).rglob('*'):
        if p.is_file(): hashes[str(p.relative_to(root))]=digest(p)
for rel in ['build/soil_solver','build/compiled_main.cc','build/unit_test','repo/model/p2_f_analytical_growth/cpp/analyticalgrowth.hh',
 'repo/model/p2_h_ale_coupling/soil_methods.inc','repo/model/p2_h_ale_coupling/grid_diagnostics.inc','repo/model/p2_h_ale_coupling/progressive.hh',
 'repo/model/p2_d_coupling/cpp/pratcrust.hh','repo/model/p2_g_surface_precipitation/hpcc.opts']:
    hashes[rel]=digest(root/rel)
for p in root.glob('prat_source*.tar.gz'): hashes[p.name]=digest(p)
provenance={'base_commit':'fdd4d1f306c09caf3f44cac39f9cc7431aadfaa0','branch':'research/prat-analytical-boundary-finite-da',
 'remote_root':str(root),'postprocess_job':os.environ['SLURM_JOB_ID'],'postprocess_node':os.uname().nodename,
 'sha256':hashes,'history_sampling':'Retained accepted records every >=100 s, both sides of activation, first and last. No interpolation.',
 'figures':'Full accepted histories; stopped cases end at actual accepted time. Evaporation excludes dt=0 initialization placeholder (no physical interval). No crust fields exist.',
 'raw':'Full step histories, VTK, source tarballs, prior build attempts and binaries retained under remote_root; not added to Git.'}
(out/'provenance.json').write_text(json.dumps(provenance,indent=2))
manifest={str(p.relative_to(out)):digest(p) for p in out.rglob('*') if p.is_file() and p.name!='checksums.json'}
(out/'checksums.json').write_text(json.dumps(manifest,indent=2))
print(json.dumps(assessment,indent=2))
