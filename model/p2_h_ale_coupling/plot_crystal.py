"""Plots use accepted states only; absent fields remain explicit."""
import csv,json,pathlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
root=pathlib.Path.cwd();out=root/'audit/figures';out.mkdir(exist_ok=True)
d=json.loads((root/'audit/verification.json').read_text());by={c['name']:c for c in d['cases']}
ns=(40,80,160,320);colors=dict(zip(ns,plt.get_cmap('tab10').colors[:4]))
plt.rcParams.update({'font.size':11,'svg.fonttype':'none','savefig.dpi':180})
def read(n,f='coupled_history.csv'):
    with (root/'matrix'/f'crystal_n{n}'/f).open() as h:r=list(csv.DictReader(h))
    return {k:np.array([float(x[k]) for x in r]) for k in r[0]} if r else {}
def save(fig,name):
    for ext in ('png','svg'):fig.savefig(out/(name+'.'+ext),bbox_inches='tight')
    plt.close(fig)
fig,ax=plt.subplots(figsize=(7,4),layout='constrained')
for key,label in [('precipitation_s','First precipitation'),('first_attempt_s','First attempt'),('activation_s','Accepted activation')]:
    ax.plot(ns,[by[f'crystal_n{n}'][key] if by[f'crystal_n{n}'][key] is not None else np.nan for n in ns],'o-',label=label)
ax.set(xlabel='Soil vertical cells',ylabel='Physical time (s)',xticks=ns,title='3.5 micrometre default; 5 mm donor limit');ax.legend();ax.grid(alpha=.2)
save(fig,'figure1_events')
fig,axs=plt.subplots(1,3,figsize=(13,4),layout='constrained')
for n in ns:
    r=read(n)
    for ax,key,factor in zip(axs,('h_m','h_mass_m','internal_solid_kg_m2'),(1e6,1e6,1000)):
        ax.plot(r['time_s'],r[key]*factor,c=colors[n],label=f'{n} cells');ax.grid(alpha=.2);ax.set_xlabel('Time (s)')
for ax,label in zip(axs,('ALE height (micrometres)','Solid-equivalent height (micrometres)','Internal solid (g/m2)')):ax.set_ylabel(label)
axs[0].legend();fig.suptitle('Accepted histories; zero ALE height before actual activation')
save(fig,'figure2_evolution')
fields={n:read(n,'crust_fields.csv') for n in ns};active={n:p for n,p in fields.items() if p}
common=set.intersection(*(set(p['time_s']) for p in active.values())) if active else set();t=max(common) if common else None
if t is None and any(2500 in set(p['time_s']) for p in active.values()):t=2500
fig,axs=plt.subplots(1,2,figsize=(10,5),layout='constrained')
for n in ns:
    p=fields[n]
    if not p or t not in set(p['time_s']):axs[0].plot([],[],c=colors[n],label=f'{n}: no field at {t} s');continue
    m=p['time_s']==t;z=(p['z_bottom_m'][m]+p['z_top_m'][m])/2;h=max(p['z_top_m'][m])
    axs[0].plot(p['epsilon'][m],z*1e6,c=colors[n],label=f'{n} cells');axs[1].plot(p['epsilon'][m],z/h,c=colors[n])
for ax in axs:ax.set_xlabel('ALE porosity');ax.grid(alpha=.2)
axs[0].set_ylabel('Height above soil (micrometres)');axs[1].set_ylabel('z/h');axs[0].legend();fig.suptitle(f'Common stored active time: {t} s; missing grids do not pass')
save(fig,'figure3_porosity')
fig,axs=plt.subplots(1,2,figsize=(11,4),layout='constrained')
for n in ns:
    r=read(n)
    for ax,k,tot in zip(axs,('salt_residual_kg','water_residual_kg'),('soil_salt_kg','soil_water_kg')):ax.plot(r['time_s'],r[k]/r[tot][0],c=colors[n],label=f'{n} cells')
for ax,label in zip(axs,('Salt relative residual','Water relative residual')):
    ax.set(xlabel='Time (s)',ylabel=label,ylim=(-1.5e-8,1.5e-8));ax.set_yscale('symlog',linthresh=1e-14);ax.set_yticks([-1e-8,-1e-10,-1e-12,-1e-14,0,1e-14,1e-12,1e-10,1e-8]);ax.grid(alpha=.2)
    ax.axhline(1e-8,c='k',ls=':');ax.axhline(-1e-8,c='k',ls=':')
axs[0].legend();fig.suptitle('Combined soil + crust + evaporated water budgets')
save(fig,'figure4_conservation')
fig,axs=plt.subplots(1,4,figsize=(13,4),sharey=True,layout='constrained')
for ax,n in zip(axs,ns):
    c=by[f'crystal_n{n}'];p=read(n,'donor_profiles.csv');a=c.get('allocation')
    if not a:ax.text(.1,.5,'No allocation',transform=ax.transAxes);continue
    m=p['attempt']==a['attempt'];z=p['depth_m'][m]*1000;order=np.argsort(z)
    ax.plot(p['sw_before'][m][order],z[order],'.-',label='Before');ax.plot(p['sw_after'][m][order],z[order],'.-',label='After')
    ax.axhspan(0,a['used_depth_m']*1000,color='C2',alpha=.12);ax.axvline(.99,c='k',ls=':');ax.set(title=f'{n} cells',xlabel='Saturation',xlim=(.989,1.001),ylim=(5,0));ax.grid(alpha=.2)
axs[0].set_ylabel('Depth below surface (mm)');axs[0].legend();fig.suptitle('Actual activation transfer; shaded: contributing BOX extent')
save(fig,'figure5_donors')
# Four physical-height columns, preferably the finest activated grid.
if active:
    eligible=[n for n,p in active.items() if len(set(p['time_s']))>=4]
    n=max(eligible) if eligible else max(active);p=active[n];ts=np.array(sorted(set(p['time_s'])));sel=ts[np.linspace(0,len(ts)-1,4).round().astype(int)]
    fig,axs=plt.subplots(1,4,figsize=(10,5),sharey=True,layout='constrained');vmin=min(p['epsilon']);vmax=.43
    for ax,t in zip(axs,sel):
        m=p['time_s']==t;polys=[[(0,a*1e6),(1,a*1e6),(1,b*1e6),(0,b*1e6)] for a,b in zip(p['z_bottom_m'][m],p['z_top_m'][m])]
        pc=PolyCollection(polys,array=p['epsilon'][m],cmap='viridis',clim=(vmin,vmax),edgecolors='none');ax.add_collection(pc);ax.set(xlim=(-.1,1.1),ylim=(0,max(p['z_top_m'])*1.04e6),xticks=[],title=f'{t:g} s')
    axs[0].set_ylabel('Height above soil (micrometres)');fig.colorbar(pc,ax=axs,label='Porosity');fig.suptitle(f'{n} cells; physical-height columns, display width arbitrary')
    save(fig,'figure6_columns')
profiles={n:read(n,'soil_fields.csv') for n in ns};times=set.intersection(*(set(p['time_s']) for p in profiles.values()));t=max(times)
if any(2500 in set(p['time_s']) for p in profiles.values()):t=2500
fig,axs=plt.subplots(1,2,figsize=(10,5),sharey=True,layout='constrained');ceq=.1*.05844/(.1*.05844+.9*.01801528)
for n,p in profiles.items():
    if t not in set(p['time_s']):axs[0].plot([],[],c=colors[n],label=f'{n}: no field at {t} s');continue
    m=p['time_s']==t;z=(.2-p['y_m'][m])*1000
    axs[0].plot(p['saturation'][m],z,'.-',ms=3,c=colors[n],label=f'{n} cells');axs[1].plot((p['mass_fraction'][m]-ceq)*1e6,z,'.-',ms=3,c=colors[n])
for ax in axs:ax.set_ylim(10,0);ax.grid(alpha=.2)
axs[0].set(xlabel='Saturation',ylabel='Depth below soil (mm)');axs[0].ticklabel_format(useOffset=False);axs[0].locator_params(axis='x',nbins=3);axs[1].set_xlabel('C - C_eq (mg/kg brine)');axs[0].legend();fig.suptitle(f'Available soil fields at {t:g} s; missing grids do not pass')
save(fig,'figure7_soil')
