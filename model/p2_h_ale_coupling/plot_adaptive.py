"""Six requested figures from accepted states, with explicit unavailable cases."""
import csv,json,pathlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
root=pathlib.Path.cwd(); out=root/'audit/figures'; out.mkdir(exist_ok=True)
result=json.loads((root/'audit/verification.json').read_text())
by={r['name']:r for r in result['cases'] if r['stage']=='matrix'}
plt.rcParams.update({'font.size':11,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','savefig.dpi':180})
colors=dict(zip((40,80,160,320),plt.get_cmap('tab10').colors[:4]))
def read(name,file='coupled_history.csv',stage='matrix'):
    with (root/stage/name/file).open() as f: rows=list(csv.DictReader(f))
    return {k:np.array([float(r[k]) for r in rows]) for k in rows[0]} if rows else {}
def save(fig,name):
    fig.savefig(out/(name+'.png'),bbox_inches='tight'); fig.savefig(out/(name+'.svg'),bbox_inches='tight'); plt.close(fig)
fig,axs=plt.subplots(2,2,figsize=(11,8),sharex=True,sharey=True,layout='constrained')
for ax,h in zip(axs.flat,(1,2,5,10)):
    for d,style in zip((2,5,10),('-', '--', ':')):
        for n in (40,80,160,320):
            r=by[f'new_n{n}_h{h}_d{d}']; t=r['activation_s']
            ax.scatter(n,t if t is not None else r['accepted_time_s'],marker='o' if t is not None else 'x',color={2:'C0',5:'C1',10:'C2'}[d],s=30)
        times=[by[f'new_n{n}_h{h}_d{d}']['activation_s'] for n in (40,80,160,320)]
        ax.plot((40,80,160,320),[np.nan if t is None else t for t in times],style,color={2:'C0',5:'C1',10:'C2'}[d],label=f'Adaptive {d} mm')
    oldnames=[f'old_n{n}' if h==1 else f'old_n80_h{h}' for n in ((40,80,160,320) if h==1 else (80,))]
    for i,name in enumerate(oldnames):
        r=by[name]; t=r['activation_s']; ax.scatter(r['n'],t if t is not None else r['accepted_time_s'],c='k',marker='s' if t is not None else 'x',label='Top only' if i==0 else None)
    ax.set_title(f'Activation threshold {h} '+r'$\mu$m'); ax.set_xticks([40,80,160,320]); ax.grid(alpha=.2)
axs[0,0].legend(fontsize=9); fig.supxlabel('Soil vertical cells'); fig.supylabel('Accepted activation time (s)')
fig.suptitle('Cross: no activation, shown at accepted endpoint; 1/2 micrometres are subcrystal demonstrations')
save(fig,'figure1_activation')

fig,axs=plt.subplots(1,3,figsize=(13,5),sharey=True,layout='constrained')
for ax,n in zip(axs,(80,160,320)):
    name=f'new_n{n}_h1_d5'; r=by[name]; p=read(name,'donor_profiles.csv')
    if not r.get('allocation') or not p:
        ax.set(title=f'{n} cells',xlim=(.989,1.001),ylim=(5,0));ax.text(.1,.5,'No actual activation\nprofile available',transform=ax.transAxes);continue
    a=r['allocation']['attempt']; mask=p['attempt']==a
    # Equal-depth lateral DOFs overlap in Sw; sum their contributed mass.
    depths=sorted(set(p['depth_m'][mask])); before=[]; after=[]; mass=[]
    for z in depths:
        m=mask & (p['depth_m']==z); before.append(p['sw_before'][m].mean());after.append(p['sw_after'][m].mean());mass.append(p['extracted_kg'][m].sum()*1e6)
    z=np.array(depths)*1000
    ax.plot(before,z,'o-',ms=3,label='Before');ax.plot(after,z,'o-',ms=3,label='After')
    ax.axvline(.99,c='k',ls=':',lw=1);ax.axhspan(0,r['allocation']['used_depth_m']*1000,color='C2',alpha=.12)
    second=ax.twiny();second.plot(mass,z,'C2.-',label='Transferred');second.set_xlabel('Transferred liquid (mg)',color='C2')
    second.set_xlim(-.025,.55)
    ax.set(xlabel='Liquid saturation',title=f'{n} cells',xlim=(.989,1.001),ylim=(5,0));ax.grid(alpha=.2)
axs[0].set_ylabel('Depth below soil surface (mm)');axs[0].legend(loc='lower left')
fig.suptitle('1 micrometre / 5 mm limit; shaded: used BOX extent; allocation trial shown')
save(fig,'figure2_donor_profiles')

fig,axs=plt.subplots(2,2,figsize=(11,8),sharex=True,sharey=True,layout='constrained')
for ax,h in zip(axs.flat,(1,2,5,10)):
    for n in (40,80,160,320):
        r=read(f'new_n{n}_h{h}_d5');mask=r['active']>0
        if mask.any():ax.plot(r['time_s'][mask],r['h_m'][mask]*1e6,c=colors[n],label=f'{n} cells')
        else:ax.plot([],[],c=colors[n],label=f'{n}: no activation')
    ax.set_title(f'{h} '+r'$\mu$m activation');ax.grid(alpha=.2);ax.set_xlim(0,2500)
axs[0,0].legend();fig.supxlabel('Time (s)');fig.supylabel(r'ALE height ($\mu$m)');fig.suptitle('Adaptive donor, 5 mm limit; curves end at actual accepted times')
save(fig,'figure3_height')

fig,axs=plt.subplots(1,2,figsize=(10,5),layout='constrained')
common_profiles={}; timesets=[]
for n in (40,80,160,320):
    p=read(f'new_n{n}_h1_d5','crust_fields.csv');common_profiles[n]=p
    if p:timesets.append(set(p['time_s']))
common=set.intersection(*timesets) if timesets else set()
t=max(common) if common else None
for n,p in common_profiles.items():
    if t is None or not p:
        axs[0].plot([],[],c=colors[n],ls=':',label=f'{n}: no active field');continue
    m=p['time_s']==t;z=(p['z_bottom_m'][m]+p['z_top_m'][m])/2
    axs[0].plot(p['epsilon'][m],z*1e6,c=colors[n],label=f'{n} cells')
    axs[1].plot(p['epsilon'][m],z/p['z_top_m'][m].max(),c=colors[n])
for ax in axs:ax.set_xlabel('Crust porosity');ax.grid(alpha=.2)
axs[0].set_ylabel(r'Height above soil ($\mu$m)');axs[1].set_ylabel('Mapped height z/h');axs[0].legend()
fig.suptitle(f'1 micrometre / 5 mm: common time among available active grids = {t} s')
if len(timesets)<4:axs[1].text(.04,.04,'Missing grids do not pass\nthe four-grid convergence gate',transform=axs[1].transAxes,fontsize=9)
if t is None:axs[0].text(.1,.5,'No common active fields',transform=axs[0].transAxes)
save(fig,'figure4_porosity_profiles')

name='new_n160_h5_d5';p=read(name,'crust_fields.csv')
if not p:name='new_n160_h1_d5';p=read(name,'crust_fields.csv')
fig,axs=plt.subplots(1,4,figsize=(10,5),sharey=True,layout='constrained')
if p:
    times=np.unique(p['time_s']); selected=times[np.round(np.linspace(0,len(times)-1,4)).astype(int)]
    vmax=max(.43,max(p['epsilon']));vmin=min(p['epsilon']);top=max(p['z_top_m'])*1e6
    for ax,t in zip(axs,selected):
        m=p['time_s']==t; polygons=[[(0,b*1e6),(1,b*1e6),(1,e*1e6),(0,e*1e6)] for b,e in zip(p['z_bottom_m'][m],p['z_top_m'][m])]
        pc=PolyCollection(polygons,array=p['epsilon'][m],cmap='viridis',clim=(vmin,vmax),edgecolors='none');ax.add_collection(pc)
        ax.set(xlim=(-.1,1.1),ylim=(0,top*1.04),xticks=[],title=f'{t:g} s');ax.axhline(0,c='k')
    fig.colorbar(pc,ax=axs,label='Porosity');axs[0].set_ylabel(r'Height above soil ($\mu$m)')
    (root/'audit/column_metadata.json').write_text(json.dumps(dict(case=name,times_s=selected.tolist(),vmin=vmin,vmax=vmax),indent=2))
else:axs[0].text(.1,.5,'No active crust fields',transform=axs[0].transAxes)
fig.suptitle(name+'\nPhysical-height columns; arbitrary display width, not a 2D field')
save(fig,'figure5_porosity_columns')

fig,axs=plt.subplots(1,2,figsize=(11,4),layout='constrained')
for n in (40,80,160,320):
    r=read(f'new_n{n}_h1_d5')
    for ax,key,total in zip(axs,('salt_residual_kg','water_residual_kg'),('soil_salt_kg','soil_water_kg')):
        ax.plot(r['time_s'],r[key]/r[total][0],c=colors[n],label=f'{n} cells');ax.grid(alpha=.2);ax.set_xlabel('Time (s)')
axs[0].set_ylabel('Signed relative salt residual');axs[1].set_ylabel('Signed relative water residual');axs[0].legend()
for ax in axs:
    ax.set_yscale('symlog',linthresh=1e-14);ax.set_ylim(-1.5e-8,1.5e-8)
    # Explicit ticks keep the linear neighbourhood of zero readable.
    ax.set_yticks([-1e-8,-1e-10,-1e-12,-1e-14,0,1e-14,1e-12,1e-10,1e-8])
    ax.axhline(1e-8,c='k',ls=':',lw=.8);ax.axhline(-1e-8,c='k',ls=':',lw=.8)
fig.suptitle('Accepted prefixes, 1 micrometre / 5 mm; dotted: budget gate, not proof of completion')
save(fig,'figure6_conservation')

fig,axs=plt.subplots(1,2,figsize=(10,5),sharey=True,layout='constrained')
profiles={n:read(f'new_n{n}_h1_d5','soil_fields.csv') for n in (40,80,160,320)}
times=set.intersection(*(set(p['time_s']) for p in profiles.values()));t=max(times)
ceq=.1*.05844/(.1*.05844+.9*.01801528)
for n,p in profiles.items():
    m=p['time_s']==t;z=(.2-p['y_m'][m])*1000
    axs[0].plot(p['saturation'][m],z,'.-',ms=3,c=colors[n],label=f'{n} cells')
    axs[1].plot((p['mass_fraction'][m]-ceq)*1e6,z,'.-',ms=3,c=colors[n])
for ax in axs:ax.set_ylim(10,0);ax.grid(alpha=.2)
axs[0].set(xlabel='Soil liquid saturation',ylabel='Depth (mm)');axs[1].set_xlabel(r'$C-C_{eq}$ (mg/kg brine)');axs[0].legend()
axs[0].ticklabel_format(axis='x',useOffset=False);axs[0].locator_params(axis='x',nbins=3)
fig.suptitle(f'Same-time near-surface fields at {t:g} s; 1 micrometre / 5 mm')
save(fig,'figure7_soil_profiles')

available={n:p for n,p in profiles.items() if by[f'new_n{n}_h1_d5']['activation_s'] is not None}
if available and len(available)<4:
    t=max(set.intersection(*(set(p['time_s']) for p in available.values())))
    fig,axs=plt.subplots(1,2,figsize=(10,5),sharey=True,layout='constrained')
    for n,p in available.items():
        m=p['time_s']==t;z=(.2-p['y_m'][m])*1000
        axs[0].plot(p['saturation'][m],z,'.-',ms=3,c=colors[n],label=f'{n} cells')
        axs[1].plot((p['mass_fraction'][m]-ceq)*1e6,z,'.-',ms=3,c=colors[n])
    for ax in axs:ax.set_ylim(10,0);ax.grid(alpha=.2)
    axs[0].set(xlabel='Soil liquid saturation',ylabel='Depth (mm)');axs[1].set_xlabel(r'$C-C_{eq}$ (mg/kg brine)');axs[0].legend()
    fig.suptitle(f'Available activated grids at {t:g} s; missing grids have no comparable coupled field')
    save(fig,'figure7b_available_coupled_soil')

fig,axs=plt.subplots(1,2,figsize=(11,4),layout='constrained')
for n in (40,80,160,320):
    p=read(f'pre_n{n}')
    axs[0].plot(p['time_s'],p['top_solid_kg_m2']*1000,c=colors[n],label=f'{n} cells')
    axs[1].plot(p['time_s'],(p['top_soil_mass_fraction']-ceq)*1e6,c=colors[n])
axs[0].set_ylabel(r'External precipitated NaCl (g m$^{-2}$)');axs[1].set_ylabel(r'$C-C_{eq}$ (mg/kg brine)')
for ax in axs:ax.set_xlabel('Time (s)');ax.grid(alpha=.2)
axs[0].legend();fig.suptitle('Precrust-only diagnostic: no ALE activation or pore-liquid transfer')
save(fig,'figure8_precrust')

fig,axs=plt.subplots(2,4,figsize=(14,8),sharey=True,layout='constrained')
controls={r['name']:r for r in result['controlled']}
for col,n in enumerate((40,80,160,320)):
    name=f'control_n{n}_h5_d5';p=read(name,'donor_profiles.csv','controlled');a=controls[name]['allocation']
    z=p['depth_m']*1000
    axs[0,col].plot(p['sw_before'],z,'.-',ms=2,label='Before');axs[0,col].plot(p['sw_after'],z,'.-',ms=2,label='After')
    axs[0,col].set(xlim=(.989,1.001),title=f'{n} cells; {a["donor_count"]:g} DOFs',xlabel='Liquid saturation')
    axs[1,col].plot(p['p_after_pa']-p['p_before_pa'],z,'.-',ms=2,c=colors[n]);axs[1,col].set(xlim=(-1400,30),xlabel='Pressure change (Pa)')
    for row in (0,1):
        axs[row,col].axhspan(0,a['used_depth_m']*1000,color='C2',alpha=.12);axs[row,col].set_ylim(5,0);axs[row,col].grid(alpha=.2)
axs[0,0].legend();axs[0,0].set_ylabel('Depth (mm)');axs[1,0].set_ylabel('Depth (mm)')
fig.suptitle('Controlled allocation ONLY: same prescribed 5-micrometre skeleton, analytic initial soil\n5 mm bound; no transport, no production salt seed; shaded: used BOX extent')
save(fig,'figure9_controlled_profiles')
