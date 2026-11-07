"""Accepted-state scientific figures; no interpolation of missing runs."""
import csv
import pathlib
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

root=pathlib.Path.cwd()
ceq=.1*.05844/(.9*.01801528+.1*.05844)
out=root/'deliverables/figures'
out.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,
                     'svg.fonttype':'none','savefig.dpi':220})
def read(path):
    with path.open() as f:
        rows=list(csv.DictReader(f))
    return {key:np.array([float(r[key]) for r in rows]) for key in rows[0]} if rows else {}
def save(fig,name):
    fig.savefig(out/(name+'.png'),bbox_inches='tight')
    fig.savefig(out/(name+'.svg'),bbox_inches='tight')
    plt.close(fig)
def columns(path,name,times,unit,title):
    data=read(path)
    if not data: return
    available=np.unique(data['time_s'])
    chosen=[available[np.argmin(abs(available-t))] for t in times]
    if len(set(chosen))<4: return
    fig,axes=plt.subplots(1,4,figsize=(10,4.8),sharey=True,layout='constrained')
    max_height=max(data['z_top_m'])*unit
    low=min(.43,float(data['epsilon'].min()))
    for ax,t in zip(axes,chosen):
        ix=data['time_s']==t
        bottoms=data['z_bottom_m'][ix]*unit; tops=data['z_top_m'][ix]*unit
        faces=np.r_[bottoms[0],tops]
        im=ax.pcolormesh([0,1],faces,data['epsilon'][ix,None],vmin=low,vmax=.43,cmap='viridis',shading='flat')
        ax.axhline(tops[-1],color='black',lw=.7)
        ax.set(xlim=(0,1),ylim=(0,max_height*1.04),xticks=[],title=f't = {t:g} s',xlabel='1D column')
    axes[0].set_ylabel('Height above soil ('+('mm' if unit==1e3 else r'$\mu$m')+')')
    fig.colorbar(im,ax=axes,label='Crust porosity',shrink=.8)
    fig.suptitle(title+'\nDisplay width is arbitrary; colors show vertical cell averages',fontsize=11)
    save(fig,name)

columns(root/'stage_a/n640/fields.csv','m3_porosity_columns',[300,3600,10800,21600],1e3,
        'Standalone M3 regression: prescribed saturated reservoir, initial 1 mm crust')

history_path=root/'stage_c/nominal/coupled_history.csv'
if history_path.exists():
    d=read(history_path)
    fields=root/'stage_c/nominal/crust_fields.csv'
    if fields.exists():
        f=read(fields)
        if f:
            t=np.unique(f['time_s'])
            columns(fields,'coupled_porosity_columns',np.linspace(t[0],t[-1],4),1e6,
                    'Conditional saturated-crust coupling: accepted states only')
    old=read(root.parent/'p2_g_surface_precipitation_20261107/runs/n80/surface_precipitation.csv')
    fig,axes=plt.subplots(2,2,figsize=(10,7),layout='constrained')
    ax=axes[0,0]
    ax.plot(d['time_s'],d['h_m']*1e6,label='ALE domain height')
    ax.plot(d['time_s'],d['h_mass_m']*1e6,label='ALE solid-equivalent height',ls=':')
    ax.plot(d['time_s'],d['h_prat_m']*1e6,label='Prat diagnostic, ALE evaporation',ls='--')
    ax.plot(old['time_s'],old['h_mass_m']*1e6,label='P2-G solid-equivalent height',color='black',lw=1)
    ax.plot(old['time_s'],old['h_prat_m']*1e6,label='Original P2-G Prat diagnostic',color='grey',ls='-.',lw=1)
    ax.set_ylabel(r'Height ($\mu$m)'); ax.legend(fontsize=8)
    ax=axes[0,1]
    ax.plot(d['time_s'][1:],d['evaporation_kg_m2_s'][1:]*3.6e6,label='ALE exposed top')
    ax.plot(old['time_s'],old['evaporation_kg_m2_s']*3.6e6,label='P2-G bare-soil film',ls='--')
    ax.set_ylabel(r'Evaporation (g m$^{-2}$ h$^{-1}$)'); ax.legend(fontsize=8)
    common=d['time_s']==1000
    if np.any(common):
        change=(d['evaporation_kg_m2_s'][common][0]/old['evaporation_kg_m2_s'][-1]-1)*100
        ax.text(.48,.60,f'Change at 1000 s: {change:.6f}%',transform=ax.transAxes,fontsize=9,ha='center')
    ax=axes[1,0]
    ax.plot(d['time_s'],(d['top_soil_mass_fraction']-ceq)*1e6,label='With ALE')
    ax.plot(old['time_s'],(old['top_mass_fraction']-ceq)*1e6,label='P2-G',ls='--')
    ax.set_ylabel(r'Soil surface $C-C_{eq}$ (mg/kg brine)'); ax.legend(fontsize=8)
    ax=axes[1,1]
    ax.plot(d['time_s'],d['top_solid_kg_m2']*1e3,label='Top solid, including activation mass')
    ax.plot(d['time_s'],d['internal_solid_kg_m2']*1e3,label='Internal solid')
    ax.plot(old['time_s'],old['solid_kg_m2']*1e3,label='P2-G total solid',ls='--')
    ax.set_ylabel(r'Solid NaCl (g m$^{-2}$)'); ax.legend(fontsize=8)
    for ax in axes.flat:
        ax.set_xlabel('Time (s)'); ax.grid(alpha=.2)
        active=d['time_s'][d['active']>0]
        if active.size: ax.axvline(active[0],color='grey',lw=.8,ls=':')
    fig.suptitle('Conditional 1D feedback; 1-micrometre activation is below the crystal scale\nP2-G comparison ends at its actual 1000-s endpoint',fontsize=11)
    save(fig,'ale_vs_p2g')
    fig,axes=plt.subplots(2,2,figsize=(10,6),layout='constrained')
    series=[('epsilon_min','Minimum crust porosity'),('top_soil_sw','Soil surface saturation'),
            ('salt_residual_kg','Global salt residual (kg)'),('water_residual_kg','Global water residual (kg)')]
    for ax,(key,label) in zip(axes.flat,series):
        values=d[key].copy()
        if key=='epsilon_min': values[d['active']==0]=np.nan
        ax.plot(d['time_s'],values); ax.set(xlabel='Time (s)',ylabel=label); ax.grid(alpha=.2)
    save(fig,'coupled_budgets_porosity')

# Separate soil near-surface profiles at actual saved times.
soil_path=root/'stage_c/nominal/soil_fields.csv'
if soil_path.exists():
    d=read(soil_path)
    fig,ax=plt.subplots(figsize=(5,4.5),layout='constrained')
    times=np.unique(d['time_s'])
    chosen=np.unique([times[np.argmin(abs(times-t))] for t in (0,200,500,1000,2500)])
    for t in chosen:
        ix=(d['time_s']==t) & (d['y_m']>=.195)
        xy=np.unique(np.column_stack([d['y_m'][ix],d['mass_fraction'][ix]]),axis=0)
        ax.plot((xy[:,1]-ceq)*1e6,(.2-xy[:,0])*1000,label=f'{t:g} s',marker='.',ms=3)
    ax.invert_yaxis(); ax.set(xlabel=r'$C-C_{eq}$ (mg/kg brine)',ylabel='Depth below soil surface (mm)')
    ax.legend(); ax.grid(alpha=.2)
    save(fig,'soil_near_surface_profiles')
