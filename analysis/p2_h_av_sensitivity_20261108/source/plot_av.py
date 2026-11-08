"""Required P2-H av sensitivity figures from accepted states only."""
import csv,pathlib
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle

root=pathlib.Path.cwd(); analysis=root/'analysis'; results=analysis/'results'; figures=analysis/'figures'
figures.mkdir(parents=True,exist_ok=True)
colors={.1:'#377eb8',1.:'#4daf4a',10.:'#e41a1c'}
labels={.1:r'$0.1a_v$',1.:r'$a_v$',10.:r'$10a_v$'}
def read(path):
    with path.open(newline='') as handle: return [{k:(v if k=='status' else float(v)) for k,v in row.items()} for row in csv.DictReader(handle)]
def save(fig,name):
    fig.savefig(figures/(name+'.png'),dpi=220,bbox_inches='tight')
    fig.savefig(figures/(name+'.svg'),bbox_inches='tight')
    plt.close(fig)
def history(av,n=80,suffix=''):
    return read(results/'time_series'/f'n{n}_av{av:g}{suffix}.csv')
def metrics(av,n=80,suffix=''):
    return read(results/'time_series'/f'n{n}_av{av:g}{suffix}_reaction_metrics.csv')
def profile(av,n=80): return read(results/'reaction_profiles'/f'n{n}_av{av:g}.csv')

fig,ax=plt.subplots(figsize=(6.4,4.2))
for av in (.1,1.,10.):
    d=history(av); ax.plot([r['time_s'] for r in d],[r['h_m']*1e6 for r in d],color=colors[av],label=labels[av])
    ax.axvline(next(r['activation_s'] for r in d if r['active']),color=colors[av],alpha=.15,lw=.8)
ax.set(xlabel='Physical time (s)',ylabel='Crust height (μm)'); ax.legend(); save(fig,'figure1_height')

fig,ax=plt.subplots(figsize=(6.4,4.2))
for av in (.1,1.,10.):
    d=[r for r in history(av) if r['active']]; ax.plot([r['time_s'] for r in d],[r['top_growth_velocity_m_s']*1e9 for r in d],color=colors[av],label=labels[av])
ax.set(xlabel='Physical time (s)',ylabel='Top growth velocity (nm s$^{-1}$)'); ax.legend(); save(fig,'figure2_top_velocity')

fig,axs=plt.subplots(1,3,figsize=(12,3.7),sharex=True)
for av in (.1,1.,10.):
    d=[r for r in history(av) if r['active']]; t=[r['time_s'] for r in d]
    for ax,key in zip(axs,('internal_positive_kg_m2','internal_solid_kg_m2','top_new_solid_kg_m2')):
        ax.plot(t,[r[key]*1e3 for r in d],color=colors[av],label=labels[av])
for ax,title in zip(axs,('Cumulative positive internal','Net internal solid','New top solid')):
    ax.set(xlabel='Physical time (s)',ylabel='Mass (g m$^{-2}$)',title=title)
axs[0].legend(); save(fig,'figure3_solid_partitions')

fig,ax=plt.subplots(figsize=(6.4,4.2))
for av in (.1,1.,10.):
    d=[r for r in history(av) if r['active']]; ax.plot([r['time_s'] for r in d],[r['internal_fraction_net'] for r in d],color=colors[av],label=labels[av])
ax.set(xlabel='Physical time (s)',ylabel='Net internal solid fraction',ylim=(0,1)); ax.legend(); save(fig,'figure4_internal_fraction')

profiles={av:profile(av) for av in (.1,1.,10.)}
activation=min(next(r['activation_s'] for r in history(av) if r['active']) for av in (.1,1.,10.))
common_times=sorted(set.intersection(*({r['time_s'] for r in profiles[av]} for av in (.1,1.,10.))))
post=[time for time in common_times if time>activation]
targets=[post[0],post[len(post)//2],post[-1]]
fig,axs=plt.subplots(1,3,figsize=(12,3.7),sharey=True)
for ax,target in zip(axs,targets):
    for av in (.1,1.,10.):
        times=sorted({r['time_s'] for r in profiles[av]}); time=min(times,key=lambda x:abs(x-target))
        d=[r for r in profiles[av] if r['time_s']==time]
        ax.plot([r['reaction_kg_m3_s'] for r in d],[.5*(r['z_bottom_m']+r['z_top_m'])*1e6 for r in d],color=colors[av],label=labels[av])
    ax.axvline(0,color='k',lw=.6); ax.set(xlabel='Signed R (kg m$^{-3}$ s$^{-1}$)',title=f't = {target:.1f} s')
axs[0].set_ylabel('Height above soil (μm)'); axs[0].legend(); save(fig,'figure5_reaction_profiles')

fig,axs=plt.subplots(1,2,figsize=(9,3.8))
for av in (.1,1.,10.):
    d=metrics(av); t=[r['time_s'] for r in d]
    axs[0].plot(t,[r['l95_bottom_m']*1e6 for r in d],color=colors[av],label=labels[av])
    axs[0].plot(t,[r['l95_shortest_m']*1e6 for r in d],color=colors[av],ls='--')
    axs[1].plot(t,[r['l95_bottom_over_h'] for r in d],color=colors[av])
    axs[1].plot(t,[r['l95_shortest_over_h'] for r in d],color=colors[av],ls='--')
axs[0].set(xlabel='Physical time (s)',ylabel='Reaction-zone thickness (μm)',title='Solid: bottom 95%; dashed: shortest 95%')
axs[1].set(xlabel='Physical time (s)',ylabel='Thickness / h',title='Normalized reaction zone'); axs[0].legend(); save(fig,'figure6_reaction_zone')

fig,ax=plt.subplots(figsize=(6.4,4.2))
for av in (.1,1.,10.):
    d=[r for r in history(av) if r['active']]; ax.plot([r['time_s'] for r in d],[r['epsilon_min'] for r in d],color=colors[av],label=labels[av])
for threshold in (.2,.1,.02): ax.axhline(threshold,color='k',ls=':',lw=.7)
ax.set(xlabel='Physical time (s)',ylabel='Minimum porosity',ylim=(0,.45)); ax.legend(); save(fig,'figure7_minimum_porosity')

fig,axs=plt.subplots(1,3,figsize=(10.5,3.8),sharey=True)
for ax,target in zip(axs,targets):
    for av in (.1,1.,10.):
        times=sorted({r['time_s'] for r in profiles[av]}); time=min(times,key=lambda x:abs(x-target))
        d=[r for r in profiles[av] if r['time_s']==time]
        ax.plot([r['epsilon'] for r in d],[.5*(r['z_bottom_m']+r['z_top_m'])*1e6 for r in d],color=colors[av],label=labels[av])
    ax.set(xlabel='Porosity',title=f't = {target:.1f} s')
axs[0].set_ylabel('Height above soil (μm)'); axs[0].legend(); save(fig,'figure8_porosity_profiles')

fig,axs=plt.subplots(1,3,figsize=(10,4),sharey=True)
norm=plt.Normalize(.02,.43); cmap=plt.get_cmap('viridis')
for ax,av in zip(axs,(.1,1.,10.)):
    patches=[]; values=[]; times=sorted({r['time_s'] for r in profiles[av]})
    # The first three columns are directly comparable early states.  The last
    # is each trajectory's final accepted state (2500 s for A1/A2 and the
    # physical-bound stop for A3), with its own time label and actual height.
    column_targets=[activation,680.,695.,max(times)]
    for x,target in enumerate(column_targets):
        time=min(times,key=lambda value:abs(value-target)); d=[r for r in profiles[av] if r['time_s']==time]
        for r in d:
            patches.append(Rectangle((x-.32,r['z_bottom_m']*1e6),.64,(r['z_top_m']-r['z_bottom_m'])*1e6)); values.append(r['epsilon'])
    collection=PatchCollection(patches,cmap=cmap,norm=norm,edgecolor='none'); collection.set_array(np.array(values)); ax.add_collection(collection)
    ax.set_xlim(-.6,3.6); ax.set_xticks(range(4),[f'{t:.0f}' for t in column_targets],rotation=30); ax.autoscale_view(scalex=False,scaley=True)
    ax.set(xlabel='Physical time (s)',title=labels[av])
axs[0].set_ylabel('Actual crust height (μm)'); fig.colorbar(collection,ax=axs,label='Porosity',shrink=.85); save(fig,'figure9_physical_columns')

fig,axs=plt.subplots(1,2,figsize=(9,3.8),sharex=True)
for av in (.1,1.,10.):
    d=history(av); t=[r['time_s'] for r in d]
    axs[0].plot(t,[abs(r['salt_residual_kg']) for r in d],color=colors[av],label=labels[av])
    axs[1].plot(t,[abs(r['water_residual_kg']) for r in d],color=colors[av])
for ax,title in zip(axs,('NaCl','Water')): ax.set_yscale('symlog',linthresh=1e-18); ax.set(xlabel='Physical time (s)',ylabel='Absolute residual (kg)',title=title)
axs[0].legend(); save(fig,'figure10_conservation')
