"""Figures use actual accepted times; failed runs are never extended."""
import os,pathlib
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import PolyCollection
from matplotlib.colors import Normalize, SymLogNorm
from matplotlib.ticker import FormatStrFormatter, MaxNLocator
assert os.environ.get('SLURM_JOB_ID')
root=pathlib.Path.cwd();out=root/'audit/figures';out.mkdir(parents=True,exist_ok=True)
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':120,'savefig.dpi':240})
colors={40:'#0072B2',80:'#E69F00',160:'#009E73',320:'#CC79A7'}
def read(p):return pd.read_csv(p)
def hist(n,fixed=False,coupled=False):
    return read(root/(f'correction/fixed_{"coupled" if coupled else "pre"}/fixed{n}' if fixed else f'diagnostic/historical{n}')/'coupled_history.csv')
def finish(fig,name):
    fig.savefig(out/(name+'.png'),bbox_inches='tight');fig.savefig(out/(name+'.svg'),bbox_inches='tight');plt.close(fig)

fig,ax=plt.subplots(1,2,figsize=(11,4),sharey=True)
for n in colors:
    for k,fixed in enumerate((False,True)):
        h=hist(n,fixed);ax[k].plot(h.time_s,h.top_solid_kg_m2*1e3,color=colors[n],label=f'{n} layers')
for a,title in zip(ax,['Historical step-end projection','Implicit equilibrium surface flux']):
    a.set(xlabel='Time (s)',title=title,xlim=(0,2500));a.legend(frameon=False)
ax[0].set_ylabel('Surface solid NaCl (g m$^{-2}$)');finish(fig,'figure1_surface_precipitation')

fig,ax=plt.subplots(1,3,figsize=(14,4));composition_limits=[]
for n in colors:
    for fixed,style in ((False,'--'),(True,'-')):
        path=root/(f'correction/fixed_pre/fixed{n}' if fixed else f'diagnostic/historical{n}')
        f=read(path/'grid_profiles.csv');f=f[f.time_s==2500].groupby('depth_m').mean(numeric_only=True)
        near=f[f.index<=.002];composition_limits.extend([near.c.min(),near.c.max()])
        ax[int(fixed)].plot(f.c,f.index*1e3,style,color=colors[n],label=str(n))
        q=read(path/'grid_budgets.csv');q=q[q.time_s==2500].sort_values('depth_m')
        ax[2].plot(q.depth_m*1e3,q.supersaturation_excess_kg*1e6,style,color=colors[n],label=f'{n}, {"implicit" if fixed else "projection"}')
for a,title in zip(ax[:2],['Projection: 2500 s','Implicit: 2500 s']):
    a.axvline(.2649408697676898,color='k',lw=.8,label='$C_{eq}$');a.set(xlabel='Dissolved NaCl mass fraction',ylabel='Depth below soil surface (mm)',ylim=(2,0),title=title);a.legend(frameon=False)
    lo,hi=min(composition_limits),max(composition_limits);pad=max(1e-6,.05*(hi-lo));a.set_xlim(lo-pad,hi+pad)
ax[2].set(xlabel='Integration depth (mm)',ylabel='Supersaturated dissolved excess (g m$^{-2}$)',title='Actual BOX overlap inventories');ax[2].legend(fontsize=7,frameon=False)
finish(fig,'figure2_soil_profiles_inventory')

fig,ax=plt.subplots(1,2,figsize=(11,4))
for n in colors:
    for k,fixed in enumerate((False,True)):
        path=root/(f'correction/fixed_pre/fixed{n}' if fixed else f'diagnostic/historical{n}')
        q=read(path/'grid_budgets.csv');q=q[q.depth_m==.01]
        ax[k].plot(q.time_s,q.surface_solid_kg*1e6,color=colors[n],label=f'{n}: surface solid')
        ax[k].plot(q.time_s,(q.salt_after_kg-q.salt_before_kg.iloc[0])*1e6,'--',color=colors[n],label=f'{n}: dissolved change')
for a,title in zip(ax,['Projection','Implicit surface flux']):
    a.set(xlabel='Time (s)',ylabel='NaCl mass / surface area (g m$^{-2}$)',title=title+' — upper 10 mm',xlim=(0,2500));a.legend(fontsize=7,frameon=False,ncol=2)
finish(fig,'figure3_solid_dissolved')

a=read(root/'audit/activation_comparison.csv');fig,ax=plt.subplots(1,2,figsize=(10,4))
for stage,marker,label in [('ale_diagnostic','s','Historical'),('correction/fixed_coupled','o','Implicit')]:
    subset=a[a.stage.isin([stage,'historical_reference'])] if stage=='ale_diagnostic' else a[a.stage==stage]
    subset=subset[~subset['case'].str.contains('half|rollback')].copy()
    subset['n']=subset['case'].str.extract(r'(\d+)').astype(int)
    subset=subset.sort_values('n')
    ax[0].plot(subset.n,subset.c_initial,marker+'-',label=label)
    ax[1].plot(subset.n,subset.donor_depth_mm,marker+'-',label=label)
ax[0].axhline(.2649408697676898,color='k',ls=':',label='$C_{eq}$')
for p in ax:p.set_xlabel('Soil vertical layers');p.set_xticks([40,80,160,320]);p.legend(frameon=False)
ax[0].set_ylabel('Initial ALE liquid mass fraction');ax[1].set_ylabel('Deepest contributed BOX extent (mm)')
finish(fig,'figure4_initialization')

f=read(root/'ale_diagnostic/historical320_ale/ale_diagnostics.csv');accepted=f[f.status=='accepted'];start=accepted.time_s.min();last=accepted.time_s.max()
times=np.sort(accepted.time_s.unique());fig,axes=plt.subplots(1,3,figsize=(13,4),sharey=True)
for ax,key,title,cmap in zip(axes,['epsilon','c','reaction_kg_m3_s'],['Porosity','$C-C_{eq}$ (symmetric log scale)','Internal precipitation (kg m$^{-3}$ s$^{-1}$; symlog)'],['viridis','cividis','magma']):
    polygons=[];values=[]
    for k,t in enumerate(times[:-1]):
        g=accepted[accepted.time_s==t]
        for r in g.itertuples():
            polygons.append([(t-start,r.z_bottom_m*1e6),(times[k+1]-start,r.z_bottom_m*1e6),(times[k+1]-start,r.z_top_m*1e6),(t-start,r.z_top_m*1e6)])
            values.append(getattr(r,key)-(.2649408697676898 if key=='c' else 0))
    norm=SymLogNorm(linthresh=1e-7) if key=='c' else SymLogNorm(linthresh=1) if key=='reaction_kg_m3_s' else None
    coll=PolyCollection(polygons,array=np.asarray(values),cmap=cmap,norm=norm,edgecolors='none',rasterized=True)
    ax.add_collection(coll);ax.set(xlim=(0,10),ylim=(0,3.56),xlabel='Time since activation (s)',title=title)
    ax.axvspan(last-start,10,facecolor='.96',hatch='//',edgecolor='.7');ax.text(7.8,1.8,'No accepted\nstate',ha='center',fontsize=9)
    fig.colorbar(coll,ax=ax,shrink=.8)
axes[0].set_ylabel('Physical height above soil ($\mu$m)');finish(fig,'figure5_original_failure_fields')

fig,ax=plt.subplots(1,3,figsize=(12,4),sharey=True)
for status,style in [('last_accepted','-'),('first_rejected','--')]:
    f=read(root/f'audit/ale320_{status}.csv');z=(f.z_bottom_m+f.z_top_m)*.5e6
    for a,key in zip(ax,['epsilon','c','reaction_kg_m3_s']):a.plot(f[key],z,style,label=f'{status.replace("_"," ")} ({f.time_s.iloc[0]-start:.6f} s)')
for a,label in zip(ax,['Porosity','NaCl mass fraction','Internal precipitation (kg m$^{-3}$ s$^{-1}$)']):a.set_xlabel(label);a.legend(fontsize=7,frameon=False)
ax[1].xaxis.set_major_locator(MaxNLocator(3));ax[1].xaxis.set_major_formatter(FormatStrFormatter('%.6f'))
ax[0].axvline(.02,color='k',ls=':');ax[0].set_ylabel('Physical height above soil ($\mu$m)');finish(fig,'figure5b_accepted_rejected')

fig,ax=plt.subplots(1,2,figsize=(11,4))
for n in colors:
    h=hist(n,True,True);ax[0].plot(h.time_s,h.h_m*1e6,color=colors[n],label=f'{n} layers (ends {h.time_s.iloc[-1]:g} s)')
    f=read(root/f'correction/fixed_coupled/fixed{n}/ale_diagnostics.csv');f=f[(f.status=='accepted')&(f.time_s==2500)]
    if len(f):ax[1].plot(f.epsilon,(f.z_top_m+f.z_bottom_m)*.5e6,color=colors[n],label=str(n))
ax[0].set(xlabel='Time (s)',ylabel='Actual ALE height ($\mu$m)',xlim=(0,2500));ax[0].legend(fontsize=8,frameon=False)
ax[1].set(xlabel='Porosity at 2500 s',ylabel='Physical height above soil ($\mu$m)');ax[1].legend(frameon=False)
finish(fig,'figure6_coupled_height_porosity')

f=read(root/'correction/fixed_coupled/fixed320/ale_diagnostics.csv');f=f[f.status=='accepted'];available=np.sort(f.time_s.unique())
targets=[available[0],1200,1900,2500];selected=[available[np.argmin(abs(available-t))] for t in targets]
fig,ax=plt.subplots(1,4,figsize=(10,4.5),sharey=True,layout='constrained');norm=Normalize(.02,.43)
for a,t in zip(ax,selected):
    g=f[f.time_s==t];polygons=[]
    for r in g.itertuples():polygons.append([(0,r.z_bottom_m*1e6),(1,r.z_bottom_m*1e6),(1,r.z_top_m*1e6),(0,r.z_top_m*1e6)])
    coll=PolyCollection(polygons,array=g.epsilon.to_numpy(),cmap='viridis',norm=norm,edgecolors='none');a.add_collection(coll)
    a.set(xlim=(-.15,1.15),ylim=(0,max(f.h_m)*1e6*1.08),xticks=[],title=f'{t:g} s\nh = {g.h_m.iloc[0]*1e6:.3f} $\mu$m');a.axhline(0,color='#826B50',lw=3)
ax[0].set_ylabel('Physical height above soil ($\mu$m)');fig.colorbar(coll,ax=list(ax),label='Local ALE porosity',shrink=.85)
finish(fig,'figure7_physical_columns')

fig,ax=plt.subplots(1,2,figsize=(11,4))
independent=read(root/'audit/independent_residual_series.csv')
for n in colors:
    h=hist(n,True,True);a=h.iloc[0]
    for p,key,den in zip(ax,['salt_residual_kg','water_residual_kg'],[a.soil_salt_kg,a.soil_water_kg]):p.plot(h.time_s,h[key]/den,color=colors[n],label=f'Implicit {n}')
    q=independent[(independent.stage=='correction/fixed_coupled')&(independent['case']==f'fixed{n}')]
    for p,key in zip(ax,['salt_relative','water_relative']):p.scatter(q.time_s,q[key],s=12,facecolors='none',edgecolors=colors[n])
h=read(root/'ale_diagnostic/historical320_ale/coupled_history.csv');a=h.iloc[0]
for p,key,den,title in zip(ax,['salt_residual_kg','water_residual_kg'],[a.soil_salt_kg,a.soil_water_kg],['NaCl','Water']):
    p.plot(h.time_s,h[key]/den,'--',color='k',label='Historical 320 (failed)');p.ticklabel_format(axis='y',style='sci',scilimits=(0,0));p.set(xlabel='Accepted time (s)',ylabel='Relative combined mass residual',title=title);p.legend(fontsize=8,frameon=False)
finish(fig,'figure8_conservation')

fig,ax=plt.subplots(1,2,figsize=(11,4));table=read(root/'audit/ablation_comparison.csv')
q=table[(table['case'].str.startswith('eq0'))&(table.requested_s==5)].copy();labels=q['case'].str.replace('eq0_','').str.replace('replay_','').str.replace('_half',' dt/2').str.replace('_',' ')
ax[0].bar(labels,q.internal_dynamic_g_m2,color='#0072B2');ax[0].set_ylabel('Internal precipitation at 5 s (g m$^{-2}$)')
ax[1].bar(labels,q.epsilon_min,color='#009E73');ax[1].axhline(.02,color='k',ls=':');ax[1].set_ylabel('Minimum porosity at 5 s')
for a in ax:a.tick_params(axis='x',rotation=35)
finish(fig,'figure9_matched_state_ablations')

old=read(root/'audit/diagnostic_summary.csv');new=read(root/'correction/audit/fixed_pre_summary.csv')
fig,ax=plt.subplots(1,2,figsize=(11,4))
for a,token,title in zip(ax,['familyA','familyB160_top'],['Fixed overall graded mapping','Fixed deep mesh; surface refinement']):
    q=old[old['case'].str.startswith(token)].sort_values('top_cell_m')
    token2='fixedA' if token=='familyA' else 'fixedB160_top'
    r=new[new['case'].str.startswith(token2)].sort_values('top_cell_m')
    a.semilogx(q.top_cell_m*1e6,q.surface_g_m2,'o--',label='Projection')
    a.semilogx(r.top_cell_m*1e6,r.surface_g_m2,'s-',label='Implicit surface flux')
    a.set(xlabel='Top soil cell thickness ($\mu$m)',ylabel='Surface NaCl at 2500 s (g m$^{-2}$)',title=title);a.legend(frameon=False)
finish(fig,'figure10_controlled_meshes')
