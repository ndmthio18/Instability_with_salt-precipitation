"""Actual simulation figures; model comparisons are conditional until all gates pass."""
from pathlib import Path
import sys,json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from analyze import ROOT,R,read as read_raw,PHI,CR,ELL
sys.path.insert(0,str(ROOT/'qa_tools'))
from audit_panel_alignment import require_matplotlib_panel_alignment
F=ROOT/'figures';F.mkdir(exist_ok=True)
plt.rcParams.update({'font.family':'serif','font.serif':['Times New Roman'],'mathtext.fontset':'custom','mathtext.rm':'Times New Roman','mathtext.it':'Times New Roman:italic','mathtext.bf':'Times New Roman:bold','font.size':8,'axes.labelsize':9,'legend.fontsize':8,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none','pdf.fonttype':42,'axes.linewidth':.7,'lines.linewidth':1.2})
colors=['#0072B2','#D55E00','#009E73','#CC79A7','#555555']

def layout(nr=2,nc=2,height=135):
 height=max(height,100) if nr==1 else height
 fig,ax=plt.subplots(nr,nc,figsize=(180/25.4,height/25.4),squeeze=False)
 fig.subplots_adjust(left=.13,right=.975,bottom=.25 if nr==1 else .17,top=.90 if nr==1 else .94,wspace=.47,hspace=.52)
 for i,a in enumerate(ax.flat):
  a.annotate(chr(97+i),xy=(0,1),xycoords='axes fraction',xytext=(-24,7),textcoords='offset points',fontweight='bold',fontsize=10,annotation_clip=False)
  a.tick_params(direction='out',length=3,width=.6);a.grid(alpha=.16,linewidth=.5)
 return fig,ax

def save(fig,name,handles=None,labels=None,ncol=3):
 if name.startswith('figure6') and not (ROOT/'raw/p2_u06_n256_k0.5_history.csv').exists():plt.close(fig);return
 if name.startswith('figure8') and not (ROOT/'raw'/f'{primary("grid_u06_n512")}_history.csv').exists():plt.close(fig);return
 if handles:fig.legend(handles,labels,loc='lower center',bbox_to_anchor=(.52,.015),frameon=False,ncol=ncol)
 fig.canvas.draw()
 require_matplotlib_panel_alignment(fig,json_out=str(F/(name+'.alignment.json')),tolerance_pt=1.5,gutter_tolerance_pt=1.5,strict=True)
 fig.savefig(F/(name+'.png'),dpi=600)
 fig.savefig(F/(name+'.pdf'))
 fig.savefig(F/(name+'.svg'))
 plt.close(fig)

def primary(c):
 candidate='fast_'+c.removeprefix('setup_')
 return candidate if (ROOT/'raw'/f'{candidate}_history.csv').exists() else c

def read(c):return read_raw(primary(c))

def choose(lab):
 for c in [f'production_{lab}_n512',f'grid_{lab}_n512',f'p2_{lab}_n256_k1',f'setup_p0_ld_{lab}_r1e-10',f'p0_ld_{lab}_r1e-10']:
  if (ROOT/'raw'/f'{primary(c)}_history.csv').exists():return primary(c)
 return None

def main():
 fig,ax=layout();handles=[];labels=[]
 names=['legacy_r9','legacy_r10','ablate_physical','fast_p0_ld_u06_r1e-09']
 for k,c in enumerate(names):
  p=ROOT/'raw'/f'{c}_diagnostics.csv'
  if not p.exists():continue
  d=pd.read_csv(p);d=d[(d.internal_t>0)&(d.internal_t<=1000)]
  for a,col,fac in zip(ax.flat,['last_step','worst_weighted_error','jacobians','liquid_cell_min'],[1,1,1,1e9]):
   h,=a.plot(d.internal_t,d[col]*fac,color=colors[k],marker='o',markersize=2);a.set_xscale('log');a.set_yscale('log');a.set_xlabel('Internal time (s)')
   if (d.flag<0).any():a.scatter(d.internal_t.iloc[-1],d[col].iloc[-1]*fac,color=colors[k],marker='x',s=35,zorder=4)
  handles.append(h);labels.append(['Original rtol 1e−9','Original rtol 1e−10','Full-Q error weights','Improved rtol 1e−9'][k])
 for a,label in zip(ax.flat,['Last accepted step (s)','Largest weighted component','Jacobian evaluations','Smallest liquid cell (nm)']):a.set_ylabel(label)
 save(fig,'figure1_numerical_stability',handles,labels,2)

 fig,ax=layout(3,2,180)
 for j,lab in enumerate(['u06','x008']):
  for i,key in enumerate(['Fplus','Fminus','Cb']):
   for k,n in enumerate([128,256,512]):
    for suffix,ls in [('b1e-09','-'),('tight','--')]:
     c=f'setup_early_{lab}_n{n}_{suffix}'
     if not (ROOT/'raw'/f'{c}_history.csv').exists():c=c.removeprefix('setup_')
     h=read(c);h=h[h.t>0];v=h[key]*1e6 if key!='Cb' else (h.Cb-.264940869768)*1e6
     ax[i,j].plot(h.t,v,color=colors[k],ls=ls);ax[i,j].set_xscale('log');ax[i,j].set_xlabel('Time (s)')
   ax[i,j].set_ylabel({'Fplus':'Crust diffusion (mg m⁻² s⁻¹)','Fminus':'Liquid diffusion (mg m⁻² s⁻¹)','Cb':'Bottom offset (10⁻⁶)'}[key])
  ax[0,j].set_title('u = 0.6' if lab=='u06' else 'x008',fontsize=9)
 handles=[Line2D([],[],color=colors[i],label=f'{n} cells') for i,n in enumerate([128,256,512])]+[Line2D([],[],color='k',ls='--')]
 save(fig,'figure2_early_flux_convergence',handles,['128 cells','256 cells','512 cells','Refined time integration'],4)

 cases=[(lab,choose(lab)) for lab in ['u00','u03','u06','u09','x008']];cases=[x for x in cases if x[1]]
 fig,ax=layout(1,3,85);handles=[];labels=[]
 for k,(lab,c) in enumerate(cases):
  h=read(c)
  for a,key in zip(ax.flat,['zt','zb','hc']):line,=a.plot(h.t/3600,h[key]*1e6,color=colors[k]);a.set_xlabel('Time (h)')
  ax[0,2].plot(h.t/3600,h.hc_corr*1e6,color=colors[k],ls='--')
  handles.append(line);labels.append({'u00':'u = 0','u03':'u = 0.3','u06':'u = 0.6','u09':'u = 0.9','x008':'x008'}[lab])
 for a,label in zip(ax.flat,['Top position (µm)','Bottom position (µm)','Geometric thickness (µm)']):a.set_ylabel(label)
 save(fig,'figure3_double_fronts',handles+[Line2D([],[],color='k',ls='--')],labels+['Analytical thickness'],3)

 fig,ax=layout();handles=[];labels=[]
 for k,(lab,c) in enumerate(cases):
  h=read(c)
  for j,key in enumerate(['vt','vnet']):
   line,=ax[0,j].plot(h.t/3600,h[key]*1e9,color=colors[k]);ax[0,j].plot(h.t/3600,h[key+'_corr']*1e9,color=colors[k],ls='--')
   ax[1,j].plot(h.t/3600,(h[key]-h[key+'_corr'])/h.vref*100,color=colors[k])
  handles.append(line);labels.append(lab)
 for i,a in enumerate(ax.flat):a.set_xlabel('Time (h)');a.set_ylabel(['Top speed (nm s⁻¹)','Net speed (nm s⁻¹)','Top error / reference (%)','Net error / reference (%)'][i])
 save(fig,'figure4_growth_comparison',handles+[Line2D([],[],color='k',ls='--')],labels+['Analytical'],3)

 fig,ax=layout(1,3,85)
 for a,lab in zip(ax.flat,['u06','u09','x008']):
  h=read(choose(lab))
  for k,key in enumerate(['delta_top','delta_bottom','delta_net']):a.plot(h.t/3600,h[key]/h.vref*100,color=colors[k])
  a.set_xlabel('Time (h)');a.set_ylabel('Signed discrepancy / reference (%)');a.set_title({'u06':'u = 0.6','u09':'u = 0.9','x008':'x008'}[lab],fontsize=9)
 save(fig,'figure5_error_cancellation',[Line2D([],[],color=colors[i]) for i in range(3)],['Top error','Bottom diagnostic','Net error'],3)

 fig,ax=layout(2,3,135)
 for row,lab in enumerate(['u06','u09']):
  for k,kappa in enumerate([.5,1,2]):
   c=f'p2_{lab}_n256_k{kappa:g}'
   if not (ROOT/'raw'/f'{c}_history.csv').exists():continue
   h=read(c)
   for a,key,scale in zip(ax[row],['vb','zb','hc'],[1e9,1e6,1e6]):a.plot(h.t/3600,h[key]*scale,color=colors[k]);a.set_xlabel('Time (h)')
  for a,label in zip(ax[row],['Bottom speed (nm s⁻¹)','Bottom position (µm)','Thickness (µm)']):a.set_ylabel(label)
  ax[row,0].set_title({'u06':'u = 0.6','u09':'u = 0.9','x008':'x008'}[lab],fontsize=9)
 save(fig,'figure6_bottom_kinetics',[Line2D([],[],color=colors[i]) for i in range(3)],['κ = 0.5','κ = 1','κ = 2'],3)

 fig,ax=layout();handles=[];labels=[]
 for k,(lab,c) in enumerate(cases):
  h=read(c);line,=ax[0,0].plot(h.t/3600,h.phi_mean,color=colors[k]);ax[0,1].plot(h.t/3600,h.core_phi,color=colors[k]);ax[1,1].plot(h.t/3600,h.phib,color=colors[k])
  c=primary(c);pp=ROOT/'profiles_selected'/f'{c}_profiles.csv'
  if not pp.exists():pp=ROOT/'raw'/f'{c}_profiles.csv'
  if pp.exists():
   p=pd.read_csv(pp);p=p[(p.t==p.t.max())&(p.region=='crust')];center=(p.zlo+p.zhi)/2;pos=(center-p.zlo.min())/(p.zhi.max()-p.zlo.min());ax[1,0].plot(pos,p.phi,color=colors[k])
  handles.append(line);labels.append(lab)
 for a in [ax[0,0],ax[0,1],ax[1,0]]:a.axhline(PHI,color='k',ls='--',lw=1)
 for a in [ax[0,0],ax[0,1],ax[1,1]]:a.set_xlabel('Time (h)')
 ax[1,0].set_xlabel('Position within active crust / thickness')
 for a,label in zip(ax.flat,['Volume-mean solid fraction','Core-region solid fraction','Final local solid fraction','Solid fraction ahead of bottom']):a.set_ylabel(label)
 save(fig,'figure7_core_fraction',handles+[Line2D([],[],color='k',ls='--')],labels+['Analytical core'],3)

 fig,ax=layout();handles=[];labels=[]
 for k,lab in enumerate(['u06','u09','x008']):
  names=[f'setup_p0_ld_{lab}_r1e-09',f'p2_{lab}_n256_k1' if lab!='x008' else 'grid_x008_n256',f'grid_{lab}_n512']
  names=[primary(c) for c in names]
  if not all((ROOT/'raw'/f'{x}_history.csv').exists() for x in names):continue
  for j in range(2):
   a,b=map(read,names[j:j+2]);t=np.intersect1d(a.t,b.t);a=a.set_index('t').loc[t];b=b.set_index('t').loc[t]
   for axis,key in zip(ax.flat,['vt','vb','vnet','hc']):
    den=b.vref if key!='hc' else b.hc
    line,=axis.plot(t/3600,abs(a[key]-b[key])/den*100,color=colors[k],ls='-' if j==0 else '--');axis.set_xlabel('Time (h)')
  handles.append(line);labels.append(lab)
 for a,label in zip(ax.flat,['Top difference / reference (%)','Bottom difference / reference (%)','Net difference / reference (%)','Thickness difference (%)']):a.set_ylabel(label);a.axhline(1,color='0.4',ls=':',lw=.8)
 save(fig,'figure8_numerical_convergence',handles+[Line2D([],[],color='k'),Line2D([],[],color='k',ls='--')],labels+['128 → 256','256 → 512'],3)
if __name__=='__main__':main()
