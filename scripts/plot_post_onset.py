"""WP4 evidence-led figures; all curves show existing observations only."""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyze_post_onset import *
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':140,'savefig.dpi':220,'pdf.fonttype':42})
COL={'resolved_onset':'#0072B2','weak_interval_onset':'#E69F00','immediate_growth':'#009E73','right_censored':'#CC79A7'}
SHORT={'resolved_onset':'Resolved','weak_interval_onset':'Weak interval*','immediate_growth':'Immediate growth*','right_censored':'Right censored*'}
def save(fig,name,note):
 fig.text(.01,.012,note,fontsize=8,wrap=True);fig.tight_layout(rect=(0,.065,1,1));fig.savefig(OUT/(name+'.png'));fig.savefig(OUT/(name+'.pdf'));plt.close(fig)
def xy(rows,x,y):
 a=np.array([(num(r.get(x)),num(r.get(y))) for r in rows]);m=np.all(np.isfinite(a),axis=1);return a[m].T

def main():
 rows=readcsv(OUT/'run_metrics.csv');main=[r for r in rows if r['source_group']=='main'];raw=[r for r in main if r['series_status']=='hash_verified'];res=[r for r in raw if r['dumux_detection_class']=='resolved_onset'];supp=[r for r in rows if r['source_group']=='wp1' and r['dumux_detection_class']=='resolved_onset']
 fig,axs=plt.subplots(1,2,figsize=(10,4.3))
 for cls,c in COL.items():
  a=sorted(num(r['frozen_G_end']) for r in main if r['dumux_detection_class']==cls and np.isfinite(num(r['frozen_G_end'])))
  axs[0].step(a,np.arange(1,len(a)+1)/len(a),where='post',color=c,label=f'{SHORT[cls]} (n={len(a)})')
 axs[0].set(xlabel='Endpoint / detector-associated minimum',ylabel='Cumulative fraction',title='436 frozen run summaries');axs[0].legend(fontsize=8);axs[0].set_xscale('log')
 for rr,key,label,ls in [(res,'G_end','Main endpoint (234)','-'),(res,'G_max','Main post-detect maximum (234)','--'),(supp,'G_end','WP1 endpoint (70)',':')]:
  a=sorted(num(r[key]) for r in rr);axs[1].step(a,np.arange(1,len(a)+1)/len(a),where='post',label=label,ls=ls)
 axs[1].set(xlabel='Amplification factor',ylabel='Cumulative fraction',title='Full-series resolved subset');axs[1].legend(fontsize=8);axs[1].set_xscale('log')
 save(fig,'figure1_amplification','* Non-resolved ratios use exploratory minima, not resolved detection times. Main G_end and G_max overlap; WP1 is separate.')
 fig,axs=plt.subplots(1,2,figsize=(10,4.3))
 for ax,x,label in [(axs[0],'post_window_s','Post-detection window (days)'),(axs[1],'end_over_detect','End time / detected onset')]:
  for warned,color,l in [('False','#0072B2','Mean remains < 0.1'),('True','#D55E00','Mean reaches 0.1 (invalid late physics)')]:
   rr=[r for r in res if r['mean_threshold_warning']==warned];xx,yy=xy(rr,x,'G_end');ax.scatter(xx/86400 if x=='post_window_s' else xx,yy,s=18,alpha=.65,c=color,label=l)
  ax.set(xlabel=label,ylabel='G_end',title='Main: 234 resolved available traces')
 axs[0].set_xscale('log');axs[1].legend(fontsize=7)
 save(fig,'figure2_observation_window','Descriptive run-level associations; repeated physical cases and varying initialization/window policies are not independent samples.')
 fig,axs=plt.subplots(2,3,figsize=(11,7))
 for ax,(x,label,log) in zip(axs.flat,[('effective_saturation_top','Initial top effective saturation',False),('permeability_m2','Permeability (m²)',True),('evaporation_m_s','Evaporation rate (m/s)',True),('r_e','R_E',True),('nx_over_global_critical','Initialized nx / global-critical nx',True),('end_over_detect','End / detected onset',False)]):
  xx,yy=xy(res,x,'G_end');ax.scatter(xx,yy,s=15,alpha=.6,color='#0072B2');ax.set(xlabel=label,ylabel='G_end');
  if log:ax.set_xscale('log')
 save(fig,'figure3_physical_associations','Recorded-input associations only: incomplete historical executable provenance, coupled parameters, protocols and repeated cases prevent causal inference.')
 fig,ax=plt.subplots(figsize=(7,4.7))
 for r in raw:r['G_plot']=r.get('G_end') or r.get('frozen_G_end')
 for cls,c in COL.items():
  rr=[r for r in raw if r['dumux_detection_class']==cls];xx,yy=xy(rr,'G_plot','A_rel_end');ax.scatter(xx,100*yy,label=SHORT[cls],color=c,s=22,alpha=.65)
 ax.axhline(.0196581,color='.4',ls=':',label='Initial relative magnitude (~0.01966%)');ax.set(xlabel='Endpoint / detector-associated minimum',ylabel='Final sigma / mean (%)',title='All 320 available main traces: small final perturbations',ylim=(0,.021));ax.legend(fontsize=8,ncol=2)
 save(fig,'figure4_relative_magnitude','* Non-resolved ratios are exploratory. No nonlinear threshold is defined; 116 unavailable traces have no relative-magnitude estimate.')
 dense=readcsv(OUT/'dense_summary.csv');fig,axs=plt.subplots(1,2,figsize=(11,5))
 labels=[r['id'].replace('r','',1) for r in dense];y=np.arange(len(dense))
 axs[0].plot([num(r['depth_centroid_initial_H']) for r in dense],y,'o',label='Initial');axs[0].plot([num(r['depth_centroid_last_pre_H']) for r in dense],y,'x',label='Last saved pre-threshold');axs[0].set_yticks(y,labels,fontsize=8);axs[0].set(xlabel='Depth-energy centroid / H',title='14 histories: broad modes start broad');axs[0].legend(fontsize=8)
 for rid,col in [('r05_c36877_n30','#0072B2'),('r06_c36877_n60','#D55E00')]:
  pp=np.load(ROOT/'analysis/dense_history_20261102/evidence/runs'/rid/'profiles.npz');depth=1-pp['z_m']/pp['z_m'][-1]
  axs[1].plot(pp['x_sigma'][0]*1e6,depth,ls=':',color=col,label=rid+' initial');axs[1].plot(pp['x_sigma'][-1]*1e6,depth,color=col,label=rid+' final')
 axs[1].set(xlabel='Lateral sigma (10^-6 mol/mol)',ylabel='Depth below surface / H',ylim=(.02,0),title='36877: top 2% depth shown explicitly');axs[1].legend(fontsize=7)
 save(fig,'figure5_depth_structure','Energy centroid uses sigma_x(z)^2 with nodal depth weights; it is not a finger tip or a flux. Snapshots only; no interpolated dynamics.')
 times=readcsv(OUT/'dense_depth_timeseries.csv');inventory=readcsv(ROOT/'analysis/dense_history_20261102/evidence/runs/r06_c36877_n60/inventory_timeseries.csv');fig,axs=plt.subplots(1,3,figsize=(12,4.2))
 for rid,col in [('r05_c36877_n30','#0072B2'),('r06_c36877_n60','#D55E00'),('r07_c36877_n5','#009E73'),('r08_c36877_n5','#CC79A7')]:
  rr=[r for r in times if r['id']==rid];t=np.array([num(r['time_s']) for r in rr])/86400
  axs[0].plot(t,[num(r['top_sigma'])*1e6 for r in rr],label=rid,color=col);axs[1].plot(t,[num(r['depth_centroid_H']) for r in rr],color=col)
 for region in ['upper_0.9','upper_0.75','upper_0.5']:
  rr=[r for r in inventory if r['region']==region];v=np.array([num(r['inventory_salt_m']) for r in rr]);axs[2].plot(np.array([num(q['time_s']) for q in times if q['id']=='r06_c36877_n60'])/86400,(v/v[0]-1)*100,label=region)
 axs[0].set(ylabel='Top sigma (10^-6 mol/mol)',title='Separate existing runs');axs[0].legend(fontsize=7);axs[1].set(ylabel='Depth-energy centroid / H',title='Small depth changes');axs[2].set(ylabel='Inventory change from initial (%)',title='nx=60: upper inventory increases');axs[2].legend(fontsize=7)
 for ax in axs:ax.set_xlabel('Physical time (days)')
 save(fig,'figure6_dense_evolution','Different nx, meshes and end times; curves are not one trajectory. Inventory change without matched M0 cannot isolate convective transport.')
 # Honest concentration anomaly map: common absolute scale, explicit shallow crop.
 data=np.load(ROOT/'analysis/dense_history_20261102/evidence/runs/r06_c36877_n60/representative_fields.npz');fig,axs=plt.subplots(1,2,figsize=(10,4),layout='constrained');z=data['z_m'];mask=z/z[-1]>=.99;vmax=max(abs(data['x_anomaly_'+str(j)][mask]).max() for j in [0,100])*1e6
 for ax,j in zip(axs,[0,100]):
  im=ax.pcolormesh(data['xs'],1-z[mask]/z[-1],data['x_anomaly_'+str(j)][mask]*1e6,cmap='RdBu_r',vmin=-vmax,vmax=vmax,shading='auto');ax.invert_yaxis();ax.set(xlabel='Horizontal coordinate (m)',ylabel='Depth / H',title=('Initial' if j==0 else 'Final')+' — top 1% depth')
 fig.colorbar(im,ax=axs,label='x - horizontal mean (10^-6 mol/mol)',shrink=.8);fig.suptitle('36877 nx=60: persistent imposed wavelength; small absolute anomaly',fontsize=11);fig.savefig(OUT/'figure7_field_anomalies.png');fig.savefig(OUT/'figure7_field_anomalies.pdf');plt.close(fig)
 print('7 PNG/PDF pairs written')
if __name__=='__main__':main()
