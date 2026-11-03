"""Reproduce baseline Pilot-0 evidence from completed, hash-verified sal pair.
No simulation launch; no interpolation of budgets; field PVD rounding is checked
against accepted-step timestamps before using the exact budget time.
"""
import csv, hashlib, importlib.util, json, math, tarfile, xml.etree.ElementTree as ET
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyze_matched_pair import validate_pair, rows, REGIONS

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'analysis/pilot0_integrated_20261102'
RUN=OUT/'sal_migration/runs'
FIG=OUT/'figures';FIG.mkdir(exist_ok=True)

def dump(path,obj):path.write_text(json.dumps(obj,indent=2,allow_nan=False)+'\n',encoding='utf-8')
def table(path,data):
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

validate_pair(RUN/'M0',RUN/'M1')
pair=rows(OUT/'sal_baseline_pair/pair_timeseries.csv')
(OUT/'transport_timeseries.csv').write_bytes((OUT/'sal_baseline_pair/pair_timeseries.csv').read_bytes())
regions={key:[r for r in pair if r['region']==key] for key in REGIONS}
base=regions[REGIONS[0]]
t=np.array([r['time_s'] for r in base]);T=float(t[-1]);hours=t/3600
assert all(r['M0_time_s']==r['M1_time_s'] for r in pair)
source=OUT/'provenance/extract_onset_2025.py'
assert sha(source)=='51e75a7f6225b2526de584ea9865208457a3d818a7c3b381cadc09af3f67fbe0'
spec=importlib.util.spec_from_file_location('frozen_extractor',source);ext=importlib.util.module_from_spec(spec);spec.loader.exec_module(ext)
sig={role:np.loadtxt(next((RUN/role).glob('*sigmatop.csv'))) for role in ('M0','M1')}
detection=ext.extract_onset(str(next((RUN/'M1').glob('*sigmatop.csv'))))
detection.update(extractor_sha256=sha(source),source_commit='158395f3d25a',interpretation='Recomputed from actual new M1; interval is plateau localization, not confidence interval')
dump(OUT/'new_M1_detection.json',detection)
archive=json.loads((OUT/'sal_migration/field_archive.json').read_text());integrity={}
for role in ('M0','M1'):
    finished=json.loads((RUN/role/'finished.json').read_text())
    checked=[]
    for name,record in finished['files'].items():
        p=RUN/role/name
        if p.exists():assert sha(p)==record['sha256'];checked.append(name)
        else:assert archive[role]['files'][name]==record['sha256'];checked.append(name)
    assert archive[role]['verified']
    assert sha(Path(archive[role]['local_path']))==archive[role]['sha256']
    integrity[role]=dict(files_verified=len(checked),wall_s=finished['wall_s'],finished_utc=finished['finished_utc'],archive_sha256=archive[role]['sha256'])

# All saved fields, using one identical nx28 basis for M0 and M1.
profile=np.loadtxt(RUN/'M1/eigenprofile.dat');field_rows=[];depth_rows=[];selected={};time_checks=[]
for role in ('M0','M1'):
    entries=list(ET.parse(next((RUN/role).glob('*.pvd'))).iter('DataSet'))
    assert len(entries)==13
    with tarfile.open(archive[role]['local_path'],'r:gz') as tar:
        for j,item in enumerate(entries):
            name=item.get('file');index=int(name.rsplit('-',1)[1].split('.')[0]);step=index*50
            assert j==index and abs(float(item.get('timestep'))-t[step])<=.5
            time_checks.append(dict(role=role,file=name,pvd_rounded_time_s=float(item.get('timestep')),accepted_step=step,exact_time_s=float(t[step])))
            with tar.extractfile(name) as stream:root=ET.parse(stream)
            arrays={a.get('Name'):np.fromstring(a.text or '',sep=' ') for a in root.iter('DataArray') if a.get('Name') in ('Coordinates','moleFraction') and a.get('format')=='ascii'}
            xyz=arrays['Coordinates'].reshape(-1,3);order=np.lexsort((xyz[:,0],xyz[:,1]));xyz=xyz[order]
            xs=np.unique(xyz[:,0]);zs=np.unique(xyz[:,1]);x=arrays['moleFraction'][order].reshape(len(zs),len(xs))
            assert x.shape==(201,601)
            mean=x.mean(axis=1);anomaly=x-mean[:,None]
            basis=np.interp(zs,profile[:,0],profile[:,1])[:,None]*np.cos(2*np.pi*28*xs/(xs[-1]-xs[0]))[None,:]
            projection=float(abs(np.sum(anomaly*basis))/np.sum(basis*basis))
            if role=='M1':assert abs(projection-sig[role][step,5])<1e-14
            r=dict(role=role,time_s=float(t[step]),step=step,common_nx28_projection=projection,top_sigma_from_field=float(np.std(x[-1])),max_abs_horizontal_anomaly=float(abs(anomaly).max()))
            for key in REGIONS:
                cut=regions[key][0]['actual_z_over_H']
                r['max_horizontal_anomaly_below_'+key]=float(abs(anomaly[zs<=cut*zs[-1]]).max())
            field_rows.append(r)
            for k,z in enumerate(zs):depth_rows.append(dict(role=role,time_s=float(t[step]),z_m=float(z),depth_below_top_m=float(zs[-1]-z),x_row_mean=float(mean[k]),horizontal_sigma=float(np.std(x[k])),max_abs_horizontal_anomaly=float(abs(anomaly[k]).max())))
            if role=='M1' and j in (0,4,12):selected[j]=(xs,zs,x,anomaly)
table(OUT/'field_common_mode_timeseries.csv',field_rows);table(OUT/'field_depth_profiles.csv',depth_rows);table(OUT/'field_time_mapping.csv',time_checks)

summary=dict(status='baseline_complete_transport_inconclusive',integrity=integrity,detection=detection,regions={},growth={},uncertainty=dict(U_dt=None,U_z=None,U_cut=None,U_x=None,time_alignment='identical accepted times; no interpolation',budget_envelope='abs(R0)+abs(R1), algebraic closure diagnostic only; not total numerical uncertainty'),transport_t_eff_s=None)
for key,rr in regions.items():
    dc=np.array([r['delta_C_m'] for r in rr]);env=np.array([r['budget_envelope_m'] for r in rr]);iref=rr[0]['reference_inventory_m']
    summary['regions'][key]=dict(actual_z_over_H=rr[0]['actual_z_over_H'],final=rr[-1],max_abs_delta_C_m=float(abs(dc).max()),max_abs_delta_C_over_Iref=float(abs(dc).max()/iref),max_budget_envelope_m=float(env.max()),max_identity_error_m=max(abs(r['identity_error_m']) for r in rr),samples_above_instantaneous_budget_envelope=int(np.sum(abs(dc)>env)),net_flux_M0_upward_every_step=all(r['M0_total_down_m_s']<0 for r in rr),net_flux_M1_upward_every_step=all(r['M1_total_down_m_s']<0 for r in rr))
for fraction in (.25,.5):
    a=sig['M1'];b=a[a[:,0]>=T*(1-fraction)];slope=float(np.polyfit(b[:,0],np.log(b[:,2]),1)[0]);end=float(a[-1,2])
    summary['growth'][str(fraction)]=dict(fit_fraction=fraction,slope_per_s=slope,doubling_s=math.log(2)/slope,conditional_sigma_at_2T=end*math.exp(slope*T),conditional_total_time_for_10x_current_sigma_s=T+math.log(10)/slope)
summary['surface']=dict(max_domain_x=max(r['M1_x_domain_max'] for r in base),final_mean_difference=base[-1]['M1_x_top_area_mean']-base[-1]['M0_x_top_area_mean'],final_max_difference=base[-1]['M1_x_top_max']-base[-1]['M0_x_top_max'],final_sigma_over_mean=base[-1]['M1_sigma_top_node']/base[-1]['M1_x_top_area_mean'],threshold_observed=False)
summary['cost_scenarios']={str(mult):{role:mult*integrity[role]['wall_s']/3600 for role in ('M0','M1')} for mult in (2,3,4)}
summary['cost_scope']='Linear runtime scaling scenarios from full baseline, not guaranteed forecasts. No validated restart.'
dump(OUT/'baseline_analysis.json',summary)

# Provisional internal research figures, no journal-specific compliance claim.
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'pdf.fonttype':42,'savefig.facecolor':'white'})
blue,orange,green='#0072B2','#D55E00','#009E73'
exports=[]
def export(fig,name):
    for suffix in ('png','pdf'):
        path=FIG/(name+'.'+suffix);fig.savefig(path,dpi=200);exports.append(dict(path=str(path.relative_to(ROOT)),sha256=sha(path)))
    plt.close(fig)
def references(ax):
    for value,color,style in [(31409.157792677746,'#666666',':'),(45474.35459078279,'#CC79A7','--'),(detection['onset_refined_s'],'#009E73','-.')]:ax.axvline(value/3600,color=color,ls=style,lw=1)
    ax.axvspan(detection['onset_low_s']/3600,detection['onset_high_s']/3600,color='#009E73',alpha=.09)
    ax.grid(alpha=.15)
fig,ax=plt.subplots(3,1,figsize=(9,9),sharex=True,layout='constrained')
for role,color,ls in [('M0',blue,'--'),('M1',orange,'-')]:ax[0].semilogy(hours,[r[role+'_sigma_top_node'] for r in base],color=color,ls=ls,label=role)
ax[0].set_ylabel('Surface node SD (mol/mol)');ax[0].legend(loc='center right');ax[0].set_title('A  Instability detection and surface evolution')
for k,label,ls in [(2,'Surface SD','-'),(4,'Surface Fourier amplitude','--'),(5,'Canonical mode projection',':')]:ax[1].plot(hours,sig['M1'][:,k]/sig['M1'][0,k],ls,label=label)
ax[1].set_ylabel('M1 amplitude / initial value');ax[1].legend(fontsize=9)
ax[2].plot(hours,[r['M0_x_top_area_mean'] for r in base],color=blue,ls='--',label='M0 mean');ax[2].plot(hours,[r['M1_x_top_area_mean'] for r in base],color=orange,label='M1 mean');ax[2].plot(hours,[r['M1_x_top_max'] for r in base],color=green,ls=':',label='M1 maximum');ax[2].set(xlabel='Time (h)',ylabel='Surface x (mol/mol)');ax[2].legend(fontsize=9)
for a in ax:references(a)
fig.supxlabel('Vertical lines: LSA critical (gray), LSA observable (pink), new detection (green).\nShading: detection plateau interval; not a confidence interval. No t_eff defined.',fontsize=9)
export(fig,'Figure_A_timeline')
fig,axes=plt.subplots(3,3,figsize=(14,10),sharex=True,layout='constrained')
for i,(key,rr) in enumerate(regions.items()):
    get=lambda k:np.array([r[k] for r in rr]);label=f"z/H = {rr[0]['actual_z_over_H']:.6f}"
    for role,color,ls in [('M0',blue,'--'),('M1',orange,'-')]:axes[i,0].plot(hours,get(role+'_cumulative_total_down_m')*1e6,color=color,ls=ls,label=role)
    axes[i,0].set_ylabel(label+'\nC (10$^{-6}$ m)');axes[i,0].set_title('Total cumulative flux; down positive' if i==0 else '')
    axes[i,1].fill_between(hours,-get('budget_envelope_m')*1e17,get('budget_envelope_m')*1e17,color='0.85',label='Budget closure envelope')
    axes[i,1].plot(hours,get('delta_C_m')*1e17,color=orange,label='Delta C');axes[i,1].plot(hours,-get('corrected_delta_inventory_m')*1e17,color=blue,ls='--',label='- corrected Delta I');axes[i,1].set_ylabel('Difference (10$^{-17}$ m)')
    for k,label2,color,ls in [('delta_adv_C_m','Advective',orange,'-'),('delta_diff_C_m','Diffusive',blue,'--')]:axes[i,2].plot(hours,get(k)*1e17,color=color,ls=ls,label=label2)
    axes[i,2].set_ylabel('Delta C component (10$^{-17}$ m)')
    for a in axes[i]:a.axhline(0,color='0.35',lw=.6);a.grid(alpha=.15)
for a in axes[0]:a.legend(fontsize=8)
for a in axes[-1]:a.set_xlabel('Time (h)')
fig.suptitle('B  Three-depth paired transport: changes remain at closure-residual scale')
fig.supxlabel('Gray bands are algebraic budget residuals, NOT time-step/grid uncertainty; those errors remain unknown.',fontsize=10)
export(fig,'Figure_B_transport')
fig,axes=plt.subplots(2,3,figsize=(13,7),layout='constrained')
vmax=max(v[2].max() for v in selected.values());alim=max(abs(v[3]).max() for v in selected.values())
for col,(index,(xs,zs,x,anomaly)) in enumerate(selected.items()):
    im=axes[0,col].pcolormesh(xs,zs,x,shading='nearest',vmin=.003598,vmax=vmax,cmap='viridis',rasterized=True)
    axes[0,col].set_title(f't = {t[index*50]/3600:.2f} h');axes[0,col].set_ylabel('z (m)' if col==0 else '')
    mask=zs>=zs[-1]-.005
    im2=axes[1,col].pcolormesh(xs,(zs[-1]-zs[mask])*1000,anomaly[mask]*1e6,shading='nearest',vmin=-alim*1e6,vmax=alim*1e6,cmap='RdBu_r',rasterized=True)
    axes[1,col].invert_yaxis();axes[1,col].set_xlabel('Horizontal position (m)');axes[1,col].set_ylabel('Depth below top (mm)' if col==0 else '')
fig.colorbar(im,ax=axes[0,:],label='Total x (mol/mol)',shrink=.8);fig.colorbar(im2,ax=axes[1,:],label='x - horizontal row mean (10$^{-6}$ mol/mol)',shrink=.8)
fig.suptitle('C  M1 concentration fields: full domain and explicitly magnified top 5 mm')
fig.supxlabel('Native nodal fields; shared color limits within each row. Bottom row removes the row mean and magnifies small perturbations.',fontsize=9)
export(fig,'Figure_C_fields')
fig,axes=plt.subplots(2,2,figsize=(11,8),layout='constrained')
axes[0,0].plot(hours,[r['M1_x_top_area_mean'] for r in base],color=blue,label='M1 surface mean');axes[0,0].plot(hours,[r['M1_x_domain_max'] for r in base],color=orange,ls='--',label='M1 domain maximum');axes[0,0].axhline(.1,color='black',ls=':',label='Solubility stop');axes[0,0].set(xlabel='Time (h)',ylabel='x (mol/mol)',ylim=(0,.105));axes[0,0].legend(fontsize=8)
axes[0,1].plot(hours,[(r['M1_x_top_max']-r['M0_x_top_max'])*1e6 for r in base],label='Maximum difference',color=orange);axes[0,1].plot(hours,[(r['M1_x_top_area_mean']-r['M0_x_top_area_mean'])*1e6 for r in base],label='Mean difference',color=blue,ls='--');axes[0,1].set(xlabel='Time (h)',ylabel='M1 - M0 (10$^{-6}$ mol/mol)');axes[0,1].legend(fontsize=8)
for key,color,ls in zip(REGIONS,[blue,orange,green],['-','--',':']):
 rr=regions[key];axes[1,0].plot([r['M1_x_top_max']/.1 for r in rr],[r['delta_C_over_reference']*1e12 for r in rr],label=f"z/H={rr[0]['actual_z_over_H']:.3f}",color=color,ls=ls)
axes[1,0].set(xlabel='M1 surface maximum / solubility threshold',ylabel='Delta C / initial inventory (10$^{-12}$)');axes[1,0].legend(fontsize=8)
for index,ls in [(0,':'),(4,'--'),(12,'-')]:
 _,zs,_,anomaly=selected[index];axes[1,1].plot(np.std(anomaly,axis=1)*1e6,(zs[-1]-zs)*1000,ls,label=f'{t[index*50]/3600:.1f} h')
axes[1,1].set(xlabel='Horizontal row SD (10$^{-6}$ mol/mol)',ylabel='Depth below top (mm)',ylim=(25,0));axes[1,1].axhline((1-regions[REGIONS[0]][0]['actual_z_over_H'])*200,color='0.5',ls='-.',label='Shallowest budget cut');axes[1,1].legend(fontsize=8)
for a in axes.flat:a.grid(alpha=.15)
fig.suptitle('D  Surface concentration, perturbation depth and unresolved extra transport')
fig.supxlabel('Differences near 10$^{-12}$ of inventory are not evidence of resolved transport. No physical effectiveness threshold assigned.',fontsize=9)
export(fig,'Figure_D_surface_transport')
dump(OUT/'figure_manifest.json',dict(destination='Provisional internal research report; journal/width requirements pending',matplotlib=matplotlib.__version__,numpy=np.__version__,source_hashes={str(p.relative_to(ROOT)):sha(p) for p in [OUT/'transport_timeseries.csv',OUT/'field_common_mode_timeseries.csv',OUT/'field_depth_profiles.csv',OUT/'new_M1_detection.json',Path(__file__)]},transformations=['No budget interpolation','Exact accepted times replace rounded PVD labels after output-step and rounding checks','Native nodal field display; row-mean anomaly explicitly disclosed','Common nx28 basis recomputed for both legs','No smoothing except frozen onset extractor five-point edge-padded boxcar'],uncertainty='Budget closure envelope only; sensitivity errors unknown',exports=exports))
print(json.dumps(dict(status=summary['status'],detection_s=detection['onset_refined_s'],regions={k:v['max_abs_delta_C_over_Iref'] for k,v in summary['regions'].items()},growth=summary['growth']),indent=2))
