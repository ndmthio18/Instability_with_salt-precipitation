"""Scientific figures from accepted states; no interpolation used for budgets."""
import csv,json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1];OUT=ROOT/'analysis/wp6_4T_exploratory_20261103';CACHE=ROOT/'tmp/wp6_4T_exploratory_20261103';T=136423.0637723484
def rows(p):
    with p.open() as f:return list(csv.DictReader(f))
def nums(r,k):return np.array([float(v[k]) for v in r])
def load(p):return json.loads(p.read_text())
def save(fig,n):fig.savefig(OUT/(n+'.png'),dpi=180);fig.savefig(OUT/(n+'.pdf'));plt.close(fig)
def main():
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    d=OUT/'runs/WP6-M1-explore4T';w=rows(d/'post/weighted_diagnostics.csv');ev=rows(d/'surface_events.csv');h=nums(ev,'time_s')/3600;sigma=nums(w,'weighted_top_sigma');sp=rows(OUT/'spatial/spatial_metrics.csv')
    fig,ax=plt.subplots(2,2,figsize=(11,7.5),layout='constrained')
    histories=[(ROOT/'analysis/wp6_execution_20261103/runs/WP6-M1-bridge','T reference (dt)',':'),(ROOT/'analysis/wp6_timestep_20261103/runs/WP6-M1-dt_half','2T reference (dt/2)','--'),(d,'4T exploratory (dt/2)','-')]
    for path,label,style in histories:
        r=rows(path/'post/weighted_diagnostics.csv');ee=rows(path/'surface_events.csv');tt=nums(ee,'time_s')/3600
        ax[0,0].semilogy(tt,nums(r,'weighted_top_sigma'),style,label=label,zorder=2 if path==d else 3,marker=None if path==d else ('o' if style==':' else 's'),markevery=240,markersize=4,markerfacecolor='none')
    ax[0,1].semilogy(h,sigma/min(sigma),label='Sigma / minimum');ax[0,1].semilogy(h,nums(w,'weighted_fourier')/min(nums(w,'weighted_fourier')),'--',label='Fourier / minimum')
    projection=abs(nums(w,'weighted_projection'));ax[0,1].semilogy(h,projection/min(projection),':',label='Eigenmode projection / minimum')
    ax[1,0].semilogy(h,nums(w,'weighted_top_relative'),label='Sigma / surface mean')
    fits=load(OUT/'growth_rate_fits.json');rr=[r for r in fits if r['quantity']=='weighted_top_sigma' and r.get('slope_per_s') is not None]
    ax[1,1].plot([(r['first_actual_s']+r['last_actual_s'])/7200 for r in rr],[r['slope_per_s'] for r in rr],'o-');ax[1,1].set_ylabel('Descriptive log-sigma slope (1/s)')
    ax[0,0].set_ylabel('Weighted top sigma (mol/mol)');ax[0,1].set_ylabel('Amplification from minimum');ax[1,0].set_ylabel('Relative surface perturbation (dimensionless)')
    for a in ax.flat:
        a.set_xlabel('Actual physical time (h)')
        for f in [1,2,4]:a.axvline(f*T/3600,color='.8',lw=.6)
    for a in [ax[0,0],ax[0,1],ax[1,0]]:a.legend(fontsize=8)
    fig.suptitle('A  Exploratory growth: T, 2T and 4T on physical time axes');save(fig,'Figure_A_growth')
    fields=np.load(CACHE/'spatial/selected_fields.npz');audit=load(OUT/'spatial/spatial_audit.json');selected=audit['selected']
    # Four actual snapshots: initial, nearest T, nearest 2T, terminal (may be censored).
    chosen=[selected[0]]+[min(selected,key=lambda r:abs(r['actual_time_s']-f*T)) for f in [1,2]]+[selected[-1]]
    fig,ax=plt.subplots(2,4,figsize=(14,6.8),layout='constrained');xx=fields['xs']*1000;depth=(.2-fields['zs'])*1000
    for j,m in enumerate(chosen):
        an=fields['anomaly_'+str(m['frame'])];peak=float(np.max(abs(an)));im=ax[0,j].pcolormesh(xx,depth,an,cmap='RdBu_r',vmin=-peak,vmax=peak,shading='nearest')
        fig.colorbar(im,ax=ax[0,j],shrink=.7,format='%.1e',label='x - horizontal mean')
        ax[0,j].set_title(f"t={m['actual_time_s']/3600:.2f} h\npeak |anomaly|={peak:.2e}")
        im2=ax[1,j].pcolormesh(xx,depth,an/peak,cmap='RdBu_r',vmin=-1,vmax=1,shading='nearest');ax[1,j].set_title('Normalized shape')
        for i in [0,1]:ax[i,j].set(xlabel='Horizontal x (mm)',ylabel='Depth (mm)',ylim=(10,0))
    fig.colorbar(im2,ax=list(ax[1,:]),shrink=.7,label='Normalized anomaly [-1,1]')
    fig.suptitle('B  Initial-shape control: top 10 mm; independent scales on top, common normalized scale below');save(fig,'Figure_B_spatial')
    fig,ax=plt.subplots(2,2,figsize=(11,7.5),layout='constrained')
    for key,label in [('fundamental_volume_energy_fraction','Fundamental m=2'),('even_higher_harmonic_fraction','Higher even harmonics'),('odd_neumann_energy_fraction','Odd Neumann modes')]:ax[0,0].semilogy(nums(sp,'time_s')/3600,np.maximum(nums(sp,key),1e-30),label=label)
    for key,label in [('depth90_mm','90% anomaly energy'),('depth99_mm','99% anomaly energy')]:ax[0,1].plot(nums(sp,'time_s')/3600,nums(sp,key),label=label)
    for m in chosen:
        profile=fields['sigma_'+str(m['frame'])];ax[1,0].semilogx(profile,depth,label=f"{m['actual_time_s']/3600:.1f} h")
        ax[1,1].plot(profile/profile[-1],depth,label=f"{m['actual_time_s']/3600:.1f} h")
    ax[0,0].set(xlabel='Actual physical time (h)',ylabel='Volume anomaly energy fraction',ylim=(1e-16,1.2));ax[0,1].set(xlabel='Actual physical time (h)',ylabel='Depth below surface (mm)')
    ax[1,0].set(xlabel='Horizontal sigma (mol/mol)',ylabel='Depth (mm)',ylim=(200,0));ax[1,1].set(xlabel='Sigma(z) / surface sigma',ylabel='Depth (mm)',ylim=(20,0))
    ax[1,0].text(.03,.56,'Tiny deep tails are shown,\nnot interpreted as fingers',transform=ax[1,0].transAxes,fontsize=8,bbox=dict(facecolor='white',alpha=.85,edgecolor='none'))
    for a in ax.flat:a.legend(fontsize=8)
    fig.suptitle('C  Modal evolution and penetration; inherited deep initial tails are retained');save(fig,'Figure_C_modes_depth')
    fig,ax=plt.subplots(2,2,figsize=(11,7.5),layout='constrained');pair=rows(OUT/'paired/paired_transport.csv');regions=['upper_0.995000','upper_0.990000','upper_0.975000','upper_0.950000','upper_0.900000','upper_0.750000','upper_0.500000'];labels=['1.03','1.95','5.04','9.87','20.10','49.54','99.89']
    for reg,label in zip(regions,labels):
        r=[r for r in pair if r['region']==reg];tt=nums(r,'time_s')/3600
        ax[0,0].plot(tt,nums(r,'delta_C_m'),label=label+' mm');ax[0,1].plot(tt,nums(r,'relative_delta_C_initial_inventory'),label=label+' mm')
    for role,path in [('M0',OUT/'runs/WP6-M0-explore4T'),('M1',d)]:
        e=rows(path/'surface_events.csv');tt=nums(e,'time_s')/3600
        ax[1,1].plot(tt,nums(e,'x_top_area_mean'),label=role+' surface mean');ax[1,1].plot(tt,nums(e,'x_domain_max'),'--',label=role+' domain max')
    ax[1,1].axhline(.1,color='k',ls=':',label='Solubility boundary')
    r=[r for r in pair if r['region']==regions[0]]
    for key,label in [('M1_signed_cut_adv','M1 adv'),('M1_signed_cut_diff','M1 diff'),('M1_signed_cut_total','M1 total'),('M0_signed_cut_total','M0 total')]:ax[1,0].plot(nums(r,'time_s')/3600,nums(r,key),label=label)
    ax[0,0].set_ylabel('Cumulative excess delta C (m)');ax[0,0].set_yscale('symlog',linthresh=1e-16);ax[0,1].set_ylabel('Delta C / initial M0 inventory');ax[0,1].set_yscale('symlog',linthresh=1e-9)
    ax[1,0].set_ylabel('1.03 mm signed downward flux (m/s)');ax[1,1].set_ylabel('Mole fraction')
    for a in ax.flat:a.set_xlabel('Actual physical time (h)');a.legend(fontsize=7);a.axhline(0,color='.8',lw=.5)
    fig.suptitle('D  Salt redistribution and solubility: quantitative transport remains unvalidated');save(fig,'Figure_D_transport_solubility')
    fig,ax=plt.subplots(2,2,figsize=(11,7.5),layout='constrained')
    for j,(role,folder) in enumerate([('M0','spatial_M0'),('M1','spatial')]):
        aa=load(OUT/folder/'spatial_audit.json');ff=np.load(CACHE/folder/'selected_fields.npz');frame=aa['final']['frame'];an=ff['anomaly_'+str(frame)];peak=np.max(abs(an));zz=(.2-ff['zs'])*1000
        im=ax[j,0].pcolormesh(ff['xs']*1000,zz,an,cmap='RdBu_r',vmin=-peak,vmax=peak,shading='nearest');ax[j,0].set(xlabel='Horizontal x (mm)',ylabel='Depth (mm)',ylim=(10,0),title=role+' final anomaly (own color scale)');fig.colorbar(im,ax=ax[j,0],label='mol/mol',format='%.1e')
        en=ff['mode_energy_'+str(frame)];ax[0,1].bar(np.arange(1,25)+(-.18 if role=='M0' else .18),en/sum(en),width=.36,label=role)
    ax[0,1].set(xlabel='Neumann mode index m (canonical m=2)',ylabel='Final volume energy fraction');ax[0,1].legend()
    w0=rows(OUT/'runs/WP6-M0-explore4T/post/weighted_diagnostics.csv');info=load(OUT/'paired/async_stopping_checks.json');idx=info['pair']['indices'];ratio=np.array([float(w0[i]['weighted_top_sigma'])/float(w[j]['weighted_top_sigma']) for i,j in idx]);tt=np.array(info['pair']['shared_times_s'])/3600
    ax[1,1].semilogy(tt,np.where(ratio>0,ratio,np.nan),label='M0 sigma / M1 sigma');ax[1,1].axhline(1e-3,color='k',ls='--',label='C01 engineering screen');ax[1,1].set(xlabel='Common accepted physical time (h)',ylabel='Background-to-perturbed sigma ratio',ylim=(1e-9,.2));ax[1,1].legend(fontsize=8);ax[1,1].text(.03,.92,'Initial ratio = 0 (not on log axis)',transform=ax[1,1].transAxes,fontsize=8)
    fig.suptitle('E  M0 loses uniformity; short-wavelength modes dominate the final states');save(fig,'Figure_E_background_breakdown')
if __name__=='__main__':main()
