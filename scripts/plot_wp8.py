import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from analyze_wp8 import OUT,CACHE,RUNS,rows
T=136423.0637723484
def finish(fig,name):
    fig.tight_layout();fig.savefig(OUT/(name+'.png'),dpi=170);fig.savefig(OUT/(name+'.pdf'));plt.close(fig)
def main():
    rr=rows(OUT/'mean_responses.csv');fig,ax=plt.subplots(1,2,figsize=(10,4))
    for g in ['R24','R48']:
        r=[x for x in rr if x['comparison']==g+'_A'];t=np.array([float(x['time_s']) for x in r])/86400
        ax[0].plot(t,[float(x['mean_response_top']) for x in r],label=g);ax[1].semilogy(t,[max(float(x['mean_response_Linf']),1e-18) for x in r],label=g)
    for a in ax:a.set_xlabel('Physical time (days)');a.legend();a.axvline(3*T/86400,c='.5',ls=':')
    ax[0].set(ylabel='M1 - M0 surface mean (mol/mol)',yscale='symlog');ax[0].set_yscale('symlog',linthresh=1e-10);ax[1].set_ylabel('Max-depth absolute mean response (mol/mol)');finish(fig,'Figure_A_mean_response')
    fig,ax=plt.subplots(1,2,figsize=(10,4))
    for g in ['R24','R48']:
        r=rows(OUT/(g+'_M1')/'timeseries.csv');t=np.array([float(x['time_s']) for x in r])/86400
        for a,key,label in zip(ax,['k_E_rad_m','k_peak_rad_m'],['Energy-weighted k (rad/m)','Peak k (rad/m)']):a.plot(t,[float(x[key]) if x[key] else np.nan for x in r],label=g);a.set(xlabel='Physical time (days)',ylabel=label);a.legend()
    finish(fig,'Figure_B_physical_wavenumbers')
    fig,ax=plt.subplots(1,2,figsize=(10,4),sharey=True)
    for a,g in zip(ax,['R24','R48']):
        d=np.load(CACHE/(g+'_M1')/'profiles_spectra.npz');en=d['energy'];v=np.log10(np.maximum(en/en.sum(axis=1)[:,None],1e-12));im=a.pcolormesh(d['time_s']/86400,d['k_rad_m'],v.T,vmin=-8,vmax=0,shading='auto',rasterized=True);a.set(title=g+' x 200 cells, M1',xlabel='Physical time (days)',ylabel='Physical k (rad/m)');fig.colorbar(im,ax=a,label='log10 energy fraction')
        if g=='R24':a.axhspan(d['k_rad_m'][-1],48*np.pi/(.6/28),color='.85',hatch='//',label='Not represented');a.legend(fontsize=7,loc='upper right')
        a.set_ylim(0,48*np.pi/(.6/28))
    finish(fig,'Figure_C_spectral_energy')
    tr=rows(OUT/'conservative_comparisons.csv.gz');fig,ax=plt.subplots(1,2,figsize=(10,4))
    for g in ['R24','R48']:
        for depth,region in [(1,'upper_0.995000'),(2,'upper_0.990000'),(5,'upper_0.975000')]:
            r=[x for x in tr if x['comparison']==g+'_A' and x['region']==region];t=np.array([float(x['time_s']) for x in r])/86400
            for a,key,label in zip(ax,['delta_C_over_initial','delta_I_over_current'],[r'$\Delta C / I_{M0}(0)$ (%)',r'$\Delta I_c / I_{M0}(t)$ (%)']):a.plot(t,[float(x[key])*100 for x in r],label=f'{g}, {depth} mm');a.set(xlabel='Physical time (days)',ylabel=label,yscale='symlog');a.set_yscale('symlog',linthresh=1e-5)
    ax[0].legend(fontsize=7);finish(fig,'Figure_D_conservative_transport')
    fig,ax=plt.subplots(1,2,figsize=(10,4));pr=rows(OUT/'selected_mean_profiles.csv')
    for a,g in zip(ax,['R24','R48']):
        r=[x for x in pr if x['comparison']==g+'_A'];tmax=max(float(x['time_s']) for x in r);r=[x for x in r if float(x['time_s'])==tmax and float(x['depth_mm'])<=10];a.plot([float(x['delta_mean']) for x in r],[float(x['depth_mm']) for x in r]);a.invert_yaxis();a.set(title=g+' M1-M0 at 4T',xlabel='Mean concentration response (mol/mol)',ylabel='Depth (mm)');a.ticklabel_format(axis='x',style='sci',scilimits=(0,0));a.xaxis.set_major_locator(plt.MaxNLocator(5))
    finish(fig,'Figure_E_depth_profiles')
    fig,ax=plt.subplots(1,2,figsize=(10,4))
    for g in ['R24','R48']:
        for a,typ in zip(ax,['B','C']):
            r=[x for x in rr if x['comparison']==g+'_'+typ];a.plot([float(x['time_s'])/86400 for x in r],[float(x['mean_response_top']) for x in r],label=g);a.set(xlabel='Physical time (days)',ylabel='Surface mean difference (mol/mol)',title=typ+': provisional one-cell reference');a.legend()
    fig.suptitle('Background available only to T; strict lateral-uniformity gate FAILED',fontsize=10);finish(fig,'Figure_F_early_background')
    fig,ax=plt.subplots(1,2,figsize=(10,4),sharey=True);vmax=0
    dd={g:dict(np.load(CACHE/(g+'_M1')/'profiles_spectra.npz')) for g in ['R24','R48']}
    for d in dd.values():vmax=max(vmax,float(abs(d['anomaly_4800']).max()))
    for a,(g,d) in zip(ax,dd.items()):
        im=a.pcolormesh(d['xs']*1000,(.2-d['z_m'])*1000,d['anomaly_4800'],shading='auto',cmap='RdBu_r',vmin=-vmax,vmax=vmax);a.set(xlabel='Horizontal position (mm)',ylabel='Depth (mm)',title=g+' x 200, M1 at 4T',ylim=(10,0));fig.colorbar(im,ax=a,label='x - own horizontal mean (mol/mol)')
    finish(fig,'Figure_G_own_mean_anomalies')
if __name__=='__main__':main()
