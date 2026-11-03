import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from wp8_strict_campaign import OUT,T
from analyze_wp8 import rows

def finish(fig,name):
    fig.tight_layout();fig.savefig(OUT/(name+'.png'),dpi=170);fig.savefig(OUT/(name+'.pdf'));plt.close(fig)
def main():
    rr=rows(OUT/'strict_responses.csv');colors={'R24':'#0072B2','R48':'#D55E00'};styles={'M0':'--','M1':'-'}
    fig,ax=plt.subplots(1,2,figsize=(10,4))
    for grid in colors:
        for role in styles:
            r=[v for v in rr if v['run']==grid+'_'+role];t=np.array([float(v['time_s']) for v in r])/86400
            for a,col in zip(ax,['delta_mean_top','delta_mean_Linf']):a.plot(t,[float(v[col]) for v in r],color=colors[grid],ls=styles[role],label=grid+' '+role);a.set_yscale('symlog',linthresh=1e-12);a.set_xlabel('Physical time (days)');a.axvline(3*T/86400,c='.7',ls=':',lw=.5)
    ax[0].set_ylabel('Surface mean minus strict background (mol/mol)');ax[1].set_ylabel('Max-depth absolute mean difference (mol/mol)');ax[0].legend(fontsize=8);finish(fig,'Figure_A_strict_mean')
    fig,ax=plt.subplots(1,3,figsize=(13.5,4))
    for grid in colors:
        for role in styles:
            r=[v for v in rr if v['run']==grid+'_'+role];t=np.array([float(v['time_s']) for v in r])/86400
            for a,col in zip(ax,['case_surface_mean','case_surface_max','delta_surface_max']):a.plot(t,[float(v[col]) for v in r],color=colors[grid],ls=styles[role],label=grid+' '+role)
    r=[v for v in rr if v['run']=='R24_M0'];t=np.array([float(v['time_s']) for v in r])/86400
    for a in ax[:2]:a.plot(t,[float(v['background_surface_mean']) for v in r],color='black',ls=':',label='Strict background');a.axhline(.1,color='.5',ls='--',label='Applicability threshold');a.set(xlabel='Physical time (days)',ylabel='Surface concentration (mol/mol)')
    ax[2].set_yscale('symlog',linthresh=1e-10);ax[2].set(xlabel='Physical time (days)',ylabel='Maximum minus background (mol/mol)',title='Local maximum difference')
    ax[0].set_title('Horizontal mean');ax[1].set_title('Local maximum');ax[0].legend(fontsize=7);finish(fig,'Figure_B_surface_concentrations')
    tr=rows(OUT/'strict_transport.csv.gz');fig,ax=plt.subplots(1,2,figsize=(10,4))
    for grid in colors:
        for depth,region in [(1.028871,'upper_0.995000'),(1.948008,'upper_0.990000'),(5.042057,'upper_0.975000')]:
            r=[v for v in tr if v['run']==grid+'_M1' and v['region']==region];t=np.array([float(v['time_s']) for v in r])/86400
            for a,col in zip(ax,['delta_I_over_current','delta_C_over_current']):a.plot(t,[100*float(v[col]) for v in r],label=f'{grid}, {depth:.3f} mm');a.set_yscale('symlog',linthresh=1e-6);a.set_xlabel('Physical time (days)')
    ax[0].set_ylabel(r'$\Delta I_c / I_b(t)$ (%)');ax[1].set_ylabel(r'$\Delta C / I_b(t)$ (%)');ax[0].legend(fontsize=7);finish(fig,'Figure_C_strict_inventory_transport')
    pr=rows(OUT/'strict_mean_profiles.csv');fig,ax=plt.subplots(1,2,figsize=(10,4),sharey=True)
    for a,grid in zip(ax,colors):
        for role in styles:
            r=[v for v in pr if v['run']==grid+'_'+role];t=max(float(v['time_s']) for v in r);r=[v for v in r if float(v['time_s'])==t and float(v['depth_mm'])<=10];a.plot([float(v['delta_mean']) for v in r],[float(v['depth_mm']) for v in r],ls=styles[role],label=role)
        a.set(title=f'{grid} at {t/T:.4f}T (last accepted)',xlabel='Mean minus strict background (mol/mol)',ylabel='Depth (mm)',ylim=(10,0));a.ticklabel_format(axis='x',style='sci',scilimits=(0,0));a.xaxis.set_major_locator(plt.MaxNLocator(5));a.legend()
    finish(fig,'Figure_D_depth_response')
if __name__=='__main__':main()
